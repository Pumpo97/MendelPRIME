"""Offline resource preparation; no model calls, downloads or inference changes."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import threading
from datetime import date, datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
if (HERE.parent / 'src').is_dir():
    sys.path.insert(0, str(HERE.parent / 'src'))


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def copy_file(source, destination):
    if not source.is_file():
        raise ValueError(f'Missing input file: {source.name}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def text_input(path):
    with path.open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt', encoding='utf-8-sig', newline='') if compressed else path.open(encoding='utf-8-sig', newline='')


def run_script(name, *arguments):
    subprocess.run([sys.executable, '-B', str(HERE / name), *map(str, arguments)], stdout=sys.stderr, check=True)


def assert_equal_hash(path, expected):
    if sha(path) != expected:
        raise ValueError(f'Unexpected snapshot content: {path.name}')


def is_bgzf(path):
    with path.open('rb') as stream:
        header = stream.read(12)
        if len(header) != 12 or header[:3] != b'\x1f\x8b\x08' or not header[3] & 4:
            return False
        extra = stream.read(int.from_bytes(header[10:12], 'little'))
    while len(extra) >= 4:
        length = int.from_bytes(extra[2:4], 'little')
        if extra[:2] == b'BC' and length == 2:
            return True
        extra = extra[4 + length:]
    return False


def bgzf_copy(source, destination):
    import pysam
    destination.parent.mkdir(parents=True, exist_ok=True)
    if is_bgzf(source):
        copy_file(source, destination)
        return
    with source.open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    with opener(source, 'rb') as incoming, pysam.BGZFile(str(destination), 'wb') as outgoing:
        shutil.copyfileobj(incoming, outgoing, 1024 * 1024)


def clinvar(args, payload, work):
    from mendelprime.knowledge.variants import _clinvar_info, new_variant_store
    output = payload / f'variant-knowledge/clinvar/clinvar_{args.release.replace("-", "")}.sqlite'
    run_script('build_clinvar_index.py', '--variant-summary', args.summary, '--var-citations', args.citations,
               '--output', output, '--release', args.release)
    state = new_variant_store({name: output for name in ('clinvar_db', 'clingen_validity', 'clingen_vcep', 'clingen_manifest')})
    if _clinvar_info(state)[1] is not None:
        raise ValueError('ClinVar output does not satisfy the runtime schema')
    with sqlite3.connect(output) as database:
        if database.execute('PRAGMA quick_check').fetchall() != [('ok',)]:
            raise ValueError('ClinVar integrity check failed')


def clinvar_protein(args, payload, work):
    with sqlite3.connect(args.database.resolve().as_uri() + '?mode=ro', uri=True) as database:
        release = database.execute("SELECT value FROM metadata WHERE key='clinvar_release'").fetchone()[0]
    if release != args.release:
        raise ValueError('Main database and requested protein release differ')
    output = payload / f'assessment/grch37/clinvar-protein/clinvar_protein_acmg_disease_{args.release.replace("-", "")}.sqlite'
    run_script('build_clinvar_protein_acmg_index.py', '--source', args.database, '--variant-summary', args.summary, '--output', output)
    with sqlite3.connect(output) as database:
        if not database.execute('SELECT count(*) FROM protein_variants').fetchone()[0]:
            raise ValueError('No qualifying protein comparator records')
        if database.execute('PRAGMA quick_check').fetchall() != [('ok',)]:
            raise ValueError('Protein database integrity check failed')


def revel(args, payload, work):
    import pysam
    from mendelprime.knowledge.annotation import _REVEL_COLUMNS
    if not shutil.which('sort'):
        raise ValueError('GNU sort is required for bounded-memory REVEL preparation')
    raw, sorted_rows = work / 'rows.tsv', work / 'sorted.tsv'
    count = 0
    with text_input(args.source) as stream, raw.open('w', encoding='utf-8', newline='') as output:
        reader = csv.DictReader(stream)
        if reader.fieldnames != list(_REVEL_COLUMNS):
            raise ValueError('Expected the nine-column REVEL CSV, including transcript IDs')
        for row in reader:
            values = [row.get(name) for name in _REVEL_COLUMNS]
            if None in row or any(value is None or any(c in value for c in '\t\r\n') for value in values):
                raise ValueError(f'Invalid REVEL row {reader.line_num}')
            if not values[0] or not values[1].isdigit() or int(values[1]) < 1:
                raise ValueError(f'Invalid hg19 coordinate at row {reader.line_num}')
            if not re.fullmatch('[ACGT]+', values[3]) or not re.fullmatch('[ACGT]+', values[4]):
                raise ValueError(f'Invalid allele at row {reader.line_num}')
            if values[7] not in {'', '.', 'NA'} and not 0 <= float(values[7]) <= 1:
                raise ValueError(f'Invalid REVEL score at row {reader.line_num}')
            output.write('\t'.join(values) + '\n')
            count += 1
    if not count:
        raise ValueError('Empty REVEL input')
    with sorted_rows.open('wb') as output:
        subprocess.run(['sort', '-s', '-t', '\t', '-k1,1', '-k2,2n', '--buffer-size=' + args.sort_memory,
                        '--temporary-directory=' + str(work), str(raw)], stdout=output,
                       env=dict(os.environ, LC_ALL='C'), check=True)
    destination = payload / 'assessment/grch37/predictions/new_tabbed_revel.tsv.gz'
    destination.parent.mkdir(parents=True)
    with pysam.BGZFile(str(destination), 'wb') as output, sorted_rows.open('rb') as source:
        output.write(('#' + '\t'.join(_REVEL_COLUMNS) + '\n').encode())
        shutil.copyfileobj(source, output, 1024 * 1024)
    pysam.tabix_index(str(destination), seq_col=0, start_col=1, end_col=1, zerobased=False)


def coverage(args, payload, work):
    import pysam
    with text_input(args.source) as source:
        header = source.readline().rstrip('\r\n')
        fields = header.lstrip('#').split('\t')
        if fields[:2] != ['chrom', 'pos'] or 'over_20' not in fields:
            raise ValueError('Expected chrom/pos/over_20 coverage header')
        column = fields.index('over_20')
        for number, line in enumerate(source, 2):
            values = line.rstrip('\r\n').split('\t')
            if len(values) != len(fields) or not values[1].isdigit() or int(values[1]) < 1:
                raise ValueError(f'Invalid coverage row {number}')
            if not 0 <= float(values[column]) <= 1:
                raise ValueError(f'Invalid over_20 value at row {number}')
    output = payload / 'assessment/grch37/gnomad/gnomad.exomes.coverage.summary.tsv.bgz'
    bgzf_copy(args.source, output)
    pysam.tabix_index(str(output), seq_col=0, start_col=1, end_col=1,
                     line_skip=0 if header.startswith('#') else 1, zerobased=False)


def vcf_index(args, payload, work):
    import pysam
    targets = {'gnomad-sites': 'gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz',
               'spliceai-snv': 'predictions/spliceai_scores.masked.snv.hg19.vcf.gz',
               'spliceai-indel': 'predictions/spliceai_scores.masked.indel.hg19.vcf.gz'}
    with text_input(args.source) as source:
        header = next((line.rstrip('\r\n') for line in source if line.startswith('#CHROM\t')), '')
    if header.split('\t') != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']:
        raise ValueError('Expected a sites-only eight-column VCF, not patient genotypes')
    output = payload / 'assessment/grch37' / targets[args.kind]
    bgzf_copy(args.source, output)
    pysam.tabix_index(str(output), preset='vcf')


def fasta(args, payload, work):
    import pysam
    target = 'reference/hg19.fa' if args.kind == 'reference' else 'vep/final_hg19.fa'
    output = payload / 'assessment/grch37' / target
    output.parent.mkdir(parents=True)
    with args.source.open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    if compressed:
        with gzip.open(args.source, 'rb') as source, output.open('wb') as destination:
            shutil.copyfileobj(source, destination, 1024 * 1024)
    else:
        copy_file(args.source, output)
    pysam.faidx(str(output))
    with pysam.FastaFile(str(output)) as reference:
        if not reference.references:
            raise ValueError('Empty FASTA')


def clingen(args, payload, work):
    from mendelprime.knowledge.variants import CSV_SPECS, _read_rows
    output = payload / 'variant-knowledge/clingen'
    output.mkdir(parents=True)
    copy_file(args.validity, output / 'gene_disease_validity.csv')
    with text_input(args.vcep) as source:
        first = source.readline()
    if '\t' in first:
        copy_file(args.vcep, output / 'variant_pathogenicity.tsv')
    else:
        with text_input(args.vcep) as source, (output / 'variant_pathogenicity.tsv').open('w', encoding='utf-8', newline='') as destination:
            writer = csv.writer(destination, delimiter='\t', lineterminator='\n')
            writer.writerows(csv.reader(source))
    urls = {'gene_disease_validity.csv': 'https://search.clinicalgenome.org/kb/gene-validity/download',
            'variant_pathogenicity.tsv': 'https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download?type=csv'}
    for _key, name, *formatting in CSV_SPECS.values():
        if not _read_rows(output / name, *formatting):
            raise ValueError('Empty ClinGen export: ' + name)
    write_json(output / 'resource_manifest.json', {'schema_version': 'mgagent_clingen_resource_v1',
               'release_date': args.release, 'source': 'ClinGen public downloads',
               'files': {name: {'url': url, 'sha256': sha(output / name)} for name, url in urls.items()}})


def aicura(args, payload, work):
    from mendelprime.knowledge.assessment_resources import _load_protocol
    template = HERE / 'assets/ai-cura'
    output = payload / 'guidelines/ai-cura'
    manifest = json.loads((template / 'resource_manifest.json').read_text())
    expected = manifest['primary_method']
    assert_equal_hash(args.workbook, expected['tables_workbook']['sha256'])
    workbook = output / 'adz4172_tables_s1_to_s13.xlsx'
    copy_file(args.workbook, workbook)
    knowledge = output / 'ai_cura_knowledgebase.json'
    run_script('extract_ai_cura_knowledgebases.py', '--workbook', workbook, '--output', knowledge)
    assert_equal_hash(knowledge, expected['derived_knowledgebase']['sha256'])
    for name in ('AI-CURA_Protocol.md', 'resource_manifest.json'):
        copy_file(template / name, output / name)
    _load_protocol({f'ai_cura_{key}': path for key, path in {
        'protocol': output / 'AI-CURA_Protocol.md', 'workbook': workbook,
        'knowledgebase': knowledge, 'manifest': output / 'resource_manifest.json'}.items()})


def svi(args, payload, work):
    from mendelprime.knowledge.assessment_resources import _load_clingen_guidance
    spec_path = HERE / 'assets/clingen-svi/manifest.json'
    spec = json.loads(spec_path.read_text())
    output = payload / 'guidelines/clingen-svi-20261021'
    output.mkdir(parents=True)
    for row in spec['resources']:
        source, destination = args.source_dir / row['name'], output / row['name']
        assert_equal_hash(source, row['sha256'])
        copy_file(source, destination)
        assert_equal_hash(destination, row['sha256'])
        text = output / row['text_file']
        if args.text_dir:
            copy_file(args.text_dir / row['text_file'], text)
        else:
            subprocess.run(['pdftotext', '-layout', str(destination), str(text)], check=True, timeout=300)
        assert_equal_hash(text, row['text_sha256'])
    copy_file(spec_path, output / 'manifest.json')
    _load_clingen_guidance({'clingen_guideline_manifest': output / 'manifest.json'})


def omim(args, payload, work):
    from mendelprime.knowledge.literature import load_omim
    output = payload / 'core/omim'
    output.mkdir(parents=True)
    seen = set()
    with text_input(args.source) as source, (output / 'omim_clean_full.jsonl').open('w', encoding='utf-8') as destination:
        for number, line in enumerate(source, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            identifier = str(row.get('omim_id', ''))
            if not re.fullmatch(r'\d{6}', identifier) or identifier in seen:
                raise ValueError(f'Invalid/duplicate OMIM ID at line {number}')
            sections = row.get('sections', {})
            if not isinstance(sections, dict):
                raise ValueError(f'Invalid sections at line {number}')
            paragraphs = sections.get('molecular_genetics', [])
            if not isinstance(paragraphs, list) or any(not isinstance(p, dict) or not isinstance(p.get('text'), str)
                    or not isinstance(p.get('pmids', []), list)
                    or any(not str(pmid).isdigit() for pmid in p.get('pmids', [])) for p in paragraphs):
                raise ValueError(f'Invalid paragraph/citation binding at line {number}')
            seen.add(identifier)
            destination.write(line if line.endswith('\n') else line + '\n')
    if not seen:
        raise ValueError('Empty OMIM mirror')
    copy_file(args.mim2gene, output / 'mim2gene.txt')
    store = load_omim(output / 'omim_clean_full.jsonl', output / 'mim2gene.txt')
    if len(store['records']) != len(seen) or not any(set(ids) & seen for ids in store['entries'].values()):
        raise ValueError('No matching gene/MIM records')


def genereviews(args, payload, work):
    columns = ['NBK_id', 'GR_shortname', 'genesymbol']
    with text_input(args.source) as source:
        reader = csv.DictReader(io.StringIO(source.read().lstrip('#')), delimiter='\t')
        if reader.fieldnames != columns:
            raise ValueError('Unexpected GeneReviews index header')
        rows = list(reader)
    if len(set(args.genes)) != len(args.genes):
        raise ValueError('Duplicate selected gene')
    selected = []
    for gene in args.genes:
        matching = [row for row in rows if row['genesymbol'] == gene]
        if not matching or any(not re.fullmatch(r'NBK\d+', row['NBK_id']) or not row['GR_shortname'] for row in matching):
            raise ValueError('Missing or invalid chapter mapping for ' + gene)
        selected.extend(matching)
    output = payload / 'report/genereviews/genereviews_index.tsv'
    output.parent.mkdir(parents=True)
    with output.open('w', encoding='utf-8', newline='') as destination:
        destination.write('#')
        writer = csv.DictWriter(destination, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(selected)


def autopvs1(args, payload, work):
    import pysam
    expected = json.loads((HERE / 'assets/autopvs1_source.json').read_text())
    output = payload / 'assessment/grch37/autopvs1'
    for name, digest in expected['files'].items():
        source = args.source_dir / name
        actual = sha(source)
        if actual != digest and not (name == 'read_data.py' and actual == expected['patched_read_data_sha256']):
            raise ValueError('Unexpected AutoPVS1 snapshot file: ' + name)
        copy_file(source, output / name)
        assert_equal_hash(output / name, actual)
    path = output / 'read_data.py'
    if sha(path) != expected['patched_read_data_sha256']:
        old = b"genome_hg19 = Fasta(config['HG19']['genome'])\ngenome_hg38 = Fasta(config['HG38']['genome'])"
        new = b"def optional_fasta(path):\n    return Fasta(path) if os.path.isfile(path) else None\n\n\ngenome_hg19 = optional_fasta(config['HG19']['genome'])\ngenome_hg38 = optional_fasta(config['HG38']['genome'])"
        content = path.read_bytes()
        if content.count(old) != 1:
            raise ValueError('Unexpected AutoPVS1 compatibility-patch location')
        path.write_bytes(content.replace(old, new))
    assert_equal_hash(path, expected['patched_read_data_sha256'])
    reference = output / 'data/hg19.fa'
    copy_file(args.reference, reference)
    pysam.faidx(str(reference))
    from mendelprime.tools.pvs1_driver import _CHILD_BOOTSTRAP
    subprocess.run([sys.executable, '-I', '-B', '-c', _CHILD_BOOTSTRAP, str(HERE.parent / 'src'),
                    '--check-source', str(output)], check=True, timeout=300)


def curated(args, payload, work):
    from mendelprime.knowledge.assessment_resources import _curated_documents
    from mendelprime.knowledge.reference import ReferenceGenome
    if not args.ack_reviewed:
        raise ValueError('Explicit --ack-reviewed is required; biological review is not automated')
    value = json.loads(args.catalog.read_text(encoding='utf-8'))
    if set(value) != {'schema_version', 'records'} or value['schema_version'] != 'verified_original_evidence_v1' or not value['records']:
        raise ValueError('Expected a nonempty reviewed evidence catalog')
    output = payload / 'report/curated-literature/reviewed'
    output.mkdir(parents=True)
    for row in value['records']:
        for name in ('identity', 'article'):
            item = row[name]
            for field in ('original', 'text'):
                relative = Path(item[field + '_file'])
                source = (args.catalog.parent / relative).resolve()
                if relative.is_absolute() or '..' in relative.parts or not source.is_relative_to(args.catalog.parent.resolve()):
                    raise ValueError('Catalog reference escapes its directory')
                if not source.is_file():
                    raise ValueError('Missing reviewed source file')
                digest = sha(source)
                if item.get(field + '_sha256') not in (None, digest):
                    raise ValueError('Reviewed source hash mismatch')
                item[field + '_sha256'] = digest
                copy_file(source, output / relative)
    write_json(output / 'catalog.json', value)
    with ReferenceGenome(args.reference.resolve(), assembly='GRCh37') as reference:
        state = {'settings': {'post_evidence_catalog': output / 'catalog.json'},
                 'lock': threading.RLock(), 'tables': {}, 'reference': reference}
        for gene in dict.fromkeys(row['gene'] for row in value['records']):
            targets = tuple({'primary': {'key': tuple(row['variant_key'])}} for row in value['records'] if row['gene'] == gene)
            if not _curated_documents(state, gene, targets):
                raise ValueError('Reviewed catalog produced no bound evidence')


def extract_text(args, payload, work):
    import xml.etree.ElementTree as ET
    output = payload / 'prepared-text' / (args.source.stem + '.txt')
    output.parent.mkdir(parents=True)
    if args.kind == 'pdf':
        subprocess.run(['pdftotext', '-layout', str(args.source), str(output)], check=True, timeout=300)
    else:
        root = ET.parse(args.source).getroot()
        output.write_text(ET.tostring(root, encoding='unicode'), encoding='utf-8')
    if not output.read_text(encoding='utf-8').strip():
        raise ValueError('No extracted text; scanned PDFs require independently reviewed OCR')


OPERATIONS = {'clinvar': clinvar, 'clinvar-protein': clinvar_protein, 'revel': revel, 'coverage': coverage,
              'vcf-index': vcf_index, 'fasta': fasta, 'clingen': clingen, 'aicura': aicura, 'svi': svi,
              'omim': omim, 'genereviews': genereviews, 'autopvs1': autopvs1, 'curated': curated,
              'extract-text': extract_text}


def parser():
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest='command', required=True)
    required = {'clinvar': ('summary', 'citations'), 'clinvar-protein': ('database', 'summary'),
                'clingen': ('validity', 'vcep'), 'aicura': ('workbook',), 'svi': ('source-dir',),
                'omim': ('source', 'mim2gene'), 'autopvs1': ('source-dir', 'reference'),
                'curated': ('catalog', 'reference')}
    for name in OPERATIONS:
        command = commands.add_parser(name)
        command.add_argument('--output-root', type=Path, required=True)
        for flag in required.get(name, ('source',)):
            command.add_argument('--' + flag, type=Path, required=True)
        if name in {'clinvar', 'clinvar-protein', 'clingen'}:
            command.add_argument('--release', required=True, help='Verified source snapshot date, YYYY-MM-DD')
        if name == 'revel':
            command.add_argument('--sort-memory', default='256M')
        if name == 'fasta':
            command.add_argument('--assembly', choices=['GRCh37'], required=True)
            command.add_argument('--kind', choices=['reference', 'vep'], default='reference')
        if name == 'vcf-index':
            command.add_argument('--kind', choices=['gnomad-sites', 'spliceai-snv', 'spliceai-indel'], required=True)
        if name == 'svi':
            command.add_argument('--text-dir', type=Path)
        if name == 'genereviews':
            command.add_argument('--genes', nargs='+', required=True)
        if name == 'curated':
            command.add_argument('--ack-reviewed', action='store_true', required=True)
        if name == 'extract-text':
            command.add_argument('--kind', choices=['pdf', 'xml'], required=True)
    return root


def primary_output(args):
    stamp = getattr(args, 'release', '').replace('-', '')
    paths = {
        'clinvar': f'variant-knowledge/clinvar/clinvar_{stamp}.sqlite',
        'clinvar-protein': f'assessment/grch37/clinvar-protein/clinvar_protein_acmg_disease_{stamp}.sqlite',
        'revel': 'assessment/grch37/predictions/new_tabbed_revel.tsv.gz',
        'coverage': 'assessment/grch37/gnomad/gnomad.exomes.coverage.summary.tsv.bgz',
        'clingen': 'variant-knowledge/clingen/resource_manifest.json',
        'aicura': 'guidelines/ai-cura/resource_manifest.json',
        'svi': 'guidelines/clingen-svi-20261021/manifest.json',
        'omim': 'core/omim/omim_clean_full.jsonl',
        'genereviews': 'report/genereviews/genereviews_index.tsv',
        'autopvs1': 'assessment/grch37/autopvs1/read_data.py',
        'curated': 'report/curated-literature/reviewed/catalog.json',
    }
    if args.command == 'fasta':
        return 'assessment/grch37/' + ('reference/hg19.fa' if args.kind == 'reference' else 'vep/final_hg19.fa')
    if args.command == 'vcf-index':
        return 'assessment/grch37/' + {'gnomad-sites': 'gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz',
            'spliceai-snv': 'predictions/spliceai_scores.masked.snv.hg19.vcf.gz',
            'spliceai-indel': 'predictions/spliceai_scores.masked.indel.hg19.vcf.gz'}[args.kind]
    if args.command == 'extract-text':
        return 'prepared-text/' + args.source.stem + '.txt'
    return paths[args.command]


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        for name, value in vars(args).items():
            if isinstance(value, Path) and name != 'output_root':
                valid = value.is_dir() if name.endswith('_dir') else value.is_file()
                if not valid:
                    raise ValueError('Missing or invalid input: ' + name)
        if hasattr(args, 'release') and date.fromisoformat(args.release).isoformat() != args.release:
            raise ValueError('Release must be an ISO date')
        if hasattr(args, 'sort_memory') and not re.fullmatch(r'[1-9][0-9]*[KMGT]?', args.sort_memory):
            raise ValueError('Invalid sort memory limit')
        root = args.output_root.resolve()
        root.mkdir(parents=True, exist_ok=True)
        primary = root / primary_output(args)
        if primary.exists() or primary.is_symlink() or not primary.resolve().is_relative_to(root):
            raise FileExistsError('Choose a fresh output location: ' + primary_output(args))
        print(f'preparation command={args.command} stage=input_checksums', file=sys.stderr, flush=True)
        inputs = {name: {'filename': value.name, 'bytes': value.stat().st_size, 'sha256': sha(value),
                         'mtime_ns': value.stat().st_mtime_ns}
                  for name, value in vars(args).items() if isinstance(value, Path) and name != 'output_root' and value.is_file()}
        with tempfile.TemporaryDirectory(prefix='.prepare-', dir=root) as temporary:
            stage = Path(temporary)
            payload, work = stage / 'payload', stage / 'work'
            payload.mkdir()
            work.mkdir()
            print(f'preparation command={args.command} stage=build', file=sys.stderr, flush=True)
            OPERATIONS[args.command](args, payload, work)
            for name, item in inputs.items():
                status = getattr(args, name).stat()
                if (status.st_size, status.st_mtime_ns) != (item['bytes'], item['mtime_ns']):
                    raise ValueError('Input changed during preparation: ' + name)
            files = sorted(path for path in payload.rglob('*') if path.is_file())
            if not files:
                raise ValueError('Preparation produced no files')
            receipt = {'command': args.command, 'prepared_at': datetime.now(timezone.utc).isoformat(),
                       'files': [{'path': str(path.relative_to(payload)), 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in files],
                       'declared_release': getattr(args, 'release', None), 'llm_calls': 0}
            receipt['inputs'] = inputs
            receipt_path = payload / 'preparation-receipts' / (stage.name.removeprefix('.') + '.json')
            write_json(receipt_path, receipt)
            files.append(receipt_path)
            print(f'preparation command={args.command} stage=publish files={len(receipt["files"])}', file=sys.stderr, flush=True)
            for source in files:
                destination = root / source.relative_to(payload)
                if not destination.resolve().is_relative_to(root) or destination.exists() or destination.is_symlink():
                    raise FileExistsError('Refusing to overwrite or escape output: ' + str(source.relative_to(payload)))
            published = []
            try:
                for source in files:
                    destination = root / source.relative_to(payload)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    os.link(source, destination)
                    published.append((destination, source.stat().st_ino))
            except Exception:
                for destination, inode in reversed(published):
                    if not destination.is_symlink() and destination.stat().st_ino == inode:
                        destination.unlink()
                raise
        print(json.dumps({'status': 'ok', 'command': args.command, 'files': len(receipt['files']),
                          'receipt': str(root / receipt_path.relative_to(payload))}, ensure_ascii=True))
        return 0
    except Exception as error:
        print(f'Preparation failed ({type(error).__name__}): {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
