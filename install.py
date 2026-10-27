"""Install a selected capability; never run patients or call an LLM."""

import argparse
import hashlib
import json
import os
import platform
import shlex
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import tomllib
import urllib.request
import venv
import zipfile
from pathlib import Path, PurePosixPath

CODE = Path(__file__).resolve().parent
RELEASE = "https://github.com/Pumpo97/MendelPRIME/releases/download/1.0.0/"


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-type", choices=("hpo", "clinical", "vcf"), required=True)
    p.add_argument("--report", action="store_true")
    p.add_argument("--with-rarea", action="store_true", help="Optional HPO/clinical retrieval")
    p.add_argument("--work", type=Path, default=Path.home() / "mendelprime-work")
    p.add_argument("--resources", type=Path, help="Reuse or populate a resource root; defaults to WORK/resources")
    p.add_argument("--base-url", default=os.getenv("MENDELPRIME_BASE_URL"))
    p.add_argument("--model", default=os.getenv("MENDELPRIME_MODEL"))
    p.add_argument("--wire-api", choices=("chat_completions", "responses"))
    p.add_argument("--java", type=Path, help="Java 17 executable; otherwise discover it on PATH")
    p.add_argument("--ack-phenotype-upload", action="store_true", help="Authorize Phen2Gene and AMELIE uploads")
    p.add_argument("--archives", type=Path, help="Directory containing the selected Release archives")
    p.add_argument("--download-resources", action="store_true", help="Download only selected prepared archives")
    p.add_argument("--offline", action="store_true", help="Disable package indexes and resource downloads")
    p.add_argument("--wheelhouse", type=Path, help="Local dependency/build wheels for offline pip installation")
    modes = p.add_mutually_exclusive_group()
    modes.add_argument("--plan", action="store_true", help="Print the selected deployment without writing files")
    modes.add_argument("--configure-only", action="store_true", help="Generate/reuse TOML only; no installation or checks")
    return p


def selection(args):
    vcf = args.input_type == "vcf"
    if vcf and args.with_rarea:
        raise ValueError("RareA is available for HPO/clinical inputs, not VCF")
    if args.offline and args.download_resources:
        raise ValueError("--offline cannot be combined with --download-resources")
    profile = ("vcf" if vcf else "phenotype") + ("-report" if args.report else "-rank")
    if args.with_rarea:
        profile += "-rarea"
    ids = ["mendelprime-hpo-20260728-v1"]
    if vcf:
        ids += ["mendelprime-clinvar-20260117", "mendelprime-vcf-ranking-guidance-v1", "lirical-2.4.1"]
    if args.with_rarea:
        ids += ["mendelprime-similar-cases-rarea-20261027", "mendelprime-retrieval-support-20260815"]
    catalog = json.loads((CODE / "resources/ASSETS.json").read_text(encoding="utf-8"))["assets"]
    return profile, "vcf" if vcf else "hpo" if args.with_rarea else "", [catalog[key] for key in ids]


def toml_text(data):
    """Serialize the scalar/table/array-of-table subset used by the shipped profiles."""
    lines = []

    def table(values, prefix="", array=False):
        if prefix:
            lines.append(("[[" if array else "[") + prefix + ("]]" if array else "]"))
        children = []
        for key, value in values.items():
            name = f"{prefix}.{key}" if prefix else key
            if isinstance(value, dict) or isinstance(value, list) and value and isinstance(value[0], dict):
                children.append((name, value))
            else:
                lines.append(f"{key} = {json.dumps(value, ensure_ascii=False, allow_nan=False)}")
        lines.append("")
        for name, value in children:
            if isinstance(value, dict):
                table(value, name)
            else:
                for item in value:
                    table(item, name, True)

    table(data)
    output = "\n".join(lines)
    if tomllib.loads(output) != data:
        raise ValueError("Configuration could not be serialized without changes")
    return output


def configure(args, profile, work):
    destination = work / "configs" / profile
    run_path, paths_path = destination / "run.toml", destination / "resources.toml"
    if destination.exists():
        if not run_path.is_file() or not paths_path.is_file():
            raise ValueError(f"Incomplete existing configuration: {destination}")
        run = tomllib.loads(run_path.read_text(encoding="utf-8"))
        paths = tomllib.loads(paths_path.read_text(encoding="utf-8"))
        if run.get("resource_file") != "resources.toml":
            raise ValueError("Existing resource_file was customized; use main.py doctor directly")
        resources = (destination / paths["resource_root"]).resolve()
        requested = {"base_url": args.base_url, "model": args.model, "wire_api": args.wire_api}
        conflict = any(value is not None and run["llm"].get(key, "chat_completions" if key == "wire_api" else None) != value
                       for key, value in requested.items())
        conflict |= bool(args.resources and args.resources.resolve() != resources)
        conflict |= bool(args.java and paths.get("genotype", {}).get("java") != str(args.java.resolve()))
        conflict |= run["output"]["explanation"] != args.report or run["ranking"]["similar_cases"] != args.with_rarea
        conflict |= args.ack_phenotype_upload and not all(run.get("external", {}).get(k) for k in
                                                       ("phen2gene_acknowledged", "amelie_acknowledged"))
        if conflict:
            raise ValueError(f"Existing configuration is preserved. Edit {destination} or choose a new --work")
    else:
        name = "vcf_report" if args.report else "vcf_rank"
        if args.input_type != "vcf":
            name = "phenotype"
        run = tomllib.loads((CODE / f"configs/minimal/{name}.toml").read_text())
        paths = tomllib.loads((CODE / f"configs/minimal/{name}.resources.toml").read_text())
        resources = (args.resources or work / "resources").resolve()
        run["resource_file"] = "resources.toml"
        run["llm"].update(base_url=args.base_url, model=args.model, wire_api=args.wire_api or "chat_completions")
        run["output"]["explanation"] = args.report
        run["ranking"]["similar_cases"] = args.with_rarea
        paths["resource_root"] = str(resources)
        if args.input_type == "vcf":
            java = str(args.java.resolve()) if args.java else shutil.which("java")
            if java:
                paths["genotype"]["java"] = java
        else:
            for key in ("phen2gene_acknowledged", "amelie_acknowledged"):
                run["external"][key] = args.ack_phenotype_upload
        if args.with_rarea:
            expanded = tomllib.loads((CODE / "configs/hpo/resources.toml").read_text())["resources"]
            for key in ("similar_case_root", "embedding_model", "cross_encoder"):
                paths["resources"][key] = expanded[key]

    for key in ("base_url", "model"):
        value = run["llm"].get(key)
        if not isinstance(value, str) or not value.strip() or "<" in value or ".invalid" in value:
            raise ValueError(f"Supply --{key.replace('_', '-')} or edit llm.{key}; placeholders are not runnable")
    sys.path.insert(0, str(CODE / "src"))
    from mendelprime.runtime.config import load_config
    from mendelprime.domain import MendelPrimeError

    def validate(path):
        try:
            load_config(path)
        except MendelPrimeError as error:
            raise ValueError(f"Invalid configuration: {error.field}: {error.code}") from None

    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
            stage = Path(temporary) / profile
            stage.mkdir()
            (stage / "run.toml").write_text(toml_text(run), encoding="utf-8")
            (stage / "resources.toml").write_text(toml_text(paths), encoding="utf-8")
            validate(stage / "run.toml")
            stage.rename(destination)
    else:
        validate(run_path)
    return run_path, resources


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def archive_path(asset, directory, download):
    path = directory / asset["archive"]
    if not path.exists():
        if not download:
            raise FileNotFoundError(f"Missing archive: {path}; obtain it from {RELEASE}")
        directory.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {path.name} ({asset['archive_bytes']} bytes)", flush=True)
        with tempfile.TemporaryDirectory(dir=directory) as temporary:
            part = Path(temporary) / path.name
            with urllib.request.urlopen(RELEASE + path.name, timeout=90) as response, part.open("wb") as output:
                shutil.copyfileobj(response, output)
            if part.stat().st_size != asset["archive_bytes"] or sha(part) != asset["sha256"]:
                raise ValueError(f"Downloaded archive checksum mismatch: {path.name}")
            os.link(part, path)
    if not path.is_file() or path.stat().st_size != asset["archive_bytes"] or sha(path) != asset["sha256"]:
        raise ValueError(f"Archive checksum mismatch: {path.name}")
    return path


def safe_member(name):
    if PurePosixPath(name).is_absolute() or ".." in PurePosixPath(name).parts or "\\" in name or ":" in name:
        raise ValueError("Unsafe archive member")


def install_archive(archive, root):
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=root, prefix=".install-") as temporary:
        stage = Path(temporary)
        if archive.suffix == ".zip":
            with zipfile.ZipFile(archive) as package:
                for member in package.infolist():
                    safe_member(member.filename)
                    if stat.S_ISLNK(member.external_attr >> 16):
                        raise ValueError("Archive symlinks are not allowed")
                package.extractall(stage)
            payload, target = stage / "lirical-cli-2.4.1", root / "genotype/lirical"
        else:
            with tarfile.open(archive) as package:
                for member in package:
                    safe_member(member.name)
                    if not member.isfile() and not member.isdir():
                        raise ValueError("Archive links and special files are not allowed")
                package.extractall(stage, filter="data")
            payload, target = stage, root
        files = [(p, target / p.relative_to(payload)) for p in payload.rglob("*") if p.is_file()]
        if not files:
            raise ValueError("Archive contains no installable files")
        for source, destination in files:
            if not destination.resolve().is_relative_to(root.resolve()) or destination.is_symlink():
                raise ValueError("Resource destination escapes its root or is a symlink")
            if destination.exists() and (not destination.is_file() or sha(destination) != sha(source)):
                raise FileExistsError(f"Different resource exists; refusing to overwrite: {destination}")
        for source, destination in files:
            if not destination.exists():
                destination.parent.mkdir(parents=True, exist_ok=True)
                os.link(source, destination)
    print(f"Installed/verified {archive.name}", flush=True)


def command(argv, *, offline=False):
    print("+ " + shlex.join(map(str, argv)), flush=True)
    env = None
    if offline:
        env = {k: v for k, v in os.environ.items() if not k.startswith("PIP_")}
        env["PIP_CONFIG_FILE"] = os.devnull
    subprocess.run(list(map(str, argv)), check=True, cwd=CODE, env=env)


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if sys.version_info[:2] != (3, 12):
            raise ValueError("Use Python 3.12; see docs/INSTALL.md for the Conda alternative")
        profile, extra, assets = selection(args)
        if args.resources:
            args.resources = args.resources.expanduser()
        if args.java:
            args.java = args.java.expanduser()
        if args.archives:
            args.archives = args.archives.expanduser().resolve()
            if args.archives.is_relative_to(CODE):
                raise ValueError("Keep downloaded archives outside the source repository")
        work = args.work.expanduser().resolve()
        resources = (args.resources or work / "resources").expanduser().resolve()
        if work.is_relative_to(CODE) or resources.is_relative_to(CODE):
            raise ValueError("Keep work and installed resources outside the source repository")
        if args.plan:
            print(json.dumps({"profile": profile, "python_extra": extra or "base", "work": str(work),
                              "resources": str(resources), "prepared_assets": [a["archive"] for a in assets],
                              "resource_download_bytes": sum(a["archive_bytes"] for a in assets),
                              "additional_resources": "See docs/INSTALL.md#resource-to-parameter-map"}, indent=2))
            return 0
        if not args.configure_only and (platform.system() != "Linux" or platform.machine() not in {"x86_64", "amd64"}):
            raise ValueError("Installation requires Linux x86-64; --plan/--configure-only are portable")
        os.umask(0o077)
        config, resources = configure(args, profile, work)
        if resources.is_relative_to(CODE):
            raise ValueError("Existing resource_root points into the source repository")
        print(f"Configuration: {config}\nResources: {resources}", flush=True)
        if args.configure_only:
            print("Configuration only: environment, resources and services have NOT been checked.")
            return 0
        env = work / ".venv"
        python = env / "bin/python"
        if not env.exists():
            venv.EnvBuilder(with_pip=True).create(env)
        if not python.is_file():
            raise ValueError("Existing WORK/.venv is incomplete; use a new work directory")
        command([python, "-c", "import sys; assert sys.version_info[:2] == (3, 12)"])
        pip = [python, "-m", "pip", "install", "-c", CODE / "constraints/runtime.txt"]
        if args.offline:
            pip += ["--no-index"]
        if args.wheelhouse:
            pip += ["--find-links", args.wheelhouse.resolve()]
        with tempfile.TemporaryDirectory(dir=work, prefix=".build-") as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            shutil.copy2(CODE / "pyproject.toml", source / "pyproject.toml")
            shutil.copytree(CODE / "src", source / "src",
                            ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"))
            command([*pip, str(source) + (f"[{extra}]" if extra else "")], offline=args.offline)
        command([python, "-m", "pip", "check"])
        if args.archives or args.download_resources:
            directory = (args.archives or work / "downloads").resolve()
            for asset in assets:
                install_archive(archive_path(asset, directory, args.download_resources), resources)
        for name in ("cache/doctor", "outputs"):
            (work / name).mkdir(parents=True, exist_ok=True)
        print("Checking configured resources; missing tools/data must be prepared using docs/INSTALL.md.", flush=True)
        command([python, CODE / "main.py", "doctor", "--config", config,
                 "--profile", "vcf" if args.input_type == "vcf" else "hpo", "--check-tools",
                 "--work-root", work / "cache/doctor"])
        print("Local deployment checks passed. API availability and scientific performance were NOT tested.")
        return 0
    except (OSError, ValueError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        print(f"Installation incomplete: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
