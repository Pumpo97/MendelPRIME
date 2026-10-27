# Immutable preparation recipes

- ai-cura/ contains the already reviewed protocol and the exact manifest pinned
  by the scientific runtime. The original workbook remains a user-provided
  supplement; the knowledge base is reconstructed by the shipped extractor.
  Source: DOI 10.1126/scitranslmed.adz4172. These files are not a newly generated
  set of clinical rules. The public manifest uses basenames instead of machine
  paths; the reader supports this exact pinned revision and the original pin.
- clingen-svi/manifest.json identifies the three original guidance PDFs, URLs
  and expected extracted-text hashes. Originals are not embedded here.
- autopvs1_source.json records the verified upstream commit, tracked-file hashes
  and the existing compatibility-patched read_data.py hash. The installer obtains
  no source automatically; supply the upstream checkout and reference explicitly.
- genereviews_snapshot.tsv retains the current three gene/chapter mappings for
  optional report enrichment; it is not the full mapping or chapter text.

Keep these recipe bytes unchanged when reproducing the resource snapshot. Different source
versions require an explicit validated resource update, not a bypass of hashes.
Upstream material/tool attribution and conditions remain separate from the
MendelPRIME source license selection. No API credentials are stored here.
