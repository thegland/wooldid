# Help editor tooling

The artifact https://claude.ai/artifact/1MuFdJQ4pfTdhtE2UtbXBm renders wooldidinstall/wooldid.sthlp block
by block with every paragraph, heading and table cell editable in place. Edits are saved to the artifact's
database, collection `edits`, one document `b<idx>` per edited block: {ver, idx, line, kind, hash,
smcl_original, smcl, edited_at, by}. `ver` is the short git sha of the sthlp the page was built from.

- parse_sthlp.py <sthlp> <blocks.json> --ver TAG: split the help into blocks; the concatenation of every
  block's `smcl` reproduces the file byte for byte (asserted).
- build.py <sthlp> <ver> <out.html>: run the parser and inject the block list into the template.
- apply_edits.py <edits_dir> <sthlp> --ver TAG [--write]: read the record's documents (saved as JSON files
  by ArtifactData with out_dir; use a directory INSIDE the project, an outside path prompts for approval),
  match each edit to a block by hash, re-wrap the edited block at 170 columns at spaces outside braces,
  print the diff, and with --write replace the file.
- blocks_<ver>.json: the block list the published page carries, kept so edits can be applied after the
  scratch directory is gone.

To publish a new version: build from the current sthlp with its sha as the tag, then republish to the
artifact URL. Edits stored under an older `ver` are shown on the page as stale (the text is kept for
copying) and apply_edits.py skips them unless the block's hash still matches.
