#!/usr/bin/env python3
"""Generate Stata shard do-files from suite config CSVs.

Usage:
  mkrun.py --configs configs.csv [--configs configs_ext.csv ...] --preps preps.json
           --ado /abs/path/wooldid.ado --label old|new --sha <gitsha>
           --data /abs/suite/data [--data /abs/test/data ...]
           --outdir /abs/rundir [--dumpdir /abs/dumpdir]
           [--components core,trig,...] [--shard-size 15]

Emits outdir/shard_<component>[_slow|_vslow]_<n>.do plus outdir/shards.manifest (one
line per shard: dofile|nconfigs|tier, tier = fast|slow|vslow). The configs' `slow' column
routes: '' -> fast shards (600 s), '1' -> *_slow shards (1800 s), '2' -> *_vslow shards of
ONE config each (3600 s). Each shard dumps one
<config>.csv per config, in the battery %21x format plus a "prov|<label>|<sha>"
provenance line.

--dumpdir sends those dumps to a fixed absolute directory instead of the runner's cwd.
That is what lets the orchestrator give every shard its OWN cwd (SI-21): the saveplots
capture below globs the whole working directory and then ERASES what it listed, so
shards sharing a cwd could both steal and destroy each other's saveplots evidence. With
per-shard cwds the glob sees only that shard's own output, while the dumps still land
together where the comparer expects them.
"""
import sys, os, csv, json
import re as _re2

# the file-writing options and their legal abbreviations: SAVEplots() -> save..saveplots,
# PERMSAve() -> permsa..permsave, BOOTSAve() -> bootsa..bootsave. Used to pin every file a
# config can write (documentation review D-11).
_SAVEOPT = _re2.compile(
    r"\b(?:save(?:p(?:l(?:o(?:t(?:s)?)?)?)?)?"
    r"|(?:perm|boot)sa(?:v(?:e)?)?)"
    r"\s*\(([^)]*)\)"
)

DUMPRES = r'''
cap program drop dumpres
program define dumpres
  gettoken cfg 0 : 0
  gettoken rc  0 : 0
  * every remaining token is a file prefix to pin (see the files| note below)
  local svpre `0'
  tempname fh
  file open `fh' using `"{DUMPDIR}`cfg'.csv"', write replace
  file write `fh' "rc|`rc'" _n
  file write `fh' "prov|{LABEL}|{SHA}" _n
  qui ds
  file write `fh' `"varlist|`r(varlist)'"' _n
  * e(sample) digest: N + per-variable sums over the realized sample (order-independent,
  * so internal re-sorts cannot fake a match; a sample-composition change shows even
  * when b/V coincidentally agree)
  cap qui count if e(sample)
  if _rc == 0 {
    local v : di %21x r(N)
    file write `fh' "esamp|__N|`v'" _n
    foreach var of varlist _all {
      cap qui sum `var' if e(sample), meanonly
      if _rc == 0 {
        local v : di %21x r(sum)
        file write `fh' "esamp|`var'|`v'" _n
      }
    }
  }
  local escalars : e(scalars)
  foreach s of local escalars {
    local v : di %21x e(`s')
    file write `fh' "scalar|`s'|`v'" _n
  }
  local emacros : e(macros)
  foreach m of local emacros {
    file write `fh' `"macro|`m'|`e(`m')'"' _n
  }
  local emats : e(matrices)
  foreach m of local emats {
    tempname M
    cap matrix `M' = e(`m')
    local mrc = _rc
    if `mrc' != 0 {
      * an e() matrix the dump cannot read used to vanish silently (harness audit H-11,
      * 2026-09-12); record it so a one-sided unreadable matrix is a matnames diff. Two
      * fields after the name on purpose - never the RxC|rn|cn shape the comparer parses.
      file write `fh' `"matnames|`m'|UNREADABLE|rc=`mrc'"' _n
      continue
    }
    local R = rowsof(`M')
    local C = colsof(`M')
    local rn : rowfullnames `M'
    local cn : colfullnames `M'
    file write `fh' `"matnames|`m'|`R'x`C'|`rn'|`cn'"' _n
    forval r = 1/`R' {
      forval c = 1/`C' {
        local v : di %21x `M'[`r',`c']
        file write `fh' "mat|`m'|`r'|`c'|`v'" _n
      }
    }
    cap matrix drop `M'
  }
  * SI-21: this globs the WHOLE working directory and erases what it lists, so it is
  * only correct while each shard owns its cwd. run_full.sh guarantees that by giving
  * every shard a fresh private workdir and sending dumps to --dumpdir.
  * The pattern is the config's OWN saveplots() prefix when it has one, falling back to
  * "svp" as a catch-all for strays. The old hardcoded "svp*" was simultaneously too
  * broad (it swept other shards' files) and too narrow (a config saving under any other
  * prefix - e.g. saveplots(wdbt_saved) - was never checked at all).
  * D-11 (author ruling 2026-09-11): `svpre' is now a LIST of prefixes, not one. It used to
  * be the single saveplots() argument matched on that option's literal spelling, so the
  * legal abbreviation save(...) slipped the pin entirely, and permsave()/bootsave() - which
  * also write files - were never pinned at all.
  if "`svpre'" == "" local svpre svp
  local svfiles
  foreach sp of local svpre {
    local sp1 : dir . files "`sp'*"
    * `: dir' returns a QUOTED list; concatenating it with a plain `local x `x' `list''
    * re-parses the line and strips the outer quotes, so one name keeps a stray quote,
    * its erase fails (rc 198) and the files| record is mangled (D-11 defect, 2026-09-12).
    * foreach unquotes each element, so append them one at a time.
    foreach svf of local sp1 {
      local svfiles `svfiles' `svf'
    }
  }
  local svfiles : list uniq svfiles
  local svfiles : list sort svfiles
  file write `fh' `"files|`svfiles'"' _n
  foreach svf of local svfiles {
    cap erase `svf'
  }
  file write `fh' "end|`cfg'" _n
  file close `fh'
end
'''

HEADER = '''* Auto-generated by mkrun.py — do not edit. Build: {LABEL} ({SHA})
clear all
set more off
set maxvar 20000
set graphics off
* CAUTION (harness-side half of a trap the ado documents from its own side): loading with
* `do <path>' is the ONE load style under which wooldid cannot recover its own Mata kernels
* after a `clear mata'; an adopath-resolved ado recovers because Stata re-reads the file when
* a function is missing, and a do-loaded one does not (SI-17/SI-18). The recovery path in the
* ado is a `cap mata:' that EXITS on failure, so the kernel loss is swallowed rather than
* raised - the run keeps the older analytic-diagonal e(V). Because this header runs once per
* SHARD and each shard then runs up to --shard-size configs in the same Stata process, a
* config or prep that clears Mata would leave every later config in that shard silently on
* the fallback path, and the per-config scrub block below cannot undo it. None does today
* (the cleanmatrices option that did was retired 2026-09-25); isolate one that ever does, or
* load by adopath.
do "{ADO}"
{DUMPRES}
'''

def main():
    args = sys.argv[1:]
    def allvals(flag):
        return [args[i + 1] for i, a in enumerate(args) if a == flag]
    def optval(flag, default=None):
        v = allvals(flag)
        return v[0] if v else default
    cfg_files = allvals("--configs")
    preps = json.load(open(optval("--preps"))) if optval("--preps") else {}
    ado = os.path.abspath(optval("--ado"))
    label = optval("--label")
    sha = optval("--sha", "unknown")
    import re as _re
    if not _re.fullmatch(r"[0-9a-f]{7,40}(-dirty)?", sha):
        print(f"WARNING: --sha '{sha}' is not a git sha — dumps will fail the "
              f"comparer's provenance conventions in orchestrated runs")
    data_dirs = [os.path.abspath(d) for d in allvals("--data")]
    outdir = os.path.abspath(optval("--outdir"))
    comps = set(optval("--components", "").split(",")) if optval("--components") else None
    shard_size = int(optval("--shard-size", "15"))
    os.makedirs(outdir, exist_ok=True)

    rows = []
    for f in cfg_files:
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                if not r.get("config") or r["config"].startswith("#"):
                    continue
                if comps and r["component"] not in comps:
                    continue
                rows.append(r)
    seen = set()
    for r in rows:
        if r["config"] in seen:
            sys.exit(f"duplicate config name: {r['config']}")
        seen.add(r["config"])

    skip_missing = "--skip-missing-data" in args

    def resolve(ds):
        if ds.endswith(".dta"):
            ds = ds[:-4]
        for d in data_dirs:
            p = os.path.join(d, ds + ".dta")
            if os.path.exists(p):
                return p
        if skip_missing:
            return None
        sys.exit(f"dataset not found in any data dir: {ds}")

    def block(r):
        lines = [f"* ---- {r['config']} ({r['component']}) ----"]
        # BOTH SEEDS ARE RESET PER CONFIG, not once per shard, and that is the point:
        # `set seed' pins the stream bootstrap()/permute() draw from, `set sortseed' pins
        # the separate seed Stata uses to break ties in `sort'. The ado never touches
        # sortseed - this is a harness decision, so that two builds executing different
        # numbers of internal sorts still agree. Resetting per config decouples a config's
        # draws and tie order from its POSITION in the shard, which is what makes
        # --shard-size and CSV row order safe to change. The project has measured that
        # session-order effects are real (see test/README.md on cont_lat / cont_uz moving
        # at 1e-16..1e-10 in a shared session but bit-identical when run alone).
        lines.append("set seed 98765")
        lines.append("set sortseed 98765")
        lines.append("ereturn clear")
        # scrub session state a previous config in this shard may have left:
        # wooldid's own globals (new build only — build-asymmetric by construction),
        # user matrices/scalars from preps, and stored estimates
        lines.append("cap macro drop WOOLDID_LASTRC WOOLDID_CLEANUP WOOLDID_LEVEL0")
        lines.append("cap matrix drop _all")
        lines.append("cap scalar drop _all")
        lines.append("cap estimates drop _all")
        if r.get("dataset", "").strip():
            p = resolve(r["dataset"].strip())
            if p is None:
                print(f"SKIPPED (missing dataset {r['dataset']}): {r['config']}")
                return None
            lines.append(f'use "{p}", clear')
        else:
            lines.append("clear")
        prep = r.get("prep", "").strip()
        if prep:
            if prep not in preps:
                sys.exit(f"unknown prep key: {prep} (config {r['config']})")
            lines.extend(preps[prep])
        call = "cap noi wooldid " + r.get("varlist", "y i t ttre").strip()
        if r.get("weight", "").strip():
            call += f" [{r['weight'].strip()}]"
        if r.get("ifcond", "").strip():
            call += f" if {r['ifcond'].strip()}"
        call += ", " + r["options"].strip()
        lines.append(call)
        lines.append("local rc = _rc")
        # EVERY FILE-WRITING OPTION, AT EVERY LEGAL ABBREVIATION (author ruling 2026-09-11,
        # documentation review D-11). This used to be a literal match on `saveplots(', so
        # save(...) - the option's declared minimum abbreviation, which resolves uniquely -
        # slipped the pin silently and left the run pinning only the "svp" fallback, the
        # pre-SI-21 state. permsave() and bootsave() also write files and were outside the
        # pin altogether. The three declarations are SAVEplots(), PERMSAve() and BOOTSAve(),
        # which fixes the minimum abbreviations below; keep this in step with them.
        svpres = [m.group(1).strip() for m in _re2.finditer(_SAVEOPT, r["options"])]
        svpres = [p for p in dict.fromkeys(svpres) if p]
        lines.append(f"dumpres {r['config']} `rc' {' '.join(svpres)}".rstrip())
        lines.append("")
        return "\n".join(lines)

    # shard: group by (component, tier, isolated), chunk by shard_size.
    # tier: slow column '' -> fast (600 s), '1' -> slow (1800 s), '2' -> vslow (3600 s,
    # one config per shard - for the v1.3-side outliers such as si20_lat_pois at ~950 s).
    # isolated: one config per shard. The cleanmatrices configs were isolated here too until
    # that option was retired 2026-09-25 (see the HEADER caution).
    groups = {}
    for r in rows:
        sl = r.get("slow", "").strip()
        tier = "vslow" if sl == "2" else "slow" if sl == "1" else "fast"
        iso = sl == "2"
        key = (r["component"], tier, iso)
        groups.setdefault(key, []).append(r)
    manifest = []
    dumpdir = optval("--dumpdir")
    if dumpdir:
        dumpdir = os.path.abspath(dumpdir)
        os.makedirs(dumpdir, exist_ok=True)
        dumpdir = dumpdir.rstrip("/") + "/"
    else:
        dumpdir = ""
    header = (HEADER.replace("{ADO}", ado)
                    .replace("{LABEL}", label)
                    .replace("{SHA}", sha)
                    .replace("{DUMPRES}", DUMPRES.replace("{LABEL}", label)
                                                 .replace("{SHA}", sha)
                                                 .replace("{DUMPDIR}", dumpdir)))
    for (comp, tier, iso), rs in sorted(groups.items()):
        csize = 1 if iso else shard_size
        sfx = {"fast": "", "slow": "_slow", "vslow": "_vslow"}[tier]
        for i in range(0, len(rs), csize):
            chunk = rs[i:i + csize]
            blocks = [b for b in (block(r) for r in chunk) if b is not None]
            if not blocks:
                continue
            name = f"shard_{comp}{sfx}_{i // csize:02d}.do"
            with open(os.path.join(outdir, name), "w") as f:
                f.write(header)
                for b in blocks:
                    f.write(b + "\n")
                f.write('di "SHARD_DONE"\n')
                f.write("exit, clear\n")
            manifest.append(f"{name}|{len(blocks)}|{tier}")
    with open(os.path.join(outdir, "shards.manifest"), "w") as f:
        f.write("\n".join(manifest) + "\n")
    print(f"wrote {len(manifest)} shards for {len(rows)} configs to {outdir}")

main()
