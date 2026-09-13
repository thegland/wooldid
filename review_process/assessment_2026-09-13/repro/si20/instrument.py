"""Insert a READ-ONLY reference-cell probe before the two pooled joint-test margins calls
(binary att site = tag F2, continuous ame site = tag F1) of a wooldid.ado copy.
The probe is the site's own margins line minus contrast(overjoint), run under `cap qui`,
never posted; it prints REFM|<tag>|<colname>|<value> for the over(_t) column whose name
starts with 0b (the reference cell) and stores r(b)/r(V) in REPRO_ANCH_<tag>/REPRO_ANCHV_<tag>."""
import sys, re
src, dst, expect_f2, expect_f1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
lines = open(src, encoding='utf-8', errors='surrogateescape').read().split('\n')
def is_site(l, kind):
    return ("contrast(overjoint)" in l and "over(`_t')" in l and "[`weight'=`w_att']" in l
            and "subgroups" not in l and "`sgcase'" not in l
            and (("`semielasticity'(`targetflag')" in l) if kind == 'F2' else ("`contreatelasticitytype'(`contreat')" in l)))
sites = {}
for kind in ('F2', 'F1'):
    idx = [i for i, l in enumerate(lines) if is_site(l, kind)]
    assert idx, f'no {kind} site'
    sites[kind] = idx[0]
assert sites['F2'] + 1 == expect_f2, (sites['F2'] + 1, expect_f2)
assert sites['F1'] + 1 == expect_f1, (sites['F1'] + 1, expect_f1)
def probe(orig, tag):
    ind = orig[:len(orig) - len(orig.lstrip())]
    core = orig.lstrip()
    assert core.startswith('qui margins '), core[:40]
    pl = 'cap qui margins ' + core[len('qui margins '):].replace('contrast(overjoint)', '', 1)
    pl = re.sub(r'[ \t]{2,}', '  ', pl)
    body = f"""* ---- REPRO_SI20 probe (inserted by the stage-B proof agent into THIS COPY ONLY; not part of any build) ----
* the joint-test margins call below, minus contrast(overjoint): its level margins over `_t' include the
* reference cell (_t == 0). cap qui, never posted: e() is untouched and the real call follows unchanged.
cap matrix drop REPRO_ANCH_{tag}
cap matrix drop REPRO_ANCHV_{tag}
{pl}
if _rc == 0 {{
  tempname rpz_M rpz_V
  matrix `rpz_M' = r(b)
  matrix `rpz_V' = r(V)
  matrix REPRO_ANCH_{tag} = r(b)
  matrix REPRO_ANCHV_{tag} = r(V)
  local rpz_cn : colnames `rpz_M'
  local rpz_eq : coleq `rpz_M'
  local rpz_k 0
  foreach rpz_c of local rpz_cn {{
    local ++rpz_k
    local rpz_e : word `rpz_k' of `rpz_eq'
    if substr("`rpz_c'", 1, 2) == "0b" {{
      noi dis "REFM|{tag}|`rpz_e':`rpz_c'|" strtrim(string(el(`rpz_M', 1, `rpz_k'), "%21.15g"))
      noi dis "REFMV|{tag}|`rpz_e':`rpz_c'|" strtrim(string(el(`rpz_V', `rpz_k', `rpz_k'), "%21.15g"))
    }}
  }}
}}
else noi dis "REFM|{tag}|PROBE_RC|" _rc
* ---- end REPRO_SI20 probe ----"""
    return [ind + b for b in body.split('\n')]
out = []
for i, l in enumerate(lines):
    for kind, j in sites.items():
        if i == j:
            out.extend(probe(l, kind))
    out.append(l)
open(dst, 'w', encoding='utf-8', errors='surrogateescape').write('\n'.join(out))
print(f'{dst}: F2 site at source line {sites["F2"]+1}, F1 site at source line {sites["F1"]+1}; {len(out)-len(lines)} lines inserted')
