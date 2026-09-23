#!/usr/bin/env python3
"""Generate the guide's messages chapter from the ado. Usage: gen_messages_chapter.py <ado> <out.md>"""
import re, sys
ado=open(sys.argv[1],encoding="utf-8").read().split("\n")
prog=[None]*len(ado); cur=None
for i,l in enumerate(ado):
    m=re.match(r'\s*prog(?:ram)?\s+def(?:ine)?\s+([A-Za-z_0-9]+)', l)
    if m: cur=m.group(1)
    prog[i]=cur
def _braces(line):
    # brace balance of a line with string literals removed, so braces in messages do not count
    t=re.sub(r'"[^"]*"', '', line)
    return t.count('{') - t.count('}')
def _else_of(i):
    # the condition an inline "else" at line i negates: an inline if on the line above, or the
    # if-block that closes just above
    m1=re.match(r'\s*(?:else )?if (.*?)\s+(?:dis|display|local|c_local)\b', ado[i-1])
    if m1 and not ado[i-1].rstrip().endswith('{'): return "not (" + m1.group(1).strip() + ")"
    bal=0
    for j in range(i-1, max(-1,i-300), -1):
        bal+=_braces(ado[j])
        if bal==0:
            m2=re.match(r'\s*(?:\}\s*)?(?:else\s+)?if (.*?)\s*\{\s*$', ado[j])
            if m2: return "not (" + m2.group(1).strip() + ")"
        if bal>0: break
    return "(else branch of the preceding condition)"
def cond_for(i, l):
    if re.match(r'\s*else\s+(?:dis|display|local|c_local)\b', l): return _else_of(i)
    mm=re.match(r'\s*(?:else )?if (.*?)\s*\{?\s*(?:dis|display|local|c_local)\b', l)
    if mm: return mm.group(1).strip()
    # walk back to the block that ENCLOSES line i: a closed sibling block above it (an inner
    # if that errors out, say) is skipped because its braces balance out
    bal=0
    for j in range(i-1, max(-1,i-60), -1):
        # "} else {" or "} else if c {" on one line nets zero braces; when it is reached at
        # balance 0 the message sits directly in that branch
        if bal==0:
            me=re.match(r'\s*\}\s*else(?:\s+if\s+(.*?))?\s*\{\s*$', ado[j])
            if me: return me.group(1).strip() if me.group(1) else "(else branch of the preceding condition)"
        bal+=_braces(ado[j])
        if bal>0:
            m2=re.match(r'\s*(?:else )?if (.*?)\s*\{\s*$', ado[j])
            if m2: return m2.group(1).strip()
            if re.match(r'\s*(?:\}\s*)?else\s*\{?\s*$', ado[j]): return "(else branch of the preceding condition)"
            return ""
    return ""
def body(l):
    s=l.strip()
    m=re.search(r'"((?:Caution|Note)[:\s].*)"\s*$', s)   # to the LAST quote: keeps embedded expression pieces
    if m: return m.group(1)
    m=re.search(r"local \S+ (.*)$", s); return m.group(1) if m else s
caut=[]; errs=[]
for i,l in enumerate(ado):
    if re.search(r'\bdi(s|splay)?\b', l) and re.search(r'"(Caution|Note)[:\s]', l):
        caut.append((prog[i], cond_for(i,l), body(l)))
    elif re.match(r'\s*local (wdbc`wdbncaut\'|wdpc`ncaut\'|imdegmsg)\s', l):
        caut.append((prog[i], cond_for(i,l), body(l)))
    elif re.search(r'\bc_local (wdb_save_l1|wdp_save_l[123])\s+"Caution:', l):
        caut.append((prog[i], cond_for(i,l), body(l)))
    elif re.search(r'\bdi(s|splay)?\s+as\s+err(or)?\b', l):
        m=re.search(r'"([^"]*)"', l); b=m.group(1) if m else l.strip()
        if b.startswith("Note") or b.startswith("Caution"): continue
        errs.append((prog[i], cond_for(i,l), b))
def esc(s): return s.replace("|","\\|").replace("\n"," ")
out=["## 14. Messages: every caution, note and error, with its trigger","",
"This chapter is generated from the program source. Each row gives the subprogram that prints the message, the condition under which it prints (the Stata expression as written in the source, with wooldid's internal macro names), and the message text as written, with macros unexpanded (a name in backquotes is filled in at run time; a message assembled from several quoted pieces shows the pieces). Cautions never stop the run; errors do. Read the condition to learn what triggered a message you saw, and the message for what to do about it.","",
"### Cautions and notes","","| program | prints when | message |","|---|---|---|"]
for p,c,b in caut: out.append(f"| `{p}` | `{esc(c) or '(unconditional at that point)'}` | {esc(b)} |")
out+=["","### Errors and refusals","","Errors stop the run with a return code. A condition mentioning `safety` or `safetyinput` is one that `safety(off)` bypasses.","","| program | refuses when | message |","|---|---|---|"]
for p,c,b in errs: out.append(f"| `{p}` | `{esc(c) or '(see the preceding lines of the program)'}` | {esc(b)} |")
open(sys.argv[2],"w").write("\n".join(out))
print("cautions:", len(caut), "errors:", len(errs), "words:", len(" ".join(out).split()))
