#!/usr/bin/env python3
"""Generate the guide's messages chapter from the ado. Usage: gen_messages_chapter.py <ado> <out.md>

A forward pass keeps a stack of the open blocks in each program, so every message is listed with
the conjunction of the if/else conditions that enclose it (loops and quietly/capture blocks add
nothing). An else branch shows the negation of the if (or if chain) it closes."""
import re, sys
ado=open(sys.argv[1],encoding="utf-8").read().split("\n")

def code_of(line):
    """the line with string literals and comments removed, for brace counting"""
    s=line
    if re.match(r'\s*\*', s): return ""
    s=re.sub(r'"[^"]*"', '""', s)
    s=re.sub(r'\s///?.*$', '', s)
    s=re.sub(r'/\*.*?\*/', '', s)
    return s

HEAD_IF=re.compile(r'\s*(?:\}\s*)?(else\s+)?if\s+(.*?)\s*\{\s*$')
HEAD_ELSE=re.compile(r'\s*(?:\}\s*)?else\s*\{\s*$')
INLINE_IF=re.compile(r'\s*(else\s+)?if\s+(.*?)\s+(?:(?:noi(?:sily)?|qui(?:etly)?|cap(?:ture)?)\s+)*(?:di|dis|display|local|c_local|exit|error)\b')
INLINE_ELSE=re.compile(r'\s*else\s+(?:(?:noi(?:sily)?|qui(?:etly)?)\s+)*(?:di|dis|display|local|c_local)\b')

def neg(conds):
    return " & ".join("not (" + c + ")" for c in conds)

prog=[None]*len(ado); ctx=[None]*len(ado); inl=[None]*len(ado)
cur=None; stack=[]; closed={}      # closed[depth] = conditions of the if-chain that just closed there
last_inline=None                   # conditions of the inline if-chain on the lines just above
for i,l in enumerate(ado):
    m=re.match(r'\s*prog(?:ram)?\s+def(?:ine)?\s+([A-Za-z_0-9]+)', l)
    if m: cur=m.group(1); stack=[]; closed={}; last_inline=None
    prog[i]=cur
    c=code_of(l)
    # the context a message on this line inherits: the open blocks BEFORE the line's own braces,
    # unless the line itself opens a branch (then the branch applies, handled below)
    mi=HEAD_IF.match(c) and HEAD_IF.match(re.sub(r'\s///?.*$', '', l).rstrip()); me=HEAD_ELSE.match(c)
    ml=re.match(r'\s*(\}*)', c); lead=len(ml.group(1))
    for _ in range(lead):
        if stack:
            e=stack.pop(); d=len(stack)
            closed[d]=e['chain'] if e['kind']=='if' else None
    d=len(stack)
    if mi or me:
        if (mi and mi.group(1)) or me:            # else-if / else: continues the chain closed at d
            prev=closed.get(d) or []
            if mi:
                cond=mi.group(2).strip(); chain=prev+[cond]
                shown=(neg(prev)+" & " if prev else "")+cond
            else:
                chain=None; shown=neg(prev) if prev else "(else branch of the preceding condition)"
            stack.append({'kind':'if' if mi else 'else','cond':shown,'chain':chain})
        else:
            cond=mi.group(2).strip()
            stack.append({'kind':'if','cond':cond,'chain':[cond]})
        closed.pop(d,None)
        ctx[i]=[e['cond'] for e in stack if e['cond']]
        continue
    ctx[i]=[e['cond'] for e in stack if e['cond']]
    # inline if / else on a message line
    mii=INLINE_IF.match(l); mie=INLINE_ELSE.match(l)
    if mii:
        cond=mii.group(2).strip()
        if mii.group(1) and last_inline:
            inl[i]=neg(last_inline)+" & "+cond; last_inline=last_inline+[cond]
        else:
            inl[i]=cond; last_inline=[cond]
    elif mie:
        prev=last_inline or closed.get(d)
        inl[i]=neg(prev) if prev else "(else branch of the preceding condition)"; last_inline=None
    elif c.strip():
        last_inline=None
    # any other braces on the line (loops, quietly {, capture noisily {, inline { })
    rest=c[ml.end():]
    net=rest.count('{')-rest.count('}')
    for _ in range(max(net,0)): stack.append({'kind':'other','cond':None,'chain':None})
    for _ in range(max(-net,0)):
        if stack: stack.pop()
    if net!=0 or lead: closed={k:v for k,v in closed.items() if k<=len(stack)}
    if re.match(r'\s*end\s*$', l): stack=[]; closed={}

def cond_for(i):
    parts=list(ctx[i] or [])
    if inl[i]: parts.append(inl[i])
    return " & ".join("(" + p + ")" if (" & " in p or " | " in p) and len(parts)>1 else p for p in parts)

def body(l):
    s=l.strip()
    m=re.search(r'"((?:Caution|Note)[:\s].*)"\s*$', s)   # to the LAST quote: keeps embedded expression pieces
    if m: return m.group(1)
    m=re.search(r"local \S+ (.*)$", s); return m.group(1) if m else s
caut=[]; errs=[]
for i,l in enumerate(ado):
    if re.match(r'\s*\*', l): continue
    if re.search(r'\bdi(s|splay)?\b', l) and re.search(r'"(Caution|Note)[:\s]', l):
        caut.append((prog[i], cond_for(i), body(l)))
    elif re.match(r'\s*local (wdbc`wdbncaut\'|wdpc`ncaut\'|imdegmsg)\s', l):
        caut.append((prog[i], cond_for(i), body(l)))
    elif re.search(r'\bc_local (wdb_save_l1|wdp_save_l[123])\s+"Caution:', l):
        caut.append((prog[i], cond_for(i), body(l)))
    elif re.search(r'\bdi(s|splay)?\s+as\s+err(or)?\b', l):
        m=re.search(r'"([^"]*)"', l); b=m.group(1) if m else l.strip()
        if b.startswith("Note") or b.startswith("Caution"): continue
        errs.append((prog[i], cond_for(i), b))
def esc(s): return s.replace("|","\\|").replace("\n"," ")
out=["## 14. Messages: every caution, note and error, with its trigger","",
"This chapter is generated from the program source. Each row gives the subprogram that prints the message, the condition under which it prints (the Stata expressions of every if and else block that encloses it, joined with &, as written in the source with wooldid's internal macro names; an early exit above the message is not shown), and the message text as written, with macros unexpanded (a name in backquotes is filled in at run time; a message assembled from several quoted pieces shows the pieces). Cautions never stop the run; errors do. Read the condition to learn what triggered a message you saw, and the message for what to do about it.","",
"### Cautions and notes","","| program | prints when | message |","|---|---|---|"]
for p,c,b in caut: out.append(f"| `{p}` | `{esc(c) or '(unconditional at that point)'}` | {esc(b)} |")
out+=["","### Errors and refusals","","Errors stop the run with a return code. A condition mentioning `safety` or `safetyinput` is one that `safety(off)` bypasses.","","| program | refuses when | message |","|---|---|---|"]
for p,c,b in errs: out.append(f"| `{p}` | `{esc(c) or '(see the preceding lines of the program)'}` | {esc(b)} |")
open(sys.argv[2],"w").write("\n".join(out))
print("cautions:", len(caut), "errors:", len(errs), "words:", len(" ".join(out).split()))
