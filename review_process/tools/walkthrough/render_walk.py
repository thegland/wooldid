#!/usr/bin/env python3
"""Render the changes walk-through page from the workflow's sections JSON.
Usage: render_walk.py <sections.json> <out.html> [--date YYYY-MM-DD] [--range a..b]
sections.json: a list of {area, title, intro, items:[{heading, what_changed, why, decision, what_you_notice, sources, kind?}]}
"""
import html, json, re, sys

def esc(s):
    return html.escape(s or '', quote=True)

def para(s):
    """plain text with blank-line paragraphs -> <p> blocks; single newlines become spaces"""
    parts = [p.strip() for p in re.split(r'\n\s*\n', (s or '').strip()) if p.strip()]
    return ''.join('<p>' + esc(re.sub(r'\s*\n\s*', ' ', p)) + '</p>' for p in parts)

def slug(s):
    return re.sub(r'[^a-z0-9]+', '-', s.lower()).strip('-')[:60]

KIND_LABEL = {'new feature': 'new', 'changed behaviour': 'changed', 'design decision': 'decision', 'removal': 'removed', 'documentation stance': 'guidance'}

def render(sections, date, rng):
    toc = []
    body = []
    n_items = 0
    for si, sec in enumerate(sections):
        sid = 's' + str(si + 1) + '-' + slug(sec['title'])
        toc.append(f'<a href="#{sid}">{esc(sec["title"])}</a>')
        items_html = []
        for it in sec['items']:
            n_items += 1
            kind = (it.get('kind') or '').lower()
            chip = f'<span class="chip k-{slug(kind) or "x"}">{esc(KIND_LABEL.get(kind, kind))}</span>' if kind else ''
            srcs = it.get('sources') or []
            src_html = ('<details class="src"><summary>Sources</summary><ul>' + ''.join('<li>' + esc(s) + '</li>' for s in srcs) + '</ul></details>') if srcs else ''
            items_html.append(f'''<article class="item" id="{sid}-{slug(it["heading"])}">
<h3>{esc(it["heading"])} {chip}</h3>
<div class="parts">
<div class="part"><span class="lab">What changed</span>{para(it["what_changed"])}</div>
<div class="part"><span class="lab">Why</span>{para(it["why"])}</div>
<div class="part"><span class="lab">The decision</span>{para(it["decision"])}</div>
<div class="part"><span class="lab">What you will notice</span>{para(it["what_you_notice"])}</div>
</div>{src_html}
</article>''')
        body.append(f'<section id="{sid}"><h2>{esc(sec["title"])}</h2><div class="intro">{para(sec["intro"])}</div>{"".join(items_html)}</section>')
    head = f'''<title>wooldid Changes Walk-through</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:wght@500;600&family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{{--ground:#F7F5F0;--ink:#1E2126;--muted:#626870;--rule:#DAD6CC;--accent:#2F5D8A;--panel:#FFFFFF;--panel2:#F1EEE7;--focus:#2F5D8A;--chip:#E9E5DA}}
@media (prefers-color-scheme: dark){{:root:not([data-theme="light"]){{color-scheme:dark;--ground:#15181C;--ink:#E7E4DD;--muted:#9CA3AB;--rule:#30363D;--accent:#7FB0DE;--panel:#1C2126;--panel2:#232930;--focus:#7FB0DE;--chip:#2A3038}}}}
:root[data-theme="dark"]{{color-scheme:dark;--ground:#15181C;--ink:#E7E4DD;--muted:#9CA3AB;--rule:#30363D;--accent:#7FB0DE;--panel:#1C2126;--panel2:#232930;--focus:#7FB0DE;--chip:#2A3038}}
*{{box-sizing:border-box}}
body{{background:var(--ground);color:var(--ink);font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:15px;line-height:1.6;padding-inline:16px;padding-block:0 90px}}
.frame{{display:grid;grid-template-columns:16rem minmax(0,1fr);gap:2.5rem;max-width:112ch;margin:0 auto;padding-top:1.4rem}}
nav.toc{{position:sticky;top:calc(env(safe-area-inset-top,0px) + 1rem);align-self:start;max-height:calc(100vh - 3rem);overflow:auto;font-size:.86rem}}
nav.toc .lab{{font-size:.68rem;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);margin:0 0 .5rem .6rem}}
nav.toc a{{display:block;color:var(--muted);text-decoration:none;padding:.22rem 0 .22rem .6rem;border-left:2px solid transparent;line-height:1.35}}
nav.toc a:hover,nav.toc a:focus-visible{{color:var(--accent);border-left-color:var(--accent);outline:none}}
main{{min-width:0;max-width:78ch}}
h1{{font-family:"Source Serif 4",Georgia,serif;font-weight:600;font-size:2rem;line-height:1.15;margin:0 0 .4rem;text-wrap:balance}}
.byline{{color:var(--muted);font-size:.9rem;margin:0 0 1.1rem}}
.preface{{border:1px solid var(--rule);border-radius:4px;padding:.8rem 1rem;margin:0 0 1.6rem;color:var(--ink)}}
.preface p{{margin:0 0 .5rem}} .preface p:last-child{{margin:0}}
h2{{font-family:"Source Serif 4",Georgia,serif;font-weight:600;font-size:1.45rem;line-height:1.2;margin:2.6rem 0 .6rem;padding-top:1.2rem;border-top:1px solid var(--rule);text-wrap:balance}}
section:first-of-type h2{{border-top:0;padding-top:0;margin-top:1rem}}
.intro{{color:var(--muted);max-width:68ch;margin:0 0 1.2rem}} .intro p{{margin:0 0 .5rem}}
article.item{{background:var(--panel);border:1px solid var(--rule);border-radius:6px;padding:1rem 1.15rem;margin:0 0 1rem}}
h3{{font-family:"IBM Plex Sans",system-ui,sans-serif;font-weight:600;font-size:1.05rem;margin:0 0 .7rem;line-height:1.3;display:flex;flex-wrap:wrap;gap:.5rem;align-items:baseline}}
.chip{{font-size:.66rem;letter-spacing:.08em;text-transform:uppercase;font-weight:600;padding:.1rem .45rem;border-radius:3px;background:var(--chip);color:var(--muted)}}
.chip.k-removal{{color:#A63D2F}} .chip.k-new-feature{{color:#2E7D4F}} .chip.k-design-decision{{color:var(--accent)}}
@media (prefers-color-scheme: dark){{:root:not([data-theme="light"]) .chip.k-removal{{color:#E48A7C}} :root:not([data-theme="light"]) .chip.k-new-feature{{color:#7CC496}}}}
:root[data-theme="dark"] .chip.k-removal{{color:#E48A7C}} :root[data-theme="dark"] .chip.k-new-feature{{color:#7CC496}}
.parts{{display:grid;grid-template-columns:1fr;gap:.6rem}}
.part{{border-left:2px solid var(--rule);padding-left:.8rem}}
.part .lab{{display:block;font-size:.68rem;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);margin:0 0 .15rem}}
.part p{{margin:0 0 .4rem;max-width:68ch}} .part p:last-child{{margin:0}}
details.src{{margin-top:.7rem;font-size:.8rem;color:var(--muted)}} details.src summary{{cursor:pointer}} details.src ul{{margin:.3rem 0 0;padding-left:1.2rem}} details.src li{{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.74rem;overflow-wrap:anywhere}}
code{{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.9em}}
a{{color:var(--accent)}}
@media (max-width:860px){{.frame{{grid-template-columns:1fr;gap:1rem}}nav.toc{{position:static;max-height:none;border:1px solid var(--rule);border-radius:4px;padding:.6rem}}}}
@media (prefers-reduced-motion: reduce){{*{{scroll-behavior:auto!important}}}}
</style>
<div class="frame">
<nav class="toc" aria-label="Contents"><p class="lab">Contents</p>{"".join(toc)}</nav>
<main>
<h1>wooldid Changes Walk-through</h1>
<p class="byline">The substantive changes to the program since v1.3, with the reasoning behind each decision. Commits {esc(rng)}, as of {esc(date)}.</p>
<div class="preface"><p>This is a reading guide to what changed in wooldid and why, written for someone who knows applied econometrics but has not followed the work. It covers new features, changed behaviour, design decisions and removals. Plain bug fixes are left out.</p><p>Each entry has four parts: what the program now does, the econometric reason, the decision that was taken and the alternatives it was chosen over, and what a user sees differently. The sources under each entry point at the development records and commits for anyone who wants to check.</p></div>
{"".join(body)}
</main>
</div>
'''
    return head, n_items

if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]
    date = sys.argv[sys.argv.index('--date') + 1] if '--date' in sys.argv else ''
    rng = sys.argv[sys.argv.index('--range') + 1] if '--range' in sys.argv else ''
    sections = json.load(open(src, encoding='utf-8'))
    page, n = render(sections, date, rng)
    open(out, 'w', encoding='utf-8').write(page)
    print('wrote', out, len(page), 'bytes;', len(sections), 'sections,', n, 'items')
