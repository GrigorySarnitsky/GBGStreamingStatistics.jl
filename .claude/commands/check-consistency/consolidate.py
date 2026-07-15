#!/usr/bin/env python3
# Consolidate a Documenter-built site into reader-text manuscripts for the /check-consistency panel.
#
#   Usage: consolidate.py <build_dir> <out_full.md> <out_siteonly.md>
#     full     — every rendered page, INCLUDING the Design & Internals pages.
#     siteonly — only the user-facing site: the Design & Internals pages (Architecture, Formalization,
#                Testing, Performance, Internals) are withheld, so the structure / cold-read hats judge
#                self-containedness of the user-facing docstrings + API pages without the design docs.
#
# Documenter renders KaTeX client-side, so each formula's raw TeX sits in the HTML text — this emits it
# verbatim (a latex/notation hat reads it as source). Headings, code, tables, admonition titles, and
# lists are preserved so a reviewer reads it like the page. Pages are auto-discovered from the build, so
# the manuscript stays complete when pages are added; each is delimited by a "PAGE:" banner.
#
# Needs beautifulsoup4  (pip install --user beautifulsoup4).
import sys, re, glob, os
from bs4 import BeautifulSoup, NavigableString, Tag

BUILD, OUT_FULL, OUT_SITE = sys.argv[1], sys.argv[2], sys.argv[3]


def title_of(rel):
    d = os.path.dirname(rel)
    return "Home" if d in ("", ".") else d.replace(os.sep, " / ")


def discover(build):
    pages = []
    for path in glob.glob(os.path.join(build, "**", "index.html"), recursive=True):
        rel = os.path.relpath(path, build)
        top = rel.split(os.sep)[0]
        if top in ("search", "assets"):          # skip the search UI and asset dirs
            continue
        pages.append((title_of(rel), path, rel))
    pages.sort(key=lambda p: (p[2].count(os.sep), p[2]))   # Home first, then by depth then name
    return pages


def render(el, out):
    if isinstance(el, NavigableString):
        out.append(str(el)); return
    if not isinstance(el, Tag):
        return
    name, cls = el.name, el.get('class', [])
    if name in ('script', 'style'):
        return
    if name == 'a' and 'docs-heading-anchor' in cls:
        return
    if name in ('h1', 'h2', 'h3', 'h4', 'h5', 'h6'):
        txt = re.sub(r'\s+', ' ', el.get_text()).strip()
        out.append('\n\n' + '#' * int(name[1]) + ' ' + txt + '\n'); return
    if name == 'pre':
        out.append('\n```\n' + el.get_text().rstrip() + '\n```\n'); return
    if name == 'code' and el.parent.name != 'pre':
        out.append('`' + el.get_text() + '`'); return
    if name == 'table':
        out.append('\n' + render_table(el) + '\n'); return
    if name == 'br':
        out.append('\n'); return
    if name == 'li':
        out.append('\n- ')
        for c in el.children:
            render(c, out)
        return
    if name in ('p', 'div', 'section', 'article', 'ul', 'ol', 'blockquote', 'header', 'span',
                'strong', 'em', 'b', 'i', 'sup', 'sub', 'dl', 'dt', 'dd'):
        if 'admonition-title' in cls:
            out.append('\n[ADMONITION] ')
        for c in el.children:
            render(c, out)
        if name in ('p', 'div', 'section', 'blockquote', 'dt', 'dd'):
            out.append('\n')
        return
    for c in el.children:
        render(c, out)


def render_table(tbl):
    rows = []
    for tr in tbl.find_all('tr'):
        cells = tr.find_all(['th', 'td'])
        rows.append(' | '.join(re.sub(r'\s+', ' ', c.get_text()).strip() for c in cells))
    return '\n'.join(rows)


def page_text(path):
    soup = BeautifulSoup(open(path, encoding='utf-8').read(), 'html.parser')
    art = soup.find('article', class_='content')
    out = []
    if art:
        for c in art.children:
            render(c, out)
    return re.sub(r'\n{3,}', '\n\n', ''.join(out))


# The Design & Internals pages are developer-facing rationale, withheld from the SITEONLY manuscript so
# the structure / cold-read hats judge the user-facing surface without the design docs as a crutch.
DESIGN = ('architecture', 'formalization', 'testing', 'performance', 'internals')

pages = discover(BUILD)
full, site = [], []
for title, path, rel in pages:
    banner = f"\n\n===== PAGE: {title} =====\n"
    text = page_text(path)
    full.append(banner); full.append(text)
    if not any(d in rel.lower() for d in DESIGN):
        site.append(banner); site.append(text)

open(OUT_FULL, 'w', encoding='utf-8').write(''.join(full))
open(OUT_SITE, 'w', encoding='utf-8').write(''.join(site))
print(f"pages: {len(pages)}")
print(f"full     : {sum(len(x) for x in full)} chars -> {OUT_FULL}")
print(f"siteonly : {sum(len(x) for x in site)} chars (design & internals withheld) -> {OUT_SITE}")
