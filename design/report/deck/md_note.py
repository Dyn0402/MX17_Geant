#!/usr/bin/env python3
"""Render a plain markdown record as a self-contained HTML note for dylan-neff.web.cern.ch/notes.

Stdlib only (no `markdown` package on these machines). Handles the subset the records use: #-headings,
paragraphs, - / 1. lists (one level of nesting by two-space indent), | tables |, ``` code blocks,
`code`, **bold**, *italic*, [text](url), and <!--note ... --> metadata passed through.

    python3 md_note.py RECORD.md --out note.html [--title T] [--description D]
then  python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py note.html --slug S --force --deploy
"""
import argparse
import html
import re

CSS = """
:root{--bg:#fbfaf7;--fg:#1d1d1b;--mut:#5f5e58;--line:#dcd9d0;--acc:#1f5f8b;--code:#f0eee8;--hl:#fff4cc}
@media (prefers-color-scheme:dark){:root{--bg:#17181a;--fg:#e6e4de;--mut:#a3a19a;--line:#34363a;--acc:#7fb4dc;--code:#232528;--hl:#3a3420}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
main{max-width:980px;margin:0 auto;padding:24px 16px 80px}
h1{font-size:1.7em;margin:.2em 0 .3em}h2{font-size:1.3em;margin:1.8em 0 .5em;padding-top:.6em;border-top:1px solid var(--line)}
h3{font-size:1.08em;margin:1.4em 0 .4em}
p,li{max-width:76ch}a{color:var(--acc)}
code{background:var(--code);padding:.05em .3em;border-radius:3px;font-size:.88em;word-break:break-word}
pre{background:var(--code);padding:10px 12px;border-radius:5px;overflow-x:auto;font-size:.85em}
pre code{background:none;padding:0}
.tw{overflow-x:auto;margin:.8em 0}
table{border-collapse:collapse;font-size:.88em;font-variant-numeric:tabular-nums}
th,td{border:1px solid var(--line);padding:4px 8px;text-align:left;vertical-align:top}
th{background:var(--code)}
blockquote{margin:1em 0;padding:.4em 1em;border-left:4px solid var(--acc);background:var(--hl)}
.meta{color:var(--mut);font-size:.9em}
"""


def inline(s):
    s = html.escape(s, quote=False)
    s = re.sub(r'`([^`]+)`', lambda m: f'<code>{m.group(1)}</code>', s)
    s = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', s)
    s = re.sub(r'(?<![\w*])\*([^*\s][^*]*)\*(?![\w*])', r'<em>\1</em>', s)
    s = re.sub(r'\[([^\]]+)\]\(([^)\s]+)\)', r'<a href="\2">\1</a>', s)
    s = re.sub(r'(?<!href=")(https?://[^\s<)]+)', r'<a href="\1">\1</a>', s)
    return s


def render(md):
    lines = md.split('\n')
    out, i = [], 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith('```'):
            j = i + 1
            while j < len(lines) and not lines[j].startswith('```'):
                j += 1
            out.append('<pre><code>' + html.escape('\n'.join(lines[i + 1:j])) + '</code></pre>')
            i = j + 1
            continue
        m = re.match(r'(#{1,4}) (.*)', ln)
        if m:
            n = len(m.group(1))
            out.append(f'<h{n}>{inline(m.group(2))}</h{n}>')
            i += 1
            continue
        if ln.startswith('|'):
            rows = []
            while i < len(lines) and lines[i].startswith('|'):
                rows.append([c.strip() for c in lines[i].strip().strip('|').split('|')])
                i += 1
            body = [r for r in rows[1:] if not all(re.fullmatch(r':?-+:?', c) for c in r)]
            t = '<div class="tw"><table><thead><tr>' + ''.join(f'<th>{inline(c)}</th>' for c in rows[0]) + '</tr></thead><tbody>'
            t += ''.join('<tr>' + ''.join(f'<td>{inline(c)}</td>' for c in r) + '</tr>' for r in body)
            out.append(t + '</tbody></table></div>')
            continue
        if re.match(r'\s*([-*]|\d+\.) ', ln):
            # list block with one nesting level
            items = []
            while i < len(lines) and (re.match(r'\s*([-*]|\d+\.) ', lines[i]) or
                                      (lines[i].startswith('   ') and lines[i].strip() and items)):
                l2 = lines[i]
                m2 = re.match(r'(\s*)([-*]|\d+\.) (.*)', l2)
                if m2:
                    items.append((len(m2.group(1)) >= 2, m2.group(2)[0].isdigit(), m2.group(3)))
                else:
                    d, o, t = items[-1]
                    items[-1] = (d, o, t + ' ' + l2.strip())
                i += 1
            ordered = items[0][1]
            h = ['<ol>' if ordered else '<ul>']
            open_sub = False
            for k, (deep, _o, t) in enumerate(items):
                if deep and not open_sub:
                    h[-1] = h[-1].removesuffix('</li>')
                    h.append('<ul>'); open_sub = True
                if not deep and open_sub:
                    h.append('</ul></li>'); open_sub = False
                h.append(f'<li>{inline(t)}</li>')
            if open_sub:
                h.append('</ul></li>')
            h.append('</ol>' if ordered else '</ul>')
            out.append(''.join(h))
            continue
        if ln.startswith('>'):
            q = []
            while i < len(lines) and lines[i].startswith('>'):
                q.append(lines[i][1:].strip()); i += 1
            out.append('<blockquote><p>' + inline(' '.join(q)) + '</p></blockquote>')
            continue
        if not ln.strip():
            i += 1
            continue
        p = []
        while i < len(lines) and lines[i].strip() and not re.match(r'(#{1,4} |\||```|>|\s*([-*]|\d+\.) )', lines[i]):
            p.append(lines[i].strip()); i += 1
        out.append('<p>' + inline(' '.join(p)) + '</p>')
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('md')
    ap.add_argument('--out', required=True)
    ap.add_argument('--title', default=None)
    ap.add_argument('--description', default='')
    a = ap.parse_args()
    md = open(a.md).read()
    title = a.title or re.search(r'^# (.*)', md, re.M).group(1)
    doc = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
           f'<meta name="viewport" content="width=device-width,initial-scale=1">'
           f'<title>{html.escape(title)}</title><meta name="description" content="{html.escape(a.description)}">'
           f'<style>{CSS}</style></head><body><main>{render(md)}</main></body></html>')
    open(a.out, 'w').write(doc)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
