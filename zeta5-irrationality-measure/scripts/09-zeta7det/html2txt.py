"""Crude HTML->text for WordPress posts (keeps LaTeX alt text).  usage: python -I html2txt.py in.html out.txt"""
import sys, re, html

s = open(sys.argv[1], encoding="utf-8", errors="replace").read()
# keep alt text of latex images
s = re.sub(r'<img[^>]*alt="([^"]*)"[^>]*>', lambda m: " $" + m.group(1) + "$ ", s)
m = re.search(r'<article.*?</article>', s, re.S)
body = m.group(0) if m else s
# comments section too
cm = re.search(r'<ol class="comment-list".*?</ol>', s, re.S)
if cm:
    body += "\n\n===== COMMENTS =====\n" + cm.group(0)
body = re.sub(r'<(script|style)[^>]*>.*?</\1>', '', body, flags=re.S)
body = re.sub(r'<br\s*/?>|</p>|</li>|</h[1-6]>|</div>', '\n', body)
body = re.sub(r'<[^>]+>', '', body)
body = html.unescape(body)
body = re.sub(r'\n\s*\n+', '\n\n', body)
open(sys.argv[2], "w").write(body)
print(len(body))
