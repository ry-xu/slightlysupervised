import re, sys, collections, json
src = sys.argv[1]
lines = open(src, encoding='utf-8').read().split('\n')
# extract P1 cipher lines L08-L30 and P3 lines
out = []
for ln in lines:
    m = re.match(r'^(L\d\d)\s+(.*)$', ln)
    m3 = re.match(r'^(P3[a-h])\s+(.*)$', ln)
    if m:
        tag, body = m.group(1), m.group(2)
    elif m3:
        tag, body = m3.group(1), m3.group(2)
    else:
        continue
    if tag.startswith('L3'):
        if tag != 'L30': continue
    # strip trailing comments: two+ spaces followed by '(' or '{gloss'
    body = re.split(r'\s{3,}\(|\s{3,}\{gloss', body)[0]
    body = re.sub(r'\{gloss[^}]*\}', '', body)
    body = re.sub(r'\{struck[^}]*\}', '', body)
    body = re.sub(r'\[[^\]]*\]', ' | ', body)  # clear text as separator
    out.append((tag, body.strip()))
for t,b in out: print(t, '::', b)
