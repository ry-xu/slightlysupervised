# Render the curated reading (reading.tsv) interleaved with the clear words of the transcription.
import csv, re, collections, sys
exec(open('parse.py').read().split("for t,b in out")[0])
rows=list(csv.DictReader(open('reading.tsv'),delimiter='\t'))
byline=collections.defaultdict(list)
for r in rows: byline[r['line']].append(r)
raw = {}
for ln in open(sys.argv[1],encoding='utf-8').read().split('\n'):
    m = re.match(r'^(L\d\d|P3[a-h])\s+(.*)$', ln)
    if m: raw[m.group(1)] = m.group(2)
for t,b in out:
    rs = iter(byline[t])
    body = raw[t]
    body = re.split(r'\s{3,}\(|\s{3,}\{gloss', body)[0]
    body = re.sub(r'\{gloss[^}]*\}', '', body); body = re.sub(r'\{struck[^}]*\}', '', body)
    parts = re.split(r'(\[[^\]]*\])', body)
    outp=[]
    for p in parts:
        if p.startswith('['): outp.append(p); continue
        for tok in p.split():
            core = tok.strip(',.()').replace('^','')
            if not core: outp.append(tok); continue
            r = next(rs)
            rd = r['reading']
            if r['class'] in 'EG': rd = rd+'(?)' if not rd.endswith('?') else rd
            if r['class']=='U': rd = '⟨'+rd+'⟩'
            pre = tok[:len(tok)-len(tok.lstrip('('))]; post = tok[len(tok.rstrip(',.')):]
            outp.append(pre+rd+post)
    print(f"{t:4} {' '.join(outp)}")
