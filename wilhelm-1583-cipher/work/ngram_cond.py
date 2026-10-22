# Conditional 4-gram model P(c4 | c1c2c3) with simple backoff to 3-gram/2-gram, from DTA 1500-1599.
import os, re, sys, collections, unicodedata, math, pickle
root = sys.argv[1]
def norm_text(t):
    t = t.replace('ſ','s').replace('ꝛ','r').replace('aͤ','ä').replace('oͤ','ö').replace('uͤ','ü').replace('ů','u')
    t = unicodedata.normalize('NFC', t).lower()
    return t.replace('ß','ss').replace('ä','e').replace('ö','o').replace('ü','u').replace('j','i').replace('y','i').replace('v','u')
c = [collections.Counter() for _ in range(5)]
for dp, dn, fn in os.walk(root):
    for f in fn:
        if not f.endswith('.txt'): continue
        t = norm_text(open(os.path.join(dp,f), encoding='utf-8', errors='ignore').read())
        for w in re.findall(r"[a-z]+", t):
            s = '__'+w+'_'
            for n in range(1,5):
                for i in range(len(s)-n+1): c[n][s[i:i+n]]+=1
A = 'abcdefghiklmnopqrstuwxz_'
lp = {}
tot1 = sum(c[1].values())
for a in A:
    for b in A:
        for d in A:
            ctx3 = a+b+d
            for e in A:
                g4 = ctx3+e
                # interpolated estimate
                p4 = c[4][g4]/c[3][ctx3] if c[3][ctx3] else 0
                p3 = c[3][b+d+e]/c[2][b+d] if c[2][b+d] else 0
                p2 = c[2][d+e]/c[1][d] if c[1][d] else 0
                p1 = c[1][e]/tot1
                p = 0.6*p4+0.25*p3+0.1*p2+0.05*p1
                lp[g4] = math.log(p) if p>0 else -20
pickle.dump(lp, open('ngram4cond.pkl','wb'))
print(len(lp))
