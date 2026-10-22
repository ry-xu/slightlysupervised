# Build a character 4-gram model (with word-boundary '_') from DTA 1500-1599 plain texts.
import os, re, sys, collections, unicodedata, math, json, pickle
root = sys.argv[1]
def norm_text(t):
    t = t.replace('ſ','s').replace('ꝛ','r').replace('aͤ','ä').replace('oͤ','ö').replace('uͤ','ü').replace('ů','u')
    t = unicodedata.normalize('NFC', t).lower()
    t = t.replace('ß','ss').replace('ä','e').replace('ö','o').replace('ü','u').replace('j','i').replace('y','i').replace('v','u')
    return t
N=4
cnt = collections.Counter()
for dp, dn, fn in os.walk(root):
    for f in fn:
        if not f.endswith('.txt'): continue
        t = norm_text(open(os.path.join(dp,f), encoding='utf-8', errors='ignore').read())
        for w in re.findall(r"[a-z]+", t):
            s = '_'+w+'_'
            for i in range(len(s)-N+1): cnt[s[i:i+N]]+=1
            if len(s)<N: cnt[s.ljust(N,'_')]+=1
tot = sum(cnt.values())
lp = {k: math.log(v/tot) for k,v in cnt.items()}
floor = math.log(0.1/tot)
pickle.dump((N,lp,floor), open('ngram4.pkl','wb'))
print(len(lp), tot)
