# build word frequency lexicon from DTA 1500-1599 plain texts (data only; read as text)
import os, re, sys, collections, unicodedata, json
root = sys.argv[1]
cnt = collections.Counter()
for dp, dn, fn in os.walk(root):
    for f in fn:
        if not f.endswith('.txt'): continue
        t = open(os.path.join(dp,f), encoding='utf-8', errors='ignore').read()
        t = t.replace('ſ','s').replace('ꝛ','r').replace('ů','u').replace('ÿ','y').replace('aͤ','ä').replace('oͤ','ö').replace('uͤ','ü')
        t = unicodedata.normalize('NFC', t).lower()
        for w in re.findall(r"[a-zäöüß]+", t):
            cnt[w]+=1
json.dump(cnt.most_common(), open('lex_1500s.json','w'))
print(len(cnt), sum(cnt.values()))
print(cnt.most_common(80))
