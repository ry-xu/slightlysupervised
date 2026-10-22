import re, sys, collections
exec(open('parse.py').read().split("for t,b in out")[0].replace("src = sys.argv[1]","src = sys.argv[1]"))
words=[]
for t,b in out:
    for tok in re.split(r'[\s|,.()]+', b):
        tok = tok.replace('^','').replace('-','')
        if not tok or tok in '+#': 
            if tok=='#': words.append((t,'#'))
            continue
        words.append((t,tok))
c = collections.Counter(w for t,w in words)
print(len(words), 'words', sum(len(w) for t,w in words), 'signs')
print(c.most_common(40))
# letter freq
lc = collections.Counter(ch for t,w in words for ch in w)
tot=sum(lc.values())
print(sorted(((k,v,round(100*v/tot,1)) for k,v in lc.items()), key=lambda x:-x[1]))
# position stats
first=collections.Counter(w[0] for t,w in words); last=collections.Counter(w[-1] for t,w in words)
print('first', first.most_common()); print('last', last.most_common())
# endings
for n in (2,3):
    print('end',n, collections.Counter(w[-n:] for t,w in words if len(w)>n).most_common(20))
    print('beg',n, collections.Counter(w[:n] for t,w in words if len(w)>n).most_common(20))
bg=collections.Counter()
for t,w in words:
    for i in range(len(w)-1): bg[w[i:i+2]]+=1
print('bigrams', bg.most_common(40))
# what precedes / follows '/'
print('before /', collections.Counter(w[i-1] for t,w in words for i in range(1,len(w)) if w[i]=='/').most_common())
print('after /', collections.Counter((w[i+1] if i+1<len(w) else '$') for t,w in words for i in range(len(w)) if w[i]=='/').most_common())
print('lengths', sorted(collections.Counter(len(w) for t,w in words).items()))
