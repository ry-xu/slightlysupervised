# From-scratch substitution annealer, no cribs: conditional 4-gram score; homophones allowed but each plaintext
# letter limited to <=3 cipher signs; initial key by frequency rank; '/' treated as an ordinary sign.
import re, sys, random, math, pickle, collections
exec(open('parse.py').read().split("for t,b in out")[0])
lp = pickle.load(open('ngram4cond.pkl','rb'))
words=[]
for t,b in out:
    for tok in re.split(r'[\s|,.()]+', b):
        tok = tok.replace('^','').replace('-','').replace('?','')
        if tok and tok not in '+#' and tok!='Sl': words.append(tok)
freq = collections.Counter(c for w in words for c in w)
signs = [s for s,_ in freq.most_common()]
idx = {s:i for i,s in enumerate(signs)}
W = [[idx[c] for c in w] for w in words]
occ = collections.defaultdict(set)
for wi,w in enumerate(W):
    for c in w: occ[c].add(wi)
seed = int(sys.argv[2]); iters = int(sys.argv[3]); random.seed(seed)
PL = list('enirstadhulgocmbfkwzp')  # ENHG letters by frequency (approx.)
PLX = PL + list('qx')
def wscore(key, w):
    p = '__'+''.join(key[i] for i in w)+'_'
    return sum(lp[p[i:i+4]] for i in range(len(p)-3))
key = [PL[i] if i<len(PL) else random.choice(PL[:8]) for i in range(len(signs))]
random.shuffle(key) if seed>100 else None
ws = [wscore(key,w) for w in W]; cur=sum(ws); best=(cur,key[:])
def ok(k):
    c = collections.Counter(k); return max(c.values())<=3
for it in range(iters):
    T = 4.0*(1-it/iters)+0.05
    k2 = key[:]
    if random.random()<0.6:
        a,b = random.sample(range(len(signs)),2); k2[a],k2[b]=k2[b],k2[a]; ch={a,b}
    else:
        a = random.randrange(len(signs)); k2[a]=random.choice(PLX); ch={a}
        if not ok(k2): continue
    aff=set().union(*(occ[x] for x in ch)); delta=0.0; nw={}
    for wi in aff:
        s=wscore(k2,W[wi]); nw[wi]=s; delta+=s-ws[wi]
    if delta>0 or random.random()<math.exp(delta/T):
        key=k2; cur+=delta
        for wi,s in nw.items(): ws[wi]=s
        if cur>best[0]: best=(cur,key[:])
key=best[1]
print('seed',seed,'score',round(best[0],1))
print(' '.join(f"{s}={key[idx[s]]}" for s in signs))
print(' '.join(''.join(key[i] for i in w) for w in W))
