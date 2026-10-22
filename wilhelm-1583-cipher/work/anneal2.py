# From-scratch homophonic annealer with incremental scoring (no cribs, no nulls, '/' treated as an ordinary sign).
import re, sys, random, math, pickle, collections
exec(open('parse.py').read().split("for t,b in out")[0])
N, lp, floor = pickle.load(open('ngram4.pkl','rb'))
words=[]
for t,b in out:
    for tok in re.split(r'[\s|,.()]+', b):
        tok = tok.replace('^','').replace('-','').replace('?','')
        if tok and tok not in '+#' and tok!='Sl': words.append(tok)
signs = sorted(set(''.join(words)))
idx = {s:i for i,s in enumerate(signs)}
W = [[idx[c] for c in w] for w in words]
occ = collections.defaultdict(set)
for wi,w in enumerate(W):
    for c in w: occ[c].add(wi)
PL = 'abcdefghiklmnopqrstuwxz'
seed = int(sys.argv[2]); iters = int(sys.argv[3])
random.seed(seed)
def wscore(key, w):
    p = '_'+''.join(key[i] for i in w)+'_'
    return sum(lp.get(p[i:i+N], floor) for i in range(len(p)-N+1))
key = [random.choice(PL) for _ in signs]
ws = [wscore(key,w) for w in W]; cur=sum(ws); best=(cur,key[:])
for it in range(iters):
    T = 15.0*(1-it/iters)+0.3
    j = random.randrange(len(signs)); old = key[j]; new = random.choice(PL)
    if new==old: continue
    key[j]=new
    aff = occ[j]; delta=0.0; newws={}
    for wi in aff:
        s = wscore(key, W[wi]); newws[wi]=s; delta += s-ws[wi]
    if delta>0 or random.random()<math.exp(delta/T):
        for wi,s in newws.items(): ws[wi]=s
        cur+=delta
        if cur>best[0]: best=(cur,key[:])
    else:
        key[j]=old
key=best[1]
freq = collections.Counter(c for w in words for c in w)
print('seed',seed,'score',round(best[0],1))
print(' '.join(f"{s}={key[idx[s]]}" for s,_ in freq.most_common()))
print(' '.join(''.join(key[i] for i in w) for w in W))
