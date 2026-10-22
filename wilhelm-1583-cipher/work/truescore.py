import re, sys, pickle, collections
exec(open('parse.py').read().split("for t,b in out")[0])
lp = pickle.load(open('ngram4cond.pkl','rb'))
words=[]
for t,b in out:
    for tok in re.split(r'[\s|,.()]+', b):
        tok = tok.replace('^','').replace('-','').replace('?','')
        if tok and tok not in '+#' and tok!='Sl': words.append(tok)
H = dict(a='n',b='l',c='t',d='c',D='h',e='s',f='c',g='e',i='u',k='e',m='z',n='m',o='i',p='e',q='i',r='a',s='l',t='p',u='a',w='b',y='e',z='d',Z='e',l='e',**{'3':'g','4':'o','ö':'i'})
def sc(K):
    s=0
    for w in words:
        p='__'+''.join(K[c] for c in w)+'_'
        s+=sum(lp[p[i:i+4]] for i in range(len(p)-3))
    return s
best=None
for zv in 'dr':
  for sl in 'lstenrw':
    K=dict(H); K['z']=zv; K['/']=sl
    v=sc(K)
    if best is None or v>best[0]: best=(v,zv,sl)
print('hand key (single-valued, z=%s, /=%s) score %.1f'%(best[1],best[2],best[0]))
