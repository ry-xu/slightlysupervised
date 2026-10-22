import re, sys
exec(open('parse.py').read().split("for t,b in out")[0])
key = dict(k='e',a='n',z='R',o='i',e='s',**{'4':'o'},i='v',r='a',w='b',D='h',n='m',c='t',s='l',**{'3':'g'},d='c')
def dec(w):
    s=''
    for ch in w:
        if ch=='/': s+= '/'
        else: s+= key.get(ch, ch.upper() if ch.islower() else '_'+ch)
    return s
for t,b in out:
    toks = re.split(r'(\s+)', b)
    print(t, b)
    print('   ', ''.join(dec(x) if x.strip() and x not in '|,.+#()' else x for x in toks))
