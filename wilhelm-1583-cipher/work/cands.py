import re, sys, json, collections
exec(open('parse.py').read().split("for t,b in out")[0])
lex = json.load(open('lex_1500s.json'))
def norm(w):
    w = w.replace('ß','ss').replace('ä','e').replace('ö','e').replace('ü','u').replace('j','i').replace('y','i').replace('w','vv').replace('f','v').replace('u','v').replace('k','c').replace('z','Z')
    return w
nlex = collections.Counter()
for w,c in lex:
    if len(w)>=1: nlex[norm(w)] += c
nlex = sorted(nlex.items(), key=lambda x:-x[1])
# sign -> regex fragment in normalized space (v=u/v/f, c=c/k, Z=z letter)
S = {
 'a':'n','c':'t','d':'c','D':'h','e':'s','i':'v','k':'e','n':'[nm]','o':'i','r':'[ar]','s':'l','t':'p','w':'b',
 'z':'[dr]','3':'g','4':'o','m':'Z?','g':'','y':'(?:vv)?','p':'.?','q':'.?','f':'.?','b':'.?','d$':'c?','u':'.?','ö':'i','?':'.?','Z':'.?','l':'.?'
}
def rx(cw):
    out=''; prev=''
    for i,ch in enumerate(cw):
        if ch=='/':
            out += '(?:%s)'%prev if prev else ''
            continue
        if ch=='d' and i==len(cw)-1: frag = S['d$']
        else: frag = S.get(ch, '.?')
        out += frag; prev = frag
    return re.compile('^'+out+'$')
words=[]
for t,b in out:
    for tok in re.split(r'[\s|,.()]+', b):
        tok = tok.replace('^','').replace('-','')
        if tok and tok not in '+#': words.append((t,tok))
for t,w in words:
    r = rx(w)
    m = [(x,c) for x,c in nlex if r.match(x)][:6]
    print(t, w, '->', ' '.join(f"{x}:{c}" for x,c in m))
