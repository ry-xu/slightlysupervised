# Mechanical decode: per cipher word, enumerate all readings under KEY and pick the best by lexicon frequency
# (fallback: first-choice letters). Prints line-by-line.
import re, sys, json, itertools, math, collections
exec(open('parse.py').read().split("for t,b in out")[0])
exec(open("key.py").read()); KEY0=dict(KEY)
lex = dict(json.load(open('lex_1500s.json')))
extra = {}  # Bavarian / Latin forms not in the printed corpus
for w in ("ain allain wais chomen chain kainer guetten guetherzige beschaidenhait abzehalten gresten "
          "belaitten belaiten alienirt dissimulationes cunctationes indignationen indignatione coloniensium leodiensium "
          "quam tam cum magna scandalum suspicion suspition confusion resigniren thumbprobstei wirzburg hildeshaim "
          "lüttich luttich cölnische colnische wartt fueg darvon einnemen fruchts phillipp regirt orthen namen herrn").split():
    extra[w]=50
def nlz(w):
    return w.replace('ü','u').replace('ö','e').replace('ä','e').replace('j','i').replace('y','i')
NL = collections.Counter()
for w,c in list(lex.items())+list(extra.items()):
    NL[nlz(w)] += c
def readings(cw, cap=20000):
    opts=[]
    for i,ch in enumerate(cw):
        if ch=='/':
            opts.append(['DOUBLE','w_']) ; continue
        opts.append(KEY.get(ch, ['?']))
    res=set()
    for n,combo in enumerate(itertools.product(*opts)):
        if n>cap: break
        s=''; prev=''
        for c in combo:
            if c=='DOUBLE': s+=prev
            elif c=='w_':
                if prev in ('u','v','f'): s=s[:-1]+'w'
                else: s+=prev
            else:
                s+=c; prev=c if c else prev
        res.add(s)
    return res
def best(cw):
    R = readings(cw)
    scored = sorted(((NL.get(nlz(r),0), r) for r in R), reverse=True)
    if scored and scored[0][0]>0: return scored[0][1], True
    # fallback first choices
    s=''; prev=''
    for ch in cw:
        if ch=='/': s+=prev; continue
        c = KEY.get(ch,['?'])[0]; s+=c; prev=c or prev
    return s.upper(), False
if __name__=='__main__':
    tot=hit=0; sig_tot=sig_hit=0
    for t,b in out:
        toks = re.split(r'(\s+|\|)', b)
        line=[]
        for x in toks:
            y = x.strip()
            if not y or y in '|,.+#()' :
                line.append(x); continue
            core = y.strip(',.()').replace('^','').replace('-','')
            if not core: line.append(x); continue
            r, ok = best(core)
            tot+=1; hit+=ok; sig_tot+=len(core); sig_hit+=len(core)*ok
            line.append(r if ok else r)
        print(f"{t:4} {''.join(line)}")
    print(f"\nwords decoded to a lexicon word: {hit}/{tot} = {hit/tot:.1%}; by signs {sig_hit}/{sig_tot} = {sig_hit/sig_tot:.1%}")
