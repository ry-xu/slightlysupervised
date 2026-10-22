# Control: lexical coverage (words of >=4 signs decoding to a 16th-c. lexicon word) under the recovered key
# versus keys made by shuffling the same value-lists among the cipher signs (same ambiguity budget).
import random, re, sys
src = sys.argv[1]
exec(open('decode.py').read().split("if __name__=='__main__':")[0])
words=[]
for t,b in out:
    for tok in re.split(r'[\s|,.()]+', b):
        tok = tok.replace('^','').replace('-','')
        if tok and tok not in '+#': words.append(tok)
longw = [w for w in words if len(w.replace('/',''))>=4]
def coverage(K):
    global KEY
    KEY = K
    hit=0
    for w in longw:
        R = readings(w, cap=3000)
        if any(NL.get(nlz(r),0)>0 for r in R if len(r)>=3): hit+=1
    return hit/len(longw)
true = coverage(dict(KEY0))
print('true key coverage of', len(longw), 'words (>=4 signs):', f"{true:.1%}")
random.seed(7)
signs = [s for s in KEY0 if s!='/']
vals = [KEY0[s] for s in signs]
res=[]
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 30):
    v = vals[:]; random.shuffle(v)
    K = dict(zip(signs, v)); K['/']=KEY0['/']
    res.append(coverage(K))
res.sort()
print('shuffled keys: min', f"{res[0]:.1%}", 'median', f"{res[len(res)//2]:.1%}", 'max', f"{res[-1]:.1%}")
