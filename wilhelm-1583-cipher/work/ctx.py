import re, sys, collections
exec(open('parse.py').read().split("for t,b in out")[0])
words=[]
for t,b in out:
    for tok in re.split(r'[\s|,.()]+', b):
        tok = tok.replace('^','').replace('-','')
        if tok and tok not in '+': words.append((t,tok))
for sign in sys.argv[2:]:
    print(sign, ':', ' '.join(f"{w}" for t,w in words if sign in w))
