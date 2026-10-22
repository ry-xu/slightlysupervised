import re, sys, json
lex = json.load(open('lex_1500s.json'))
pat = re.compile(sys.argv[1])
n=0
for w,c in lex:
    if pat.search(w):
        print(w,c); n+=1
        if n>40: break
