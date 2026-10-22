# Recovered key for R9424 (cipher sign -> plaintext candidates, in order of preference).
# '' = null. Ambiguity reflects either genuine homophony/polyphony or transcription merging of glyphs.
KEY = {
 'a': ['n'],
 'b': [''],            # only at word start before l-words (Lüttich, Leodiensium) and in Philipp...: null
 'c': ['t'],
 'd': ['c','k',''],    # c/k ; dD = ch ; word-final d behaves as a null
 'D': ['h'],
 'e': ['s'],
 'f': ['c','h','b'],   # rare (3x): cunctationes, Hildesheim(H), ?
 'g': [''],            # null
 'i': ['u','v','f'],   # u/v/f ; 'i/' = w (uu) or ff/uv ; 'ii' = fu/w
 'k': ['e'],           # e (incl. Bavarian unrounded ö)
 'm': ['z','m',''],    # z in German words, m in Latin words, null at word edges
 'n': ['m','n','u'],   # mostly m, sometimes n (noch, nit, mangel), u (quam, darvon)
 'o': ['i','y','c'],
 'p': [''],            # null (5 of 6 occurrences)
 'q': ['y','q'],
 'r': ['a','r'],       # mostly a; r in kriegswesen, zuverhindern
 's': ['l'],
 't': ['p'],
 'u': ['a'],
 'w': ['b'],
 'y': [''],            # null (mostly word-initial)
 'z': ['d','r'],       # d and r (probably two glyphs merged in the transcription)
 '3': ['g'],
 '4': ['o'],
 'ö': ['i'],
 '/': ['DOUBLE'],      # doubles preceding letter: s/=ll, c/=tt, e/=ss, a/=nn, n/=mm, z/=rr, i/=w
}
