# Reading a 1583 cipher letter to Duke Wilhelm V of Bavaria

In the Bavarian State Archives there is a letter to Duke Wilhelm V, dated 12 June 1583. It starts out in plain German,
and then at line 8 it switches to cipher. It was written during the Cologne War by someone at the court of Wilhelm's
brother Ernst, three weeks after Ernst was elected Archbishop of Cologne.

David Bourdeau maintains a [ledger](https://github.com/dbourdeau/cyphersolver) of historical ciphers, along with
transcriptions of them. This letter (DECODE record R9424) was listed there as "attempted, not solved". The cipher part
had never been read.

We crack the cipher and find that it's a simple letter-substitution cipher with a few quirks. Claude Opus 5.5 recovered
the key overnight, and 86% of the cipher words now read directly.

None of this has been checked against the original images yet. Rather, we depended on the ledger's transcription of the
cipher symbols, made from photos in the DECODE database.

### What does it say?

In short (full English translation and line-by-line reading in [`reading.md`](reading.md)):
- Three men, "Ho., Fe. und Ni.", show *grosse dissimulationes und cunctationes*. None will give way to the others, and
  each wants to command the field himself.
- "Si." rules alone at court, *cum magna indignatione tam Coloniensium quam Leodiensium*, to the great indignation of
  the people of both Cologne and Liège.
- There is a suspicion that a "Susanna" is to be brought down to Cologne or Liège, escorted by someone from Freising.
  The writer urges that this be prevented, since *es wurdt diser orthen ain gross scandalum geben*: it would cause a
  great scandal around here.
- A cathedral provostship at Würzburg is to be resigned to "son Phillipp", presumably Wilhelm's son Philipp Wilhelm.

### What was the cipher?

Each letter is swapped for another, with three quirks:
- a few dummy letters (nulls) at the start or end of words
- a slash meaning "double the previous letter"
- some look-alike symbols merged in the transcription, so a few symbols seem to stand for two letters. This is probably
  why the earlier attempt failed.

### How was it read?

Word endings gave a hint. One symbol pair ends about a fifth of all words, and in German that is almost certainly *-en*.
Short words followed (*den*, *aber*, *die*, *ich*, *von*, *als*, *nit*), and then long words that come out
letter-perfect: *dissimulationes*, *indignationen*, *beschaidenhait*, *scandalum*, *zuverhindern*, *coloniensium*.

Someone at the time also wrote decodings above a few of the cipher words on the second page: "Wirtzburg", "Hildesheim",
"thumbrobst". Where we can check them, they match. The exception is one word glossed "Aichstet", which doesn't decipher
to that.

### Files

- `reading.md`: the English translation, the key, and the line-by-line reading against the cipher
- `notes.md`: the agent's full write-up (method, key, confidence, controls, open points, context)
- `transcription/`: the ledger's transcription of the cipher symbols, by D. Bourdeau (CC BY 4.0)
- `work/`: the scripts (`key.py`, `decode.py`, `control.py`, `anneal3.py`, …), the per-word reading `reading.tsv`, and a
  16th-century German word list (`lex_1500s.json`) built from the Deutsches Textarchiv

To reproduce, from `work/`:

```bash
python3 decode.py ../transcription/transcription.txt
python3 control.py ../transcription/transcription.txt 40
```

The first prints the decoded letter. The second compares the key against 40 shuffled keys.

The from-scratch search (`anneal3.py`) needs a 4-gram model built with `ngram_cond.py` from the DTA 1500–1599 download.
It isn't included because of its size.

Thanks to David Bourdeau for the ledger and the transcriptions, and to the Deutsches Textarchiv for the 16th-century
German.

Text and reading: CC BY 4.0. Code: MIT. The transcription is by D. Bourdeau (CC BY 4.0). The word list is derived from
the Deutsches Textarchiv (CC BY-SA 4.0).
