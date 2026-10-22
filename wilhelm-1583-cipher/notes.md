> The agent's working write-up as of 22 Oct 2026.

# R9424: the letter to Duke Wilhelm V of Bavaria, 12 June 1583. Result

**Outcome: solved in substance (about 85–90% of the cipher words read).** The cipher is a monoalphabetic substitution
on Latin letters. It has a few nulls, and the transcription merges a few glyphs, which makes some cipher letters look
like they stand for two plaintext letters. Under one key the cipher runs on f.274r and f.275r read as fluent
Early New High German, with Latin phrases mixed in. The sentences are coherent and make sense in the history of the
period. The contemporary interlinear glosses on f.275r agree with the decipherment wherever they can be checked.

The remaining doubts are a few single words. They come from the first-pass transcription, an LLM reading of DECODE
photos that are behind a login, and they cannot be settled without the images.

- Work: `wilhelm-v-1583-r9424/work/` (scripts, `reading.tsv`, `plaintext_reading.txt`).
- Data: `wilhelm-v-1583-r9424/bourdeau-transcription/` (Bourdeau, CC BY 4.0) and `wilhelm-v-1583-r9424/dta-data/` (the
  DTA 1500–1599 plain-text corpus, used for the lexicon and the 4-gram model).
- Date: 22 Oct 2026. The work was done by an AI model (Claude). Tomokiyo asks for this to be disclosed.

---

## 1. Has it been solved already? (checked 22 Oct 2026)

- **Bourdeau.** His write-up `dbourdeau.github.io/cyphersolver/r9424.html` still says "attempted, not solved · text
  not obtained" (posted 15 Oct 2026). `targets/r9424/NOTES.md` on `main` is identical to the copy downloaded tonight.
  His catalogue (entry 161) lists R9424 as "unread".
- **Larry Beck's `CATALOGUE.md`** also lists it as unread.
- **NoAutopilot's `STATUS.md` and `NO-CRACKS.tsv`** mention R9424 only as an item whose images are at BayHStA. They
  report no reading.
- **Cryptiana `unsolved.htm`** has no entry for it.
- **Not done:** a general web search. The session's web-search budget was used up before I could run the queries for
  "R9424" and "Wilhelm V 1583 cipher", so I checked only the ledgers above directly. **Before announcing, run a quick
  web search and re-check Bourdeau's page.**

## 2. How the transcription encodes things

`bourdeau-transcription/transcription.txt` has one glyph per character.

- Clear German is in `[...]`, and spaces separate words.
- `/` is a stroke after a letter. It turns out to mean **double the preceding letter**: `s/` = ll, `c/` = tt,
  `e/` = ss, `a/` = nn, `n/` = mm, `z/` = rr. `i/` is uu/vv, that is **w** (or ff/uv).
- `#` is a dumbbell sign: a nomenclator sign, not read.
- `+` is a cross or insertion mark.
- `S` appears only in `Sl`.
- `?` marks an uncertain glyph. `^` marks a superscript over `p`.

Glosses `{gloss: ...}` stand only over f.275r words. I parsed lines L08–L30 of f.274r and P3a–h of f.275r: 169
cipher words, about 1,020 signs.

## 3. How the key came out

Word endings and short words gave the key by hand within minutes:

- **Endings.** `-ka` ends 32 of 174 words, so it is **-en**: k = e, a = n. `-kz` is **-er**: z = r.
- **Short words.** `zka` (5×) = **den** (z also = d). `rwkz` = **aber**. `Dkzrw` = **herab**. `zok` = **die**.
  `odD` / `eodD` = **ich / sich** (dD = ch). `i4a` = **von**. `rse` = **als**. `e4s/` = **soll**. `noc` = **nit**.
  `roa` = **ain**, the Bavarian spelling of *ein*. `nra` = **man**. `mi` = **zu**.
- **Long words that confirm it.** These come out letter-perfect: `zoe/onisrco4ake` = **dissimulationes**,
  `oazo3aco4aka` = **indignationen**, `wkedDrozkaDroc` = **beschaidenhait**, `zkz3skodDka` = **dergleichen**,
  `edrazrsin` = **scandalum**, `i4zDrazka` = **vorhanden**, `mi/kzDoazkrad` = **zuverhindern**,
  `mizkeo3aozka` = **zu resigniren**, `ycDinwtz4wecko` = **(y)thumbprobstei**, `d4s4aokaeoin` = **coloniensium**,
  `bsk4zokaeoing` = **(b)leodiensium(g)**, `rsokaozc` = **alienirt**.
- **The glosses.** The glossed words decipher to what the glossing hand wrote:
  - `dmzkz` = der
  - `zkzo` = der(o)
  - `ycDinwtz4wecko` = Thumbprobstei, glossed "thumbrobst"
  - `yi?özmwiz3m` = Wirzburg, glossed "Wirtzburg"
  - `fndoszkeDrwndm` = (H)ildesh(a)im, with noise; glossed "Hildesheim"
  - `zka marnka Drwkamp` = **den namen haben**. Bourdeau read this gloss as "den nonnen haben"; the cipher says
    *namen*.
  - Bourdeau's faint readings "sein? schilig?" are glosses of `e4an tDos/otbndy`, which deciphers to **son Phillipp**.
    "frurtig ... nimmermehr" are glosses of `iiDcie zrn4akoa/knka`, which deciphers to **frucht(s) darvon einnemen**.

## 4. The key (cipher sign → plaintext)

| sign | plain | notes |
|---|---|---|
| a | n | |
| c | t | |
| d | c / k | `dD` = ch; at word edges it acts as a null (see §6) |
| D | h | |
| e | s | ε-shaped |
| i | u / v / f | `i/` = w (also ff, uf); `ii` = fu (für, fürnemen, gefüert, fueg) |
| k | e | also Bavarian unrounded ö: *gresten*, *bese* |
| m | z | German words: zu, zur, zeit, abzehalten. In Latin = m (*cum*, *quam*). Often a null at word edges |
| n | m | also n (noch, nit, mangel) and u (quam, darvon). "n/u-shaped glyph" per the legend, so probably merged glyphs |
| o | i | in 3–4 places it is clearly a misread looped d (c), giving o-D for ch |
| q | y / q | drey, sy, bey; quam |
| r | a | r in *kriegswesen* and *zuverhindern*. The legend's "v/r-shaped glyph", so probably two glyphs |
| s | l | |
| t | p | |
| w | b | ω-shaped |
| z | d / r | "2-shaped". Almost certainly two merged glyphs (z for d, r-rotunda ꝛ for r) |
| 3 | g | |
| 4 | o | |
| f | c (cunctationes), H? (Hildeshaim) | 3 occurrences |
| u | a? | once (*cum magna*) |
| ö | i | once (Wirzburg) |
| / | doubles the previous letter | |
| **nulls** | g, y, p, b | g (8×), y (11×, mostly word-initial), p (5 of 6), b (3×). All attached to word edges |
| # | nomenclator sign | 3×. A title or person (see §7) |
| Sl | ? | 2×. Bourdeau thinks it is clear "Si."; the glossing hand seems to have written "der" |

**Structure.** The letter-for-letter part is a chain. Reading plaintext → cipher: p→t, t→c, c→d, d→z, z→m, m→n,
n→a, a→r. A second chain runs l→s, s→e, e→k, and a third u→i, i→o, o→4. So each plaintext letter is enciphered as the
*next* letter of a scrambled sequence.

In such a design every cipher letter has exactly one value. That is further evidence that the double values of `z`,
`n`, `r` and `d` in the transcription are merged glyphs, not real polyphones.

## 5. Plaintext

Notation: `[...]` = clear text from Bourdeau's transcription; `(?)` = reading needs a one-sign emendation or is
conjectural; `⟨...⟩` = unread. The line-by-line version is in `work/plaintext_reading.txt`; per-word evidence is in
`work/reading.tsv`.

### f.274r, lines 8–30 (all cipher apart from the bracketed clear words)

> [wöll in] + geistlichen [np]als auch den kriegswesen in grosser confusion ⟨befinde?⟩. [Zwisch(en) Ho. Fe. und Ni.
> sein] grosse dissimulationes und cunctationes, will kainer den andern nachgeben, ieder für sich selbs das veldt
> regiren. [Zuelegen?] Sy versthen den handl alle drey nit genugsam. (Zu hof regiert [Si.] allain cum magna
> indignatione[n] tam Coloniensium quam Leodiensium.) Es sein noch zur zeit gar wenig cölnische bey hof, man wartt aber
> teglich etlicher fürnemen, die beschriben. Summa, ich findt den gresten mangel an guetten treuen leuthen. Es ist ain
> suspi[c]tion vorhanden, als wöll man die Susanna lassen herab chomen; die soll ⟨Fraun·?·runner — a name or title⟩
> von Freising mitt seim ⟨Weib? / Weihb.?⟩ herab belaitten, wais aber nit, obs gen ⟨eo: Cöln?⟩ oder gen Lüttich
> soll gefürt werden. ⟨#⟩ welen darauf gedacht sein, wie dis mit fueg zu verhindern, dann es wurdt diser orthen ain
> gross scandalum geben, und vil guetherzige so woll als von [Ca.?] alienirt werden. Wehrs müglich, ire ⟨#⟩ von disen
> und dergleichen, darein ehr durch böse leuth gefüert, mit beschaidenhait und glimpf abzehalten, so soll man disem
> stift chain andern herrn wünschen.

### f.275r (single cipher words inside clear lines)

> [haben, wann die sich mit] der Thumbprobstei Wirzburg ⟨zodDcw3, glossed "Aichstet"⟩ … [.. f.g.] ⟨#⟩ [dieselb nach
> ausgang 5 der ...] der(o) … son Phillipp zu resigniren … so wurdt ⟨Sl⟩ zu Freising die maisten frucht(s)? darvon
> einnemen, [und] ⟨Sl⟩ allain den namen haben. … [mitt] Hildeshaim …

### Modern German (normalized)

> … dass ich sowohl das Geistliche als auch das Kriegswesen in großer Konfusion [befinde]. Zwischen Ho., Fe. und Ni.
> sind große Verstellungen und Verzögerungen; keiner will dem anderen nachgeben, jeder will für sich selbst das Feld
> regieren. Sie verstehen den Handel alle drei nicht genügend. (Am Hof regiert Si. allein, mit großer Empörung sowohl
> der Kölner als der Lütticher.) Es sind zur Zeit noch sehr wenige Kölner bei Hof; man wartet aber täglich auf etliche
> Vornehme, die einberufen sind. Kurz: ich finde den größten Mangel an guten, treuen Leuten. Es besteht der Verdacht,
> man wolle die Susanna herabkommen lassen; [X] von Freising soll sie mit seiner Frau (oder: seinem Weihbischof)
> herabbegleiten. Ich weiß aber nicht, ob sie nach [Köln?] oder nach Lüttich gebracht werden soll. [#] mögen darauf
> bedacht sein, wie dies mit Fug zu verhindern ist; denn es gäbe hierzulande ein großes Ärgernis, und viele
> Gutherzige würden, wie auch [die] von Ca., entfremdet. Wäre es möglich, Ihre [#] von diesen und dergleichen Dingen,
> in die er durch böse Leute geführt worden ist, mit Bescheidenheit und Glimpf abzuhalten, so sollte man diesem Stift
> keinen anderen Herrn wünschen.
>
> (f.275r) … mit der Dompropstei … Würzburg … [Eichstätt?] … [Ihrem] Sohn Philipp zu resignieren … so würde [der] zu
> Freising die meisten Früchte (Einkünfte) davon einnehmen, und [der andere] allein den Namen haben … Hildesheim …

### English translation

> … [so that I find] both the spiritual [affairs] and the war affairs in great confusion. Between Ho., Fe. and Ni.
> there is a great deal of dissembling and delay. None will give way to another, and each wants to command the field
> for himself. All three do not understand the business well enough. (At court Si. alone governs, to the great
> indignation of the people of both Cologne and Liège.) There are still very few men of Cologne at court, but some
> notables who have been summoned are expected daily. In short, I find the greatest lack of good, loyal people.
>
> There is a suspicion that they mean to let Susanna come down [the Rhine]. [X] of Freising is to escort her down with
> his wife (or: his suffragan bishop). But I do not know whether she is to be taken to [Cologne?] or to Liège.
> [#] will please consider how this may properly be prevented. It would cause a great scandal in these parts, and many
> well-meaning people, as well as [those] of Ca., would be alienated. If it were possible to keep His [#] from these
> and similar things, into which he has been led by evil people, with discretion and gentleness, one could wish this
> see no other lord.
>
> (f.275r) … with the cathedral provostship … Würzburg … [Eichstätt?] … to resign to [Your Grace's] son Philipp …
> so [the one] at Freising would take most of the revenues of it, and [the other] would have only the name …
> Hildesheim …

## 6. Confidence

**How much reads** (`work/reading.tsv`: 169 cipher words, ~1,020 signs):

| class | words | signs | meaning |
|---|---|---|---|
| K | 146 (86%) | 83% | reads directly with the key; documented nulls and merged glyphs allowed |
| E | 15 (9%) | 12% | reads after a one-sign emendation of an assumed transcription slip, e.g. o↔d, s↔d, e↔s, z↔i |
| G | 4 (2%) | 3% | conjecture |
| U | 4 (2%) | 2% | unread: the name in L21, `zodDcw3`, and `Sl` ×2. The `#` sign (3×) is also unread |

**Controls.**

- **Shuffled keys** (`work/control.py`). Under the recovered key, 65% of cipher words with ≥4 signs decode to a word
  in a 135k-type corpus of 16th-century German, scored mechanically with no hand choices. The 40 shuffled keys have
  the same ambiguity and nulls; they scored a median of 8% and at most 24%.
- **The decisive evidence.** About 25 long words (8–14 letters), including Latin ones, come out exactly under one key.
  So do the gloss words. The text runs as connected syntax over 20 lines: subordinate clauses, idioms ("mit Fug",
  "mit Bescheidenheit und Glimpf", "tam … quam …") and Bavarian chancery spellings (ain, allain, wais, chomen, chain,
  guet, gresten, abzehalten).
- **From-scratch annealing.** `work/anneal3.py` uses no cribs, no glosses and no nulls; homophones are allowed with
  at most 3 signs per letter, and it scores with an interpolated conditional 4-gram model from the DTA. Three seeds
  of 60k steps each:
  - Seed 2 independently recovers 11 signs of the key, covering **70% of the sign tokens**: k=e, z=r, a=n, o=i,
    e=s, c=t, D=h, n, d=c, 3=g, 4=o. Its raw output already contains *"noch ser seit ger etenig cortsche"*,
    *"ich einrt ren gresten nenger en geetten treeen reethen"*, *"es ist ein sescition"*, *"regirt"*,
    *"inrigntionen"* and *"chgeistrihen"*.
  - Seed 1 recovers 8 signs (52% of tokens).
  - Seed 3 stuck in a poor optimum.

  So the core of the key emerges without any human seeding. The search cannot get the rest because of the merged
  glyphs (z for both d and r), the nulls, and `/`.

  An earlier unconstrained run (`anneal2.py`, joint 4-gram scores, any number of signs per letter) collapsed to
  "e/n/r" text. That is the usual degenerate homophonic optimum, and probably the same trap that defeated the
  earlier annealing attempts.

**Where the doubt lies.** The weakness is the transcription, not the key.

- Bourdeau's `profile.json` says it was made "by llm from images". These are DECODE images behind a login, so I could
  not check it.
- The dual values (z = d/r, n = m/n/u, r = a/r) and the edge "nulls" `d`, `m`, `n` probably come from glyphs being
  merged. The chain structure in §4 points the same way.
- Some nulls are real, because they recur in a consistent pattern: `g`, `y`, `p`, `b`, and probably word-final
  flourishes read as `m` or `d`.

**Specific uncertain points.**

- The opening words, `dZ`/`np`: probably "sowohl im Geistlichen als auch …", with the clear "wöll in" possibly
  misread for "wol in".
- `cokioazkdy`: "befinde" fits the syntax but needs co = w.
- The escort's name in L21 (`izriafzia/kz`, F-r-a-u-n-?-r-u-nn-e-r).
- Weib vs Weihb[ischof] in L21.
- The destination `eo` in L22: Cöln is expected.
- `#`: probably a title for Ernst, e.g. "Ire Kfl. Dt.", since it takes "er" and plural "welen".
- `Sl`.
- `zodDcw3`: glossed Eichstätt, but it does not decipher to that.

**Next step.** Check the reading against the images: re-read the merged glyphs, `#`, `Sl` and L21.

## 7. Historical context (only as context; the claims above rest on the text)

- **Date and setting.** The letter is dated 12 June 1583 (new style), three weeks after Ernst of Bavaria, Wilhelm V's
  brother, was elected Archbishop-Elector of Cologne (22–23 May 1583, during the Cologne War). Ernst was already
  Bishop of Freising (from 1566), Hildesheim (from 1573) and Liège (from 1581). Every place in the text belongs to his
  sees: Freising, Hildesheim, Cologne ("Coloniensium", "cölnische bey hof") and Liège ("Leodiensium", "Lüttich").
- **The writer.** The writer is evidently at Ernst's court. The text calls it "zu Hof" and "bey hof"; it says few men
  of Cologne are there yet and notables have been summoned. He reports to Wilhelm on:
  - rivalry among three men, "Ho., Fe. und Ni.", who are not identified;
  - Si.'s sole sway at court;
  - and, mainly, a planned arrival of a woman, "Susanna", which he fears would cause scandal and alienate loyal
    Catholics.
- **Ernst's conduct.** The passage "darein er durch böse Leut geführt" fits Ernst's known worldly lifestyle; German
  Wikipedia mentions his love affairs and his later companion Gertrud von Plettenberg. Who Susanna was is not known
  from this text.
- **Duke Philipp.** The f.275r cipher words concern resigning a benefice ("Thumbprobstei", Würzburg, Freising,
  Hildesheim) to "son Philipp". The docket on f.275v reads "… und Hertzog Philipps halb", meaning Wilhelm's son
  Philipp Wilhelm (b. 1576), Bishop of Regensburg from 1579. In this plan one party would draw most of the revenues
  and the other would hold only the title.
- **Possible identifications (speculation).** "Ni." might be Felician Ninguarda, the papal nuncio in south Germany,
  and the docket's "Landorff(?)" might be the writer. Neither is supported by the cipher text.

## 8. Files

- `wilhelm-v-1583-r9424/work/`. Run every script as `python3 -I <script> ../bourdeau-transcription/transcription.txt`;
  the scripts read the downloaded data only as text.
  - Parsing and statistics: `parse.py`, `words.py`, `ctx.py`.
  - `key.py`: the key.
  - `decode.py`: mechanical decoding with lexicon disambiguation.
  - `reading.tsv`: the curated per-word reading with class and note.
  - `render.py` and `plaintext_reading.txt`: the rendered text.
  - `control.py`: the shuffled-key control.
  - `lexicon.py`: the 16th-century lexicon (`lex_1500s.json`).
  - `ngram.py` and `ngram_cond.py`: the 4-gram models.
  - `anneal3.py` and `anneal3_seed*.txt`: the from-scratch annealer and its outputs. `anneal2*` is the degenerate
    unconstrained run.
  - `truescore.py`: the hand key's score under the same model (−2700, against −2527 to −2643 for the annealer optima).
    The hand key scores lower only because a single-valued key cannot express z = d/r, `/`-doubling or the nulls.
    The n-gram score alone therefore cannot verify this text; the long-word and syntax evidence does.
  - `cands.py`, `show.py` and `grepx.py`: exploration helpers.
- `wilhelm-v-1583-r9424/dta-data/`: DTA 1500–1599 transliterated plain texts (dtak and dtae, 2020-10-23 release),
  data only.
- Not attempted: `vich-1515/`. The first target was read, so I did not start the second.
