"""Extract text from a PDF with pypdf.  usage: python -I pdf2txt.py in.pdf out.txt"""
import sys
from pypdf import PdfReader

r = PdfReader(sys.argv[1])
with open(sys.argv[2], "w") as fh:
    for i, p in enumerate(r.pages):
        fh.write(f"\n\n===== PAGE {i+1} =====\n")
        try:
            fh.write(p.extract_text() or "")
        except Exception as e:  # noqa
            fh.write(f"[extract error {e}]")
print(len(r.pages), "pages")
