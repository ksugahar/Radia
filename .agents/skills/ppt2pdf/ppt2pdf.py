import os
import sys
import glob
import win32com.client
import winerror
from win32com.client.dynamic import ERRORS_BAD_CONTEXT
ERRORS_BAD_CONTEXT.append(winerror.E_NOTIMPL)

def convert_file(ppt, filepath):
    filepath = os.path.abspath(filepath)
    print(f"  Converting: {os.path.basename(filepath)}")
    ppt.Presentations.Open(filepath)
    pdf_path = os.path.splitext(filepath)[0] + '.pdf'
    ppt.ActivePresentation.SaveAs(pdf_path, 32)
    ppt.ActivePresentation.Close()
    print(f"  -> {os.path.basename(pdf_path)}")
    return pdf_path

if len(sys.argv) < 2:
    print("Usage: python ppt2pdf.py <file.pptx|directory> [file2.pptx ...]")
    print("  file.pptx   - Convert a specific pptx file to PDF")
    print("  directory    - Convert all pptx files in the directory")
    sys.exit(1)

PPT = win32com.client.Dispatch("PowerPoint.Application")
converted = 0

for arg in sys.argv[1:]:
    arg = os.path.abspath(arg)
    if os.path.isdir(arg):
        print(f"Converting all pptx in: {arg}")
        for filepath in glob.glob(os.path.join(arg, '*.ppt*'), recursive=False):
            convert_file(PPT, filepath)
            converted += 1
    elif os.path.isfile(arg):
        convert_file(PPT, arg)
        converted += 1
    else:
        print(f"  Warning: not found: {arg}")

PPT.Quit()
print(f"Done. {converted} file(s) converted.")
