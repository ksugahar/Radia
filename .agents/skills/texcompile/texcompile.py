"""
TeX compile with automatic PDF lock release.
Closes any process locking the output PDF, then runs pdflatex + bibtex.
"""
import os
import sys
import subprocess
import ctypes
import time

def find_texlive():
    """Find TeXLive bin directory."""
    for year in range(2030, 2020, -1):
        path = f"C:\\texlive\\{year}\\bin\\windows"
        if os.path.isdir(path):
            return path
    return None

def unlock_file(filepath):
    """Try to unlock a file by closing processes that hold it open."""
    filepath = os.path.abspath(filepath)
    if not os.path.exists(filepath):
        return True

    # Try to open exclusively to test if locked
    try:
        with open(filepath, 'r+b'):
            pass
        return True  # Not locked
    except (PermissionError, OSError):
        pass

    print(f"  PDF is locked: {os.path.basename(filepath)}")

    # Try closing known PDF viewers by window title
    ps_cmd = (
        "Get-Process | Where-Object {$_.MainWindowTitle -ne ''} | "
        "Select-Object Id,ProcessName,MainWindowTitle | Format-Table -AutoSize"
    )
    result = subprocess.run(
        ["powershell", "-Command", ps_cmd],
        capture_output=True, text=True, encoding='utf-8', errors='replace'
    )

    basename = os.path.splitext(os.path.basename(filepath))[0]
    pdf_viewers = ['Acrobat', 'AcroRd32', 'SumatraPDF', 'FoxitReader', 'FoxitPDFReader']

    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            try:
                pid = int(parts[0])
                proc_name = parts[1]
            except ValueError:
                continue

            # Kill if it's a PDF viewer or if window title matches our file
            if proc_name in pdf_viewers or basename.lower() in line.lower():
                print(f"  Closing process: {proc_name} (PID {pid})")
                try:
                    subprocess.run(
                        ["powershell", "-Command", f"Stop-Process -Id {pid} -Force"],
                        capture_output=True
                    )
                except Exception:
                    pass

    time.sleep(1)

    # Try to delete the locked file
    try:
        os.remove(filepath)
        print(f"  Removed locked PDF.")
        return True
    except (PermissionError, OSError) as e:
        print(f"  WARNING: Could not remove PDF: {e}")
        print(f"  Please close the application holding the file and retry.")
        return False

def run_cmd(cmd, cwd):
    """Run a command and return success status."""
    result = subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True, encoding='utf-8', errors='replace'
    )
    return result.returncode == 0, result.stdout, result.stderr

def compile_tex(tex_file):
    tex_file = os.path.abspath(tex_file)
    tex_dir = os.path.dirname(tex_file)
    tex_base = os.path.splitext(os.path.basename(tex_file))[0]
    pdf_file = os.path.join(tex_dir, tex_base + '.pdf')

    texlive = find_texlive()
    if texlive:
        os.environ['PATH'] = texlive + os.pathsep + os.environ.get('PATH', '')
        print(f"  TeXLive: {texlive}")

    # Unlock PDF if it exists and is locked
    if os.path.exists(pdf_file):
        if not unlock_file(pdf_file):
            return False

    # Clean aux files to avoid corruption
    for ext in ['.aux', '.out']:
        f = os.path.join(tex_dir, tex_base + ext)
        if os.path.exists(f):
            try:
                os.remove(f)
            except OSError:
                pass

    # pdflatex pass 1
    print("  [1/4] pdflatex (1st pass)...")
    ok, out, err = run_cmd(['pdflatex', '-interaction=nonstopmode', tex_file], tex_dir)
    if not ok and 'no output PDF' in (out + err):
        print(f"  ERROR: pdflatex failed. Check {tex_base}.log")
        return False

    # bibtex
    print("  [2/4] bibtex...")
    run_cmd(['bibtex', tex_base], tex_dir)

    # pdflatex pass 2
    print("  [3/4] pdflatex (2nd pass)...")
    run_cmd(['pdflatex', '-interaction=nonstopmode', tex_file], tex_dir)

    # pdflatex pass 3
    print("  [4/4] pdflatex (3rd pass)...")
    ok, out, err = run_cmd(['pdflatex', '-interaction=nonstopmode', tex_file], tex_dir)

    if os.path.exists(pdf_file):
        size_kb = os.path.getsize(pdf_file) / 1024
        # Count pages from log
        log_file = os.path.join(tex_dir, tex_base + '.log')
        pages = '?'
        if os.path.exists(log_file):
            with open(log_file, 'r', encoding='utf-8', errors='replace') as f:
                for line in f:
                    if 'Output written' in line and 'pages' in line:
                        try:
                            pages = line.split('(')[1].split(' ')[0]
                        except (IndexError, ValueError):
                            pass
        print(f"  Done: {os.path.basename(pdf_file)} ({pages} pages, {size_kb:.0f} KB)")
        return True
    else:
        print(f"  ERROR: PDF not produced. Check {tex_base}.log")
        return False

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python texcompile.py <file.tex>")
        sys.exit(1)

    tex_file = sys.argv[1]
    if not os.path.isfile(tex_file):
        print(f"Error: {tex_file} not found")
        sys.exit(1)

    print(f"Compiling: {os.path.basename(tex_file)}")
    success = compile_tex(tex_file)
    sys.exit(0 if success else 1)
