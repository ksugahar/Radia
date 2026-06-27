#!/usr/bin/env python
"""
Adaptive LaTeX compilation script for URN paper.

Automatically detects when BibTeX needs to be run based on:
1. Missing .bbl file
2. Changed references.bib
3. Undefined citation warnings in log

Usage:
    python compile.py [--clean] [--view]
"""

import subprocess
import os
import sys
import hashlib
from pathlib import Path


def get_file_hash(filepath):
    """Get MD5 hash of file contents."""
    if not filepath.exists():
        return None
    return hashlib.md5(filepath.read_bytes()).hexdigest()


def check_undefined_citations(log_file):
    """Check if log file contains undefined citation warnings."""
    if not log_file.exists():
        return True  # Assume bibtex needed if no log

    content = log_file.read_text(encoding='utf-8', errors='ignore')
    return 'Citation' in content and 'undefined' in content


def check_undefined_references(log_file):
    """Check if log file contains undefined reference warnings."""
    if not log_file.exists():
        return False

    content = log_file.read_text(encoding='utf-8', errors='ignore')
    return 'Reference' in content and 'undefined' in content


def run_command(cmd, cwd):
    """Run command and return success status."""
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True
    )
    return result.returncode == 0


def compile_latex(paper_dir, basename='urn_paper', clean=False, view=False):
    """
    Compile LaTeX document with adaptive BibTeX execution.

    Returns True if compilation succeeded.
    """
    paper_dir = Path(paper_dir)
    tex_file = paper_dir / f'{basename}.tex'
    bib_file = paper_dir / 'references.bib'
    bbl_file = paper_dir / f'{basename}.bbl'
    aux_file = paper_dir / f'{basename}.aux'
    log_file = paper_dir / f'{basename}.log'
    pdf_file = paper_dir / f'{basename}.pdf'
    hash_file = paper_dir / '.bib_hash'

    if not tex_file.exists():
        print(f"Error: {tex_file} not found")
        return False

    if clean:
        print("Cleaning auxiliary files...")
        for ext in ['.aux', '.bbl', '.blg', '.log', '.out', '.toc', '.lot', '.lof']:
            f = paper_dir / f'{basename}{ext}'
            if f.exists():
                f.unlink()
        if hash_file.exists():
            hash_file.unlink()

    # Check if BibTeX is needed
    need_bibtex = False
    reason = ""

    # 1. No .bbl file
    if not bbl_file.exists():
        need_bibtex = True
        reason = "no .bbl file"

    # 2. references.bib changed since last compile
    if not need_bibtex and bib_file.exists():
        current_hash = get_file_hash(bib_file)
        stored_hash = hash_file.read_text().strip() if hash_file.exists() else None
        if current_hash != stored_hash:
            need_bibtex = True
            reason = "references.bib modified"

    # First pdflatex run
    print("\n[1/N] First pdflatex run...")
    if not run_command(['pdflatex', '-interaction=nonstopmode', f'{basename}.tex'], paper_dir):
        print("  Warning: pdflatex returned non-zero (may have warnings)")

    # 3. Check for undefined citations in log
    if not need_bibtex and check_undefined_citations(log_file):
        need_bibtex = True
        reason = "undefined citations in log"

    # Run BibTeX if needed
    if need_bibtex:
        print(f"\n[2/N] Running BibTeX ({reason})...")
        if not run_command(['bibtex', basename], paper_dir):
            print("  Warning: bibtex returned non-zero")

        # Store hash of bib file
        if bib_file.exists():
            hash_file.write_text(get_file_hash(bib_file))

        # Second pdflatex run (to incorporate bibliography)
        print("\n[3/N] Second pdflatex run...")
        run_command(['pdflatex', '-interaction=nonstopmode', f'{basename}.tex'], paper_dir)

    # Check if references need another pass
    if check_undefined_references(log_file):
        print("\n[4/N] Third pdflatex run (resolving references)...")
        run_command(['pdflatex', '-interaction=nonstopmode', f'{basename}.tex'], paper_dir)

    # Final check
    if pdf_file.exists():
        print(f"\nSuccess! PDF generated: {pdf_file}")
        print(f"  Size: {pdf_file.stat().st_size / 1024:.1f} KB")

        # Count pages (rough estimate from log)
        if log_file.exists():
            log_content = log_file.read_text(encoding='utf-8', errors='ignore')
            if 'Output written on' in log_content:
                for line in log_content.split('\n'):
                    if 'Output written on' in line and 'pages' in line:
                        print(f"  {line.strip()}")
                        break

        # Check for remaining warnings
        if log_file.exists():
            log_content = log_file.read_text(encoding='utf-8', errors='ignore')
            warnings = []
            for line in log_content.split('\n'):
                if 'Citation' in line and 'undefined' in line:
                    warnings.append(line.strip())
                elif 'Reference' in line and 'undefined' in line:
                    warnings.append(line.strip())

            if warnings:
                print("\n  Warnings:")
                for w in warnings[:5]:
                    print(f"    {w}")
                if len(warnings) > 5:
                    print(f"    ... and {len(warnings)-5} more")
            else:
                print("  No undefined citations or references.")

        if view:
            print("\nOpening PDF...")
            if sys.platform == 'win32':
                os.startfile(pdf_file)
            elif sys.platform == 'darwin':
                subprocess.run(['open', pdf_file])
            else:
                subprocess.run(['xdg-open', pdf_file])

        return True
    else:
        print(f"\nError: PDF not generated")
        return False


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Compile URN paper with adaptive BibTeX')
    parser.add_argument('--clean', action='store_true', help='Clean auxiliary files before compiling')
    parser.add_argument('--view', action='store_true', help='Open PDF after compilation')
    args = parser.parse_args()

    paper_dir = Path(__file__).parent
    success = compile_latex(paper_dir, clean=args.clean, view=args.view)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
