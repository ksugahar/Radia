"""Compatibility CLI for the canonical radia_mcp.md2html converter."""

import argparse


def md_to_html(md_file, output_file=None, title=None):
    """Preserve the legacy path return value while using the maintained core."""
    from radia_mcp.md2html import md_to_html as convert

    result = convert(md_file, output_file, title)
    if result.get("log"):
        print(result["log"])
    print(f"Converted: {md_file} -> {result['output_file']}")
    return result["output_file"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output", nargs="?")
    parser.add_argument("title", nargs="?")
    args = parser.parse_args(argv)
    md_to_html(args.input, args.output, args.title)


if __name__ == "__main__":
    main()
