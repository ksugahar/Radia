"""
Markdown to HTML converter with MathJax support

Usage:
    python md2html.py <input.md> <output.html> <title>
    python md2html.py <input.md>  # output.html and title derived from input
    md2html <input.md>            # with wrapper command

Features:
    - MathJax support for LaTeX math expressions
    - Table support via markdown extension
    - Code highlighting via codehilite extension
    - Automatic ||...|| to \\Vert notation conversion for markdown table compatibility
    - Fallback encoding support (utf-8 and cp932)
    - Reference link conversion: [1] -> clickable link to #ref1
    - Reference list ID assignment for anchor links
    - Image embedding: local images converted to base64 data URIs (self-contained HTML)
"""
import markdown
import sys
import os
import re
import base64
import mimetypes

def convert_norm_notation(math_content):
    """Convert || notation to \\Vert in math expressions to avoid table conflicts"""
    result = re.sub(r'\\\|\\\|([^|]+)\\\|\\\|', r'\\Vert \1 \\Vert', math_content)
    result = re.sub(r'\|\|([^|]+)\|\|', r'\\Vert \1 \\Vert', result)
    return result

def protect_math(md_content):
    """Protect math expressions from markdown processing"""
    math_blocks = []

    def save_math(match):
        converted = convert_norm_notation(match.group(0))
        math_blocks.append(converted)
        return f"%%MATH_BLOCK_{len(math_blocks)-1}%%"

    def save_math_codeblock(match):
        math_content = match.group(1).strip()
        converted = convert_norm_notation(f"$$\n{math_content}\n$$")
        math_blocks.append(converted)
        return f"%%MATH_BLOCK_{len(math_blocks)-1}%%"

    # Protect ```math code blocks (GitHub style)
    content = re.sub(r'```math\s*\n([\s\S]*?)\n```', save_math_codeblock, md_content)
    # Protect display math ($$...$$)
    content = re.sub(r'\$\$[\s\S]*?\$\$', save_math, content)
    # Protect inline math ($...$)
    content = re.sub(r'\$[^\$\n]+?\$', save_math, content)

    return content, math_blocks

def restore_math(html_content, math_blocks):
    """Restore protected math expressions"""
    for i, math in enumerate(math_blocks):
        html_content = html_content.replace(f"%%MATH_BLOCK_{i}%%", math)
    return html_content

def convert_reference_links(html_content):
    """Convert [N] notation to clickable links pointing to #refN"""
    def replace_ref(match):
        full_match = match.group(0)
        prefix = match.group(1)
        num = match.group(2)
        if 'href=' in prefix or '<a' in prefix:
            return full_match
        return f'{prefix}<a href="#ref{num}" class="ref-link">[{num}]</a>'

    html_content = re.sub(r'((?:<[^>]*>|[^<])*)(\[\d+\])',
                          lambda m: m.group(1) + re.sub(r'\[(\d+)\]', r'<a href="#ref\1" class="ref-link">[\1]</a>', m.group(2))
                          if 'href=' not in m.group(0) and '#ref' not in m.group(0) else m.group(0),
                          html_content)

    html_content = re.sub(r'(?<!href="#ref)(?<!">)\[(\d+)\](?!</a>)',
                          r'<a href="#ref\1" class="ref-link">[\1]</a>',
                          html_content)

    return html_content

def add_reference_ids(html_content):
    """Add id="refN" to ordered list items in References section"""
    ref_section_match = re.search(r'(参考文献|References)', html_content, re.IGNORECASE)
    if not ref_section_match:
        return html_content

    ref_pos = ref_section_match.end()
    before_ref = html_content[:ref_pos]
    after_ref = html_content[ref_pos:]

    ol_match = re.search(r'<ol>([\s\S]*?)</ol>', after_ref)
    if not ol_match:
        return html_content

    ol_start = ol_match.start()
    ol_end = ol_match.end()
    ol_content = ol_match.group(1)

    li_count = [0]
    def add_li_id(_):
        li_count[0] += 1
        return f'<li id="ref{li_count[0]}">'

    ol_with_ids = re.sub(r'<li>', add_li_id, ol_content)

    after_ref = (after_ref[:ol_start] +
                '<ol class="references">\n' + ol_with_ids + '</ol>' +
                after_ref[ol_end:])

    return before_ref + after_ref

def embed_images(html_content, base_dir):
    """Replace local image src with base64 data URIs for self-contained HTML"""
    img_count = 0

    def replace_img_src(match):
        nonlocal img_count
        full_tag = match.group(0)
        src = match.group(1)

        if src.startswith(('http://', 'https://', 'data:')):
            return full_tag

        img_path = os.path.join(base_dir, src)
        if not os.path.isfile(img_path):
            print(f"  Warning: Image not found: {img_path}")
            return full_tag

        mime_type, _ = mimetypes.guess_type(img_path)
        if mime_type is None:
            mime_type = 'application/octet-stream'

        with open(img_path, 'rb') as f:
            img_data = base64.b64encode(f.read()).decode('ascii')

        img_count += 1
        size_kb = os.path.getsize(img_path) / 1024
        print(f"  Embedded: {src} ({size_kb:.0f} KB)")

        return full_tag.replace(f'src="{src}"', f'src="data:{mime_type};base64,{img_data}"')

    result = re.sub(r'<img\s[^>]*src="([^"]+)"[^>]*>', replace_img_src, html_content)
    print(f"  Images embedded: {img_count}")
    return result


def read_file_with_fallback(file_path):
    """Read file with UTF-8, fallback to cp932 for Japanese files"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return f.read()
    except UnicodeDecodeError:
        try:
            with open(file_path, 'r', encoding='cp932') as f:
                return f.read()
        except UnicodeDecodeError:
            print(f"Error: Could not decode {file_path} with utf-8 or cp932.")
            sys.exit(1)

def md_to_html(md_file, output_file=None, title=None):
    """Convert markdown file to HTML with MathJax support"""

    if output_file is None:
        output_file = md_file.replace('.md', '.html')

    md_content = read_file_with_fallback(md_file)

    # Protect math expressions before markdown processing
    protected_content, math_blocks = protect_math(md_content)

    # Convert markdown to HTML using markdown library with extensions
    html_body = markdown.markdown(
        protected_content,
        extensions=['tables', 'fenced_code', 'codehilite']
    )

    # Restore math expressions
    html_body = restore_math(html_body, math_blocks)

    # Add reference IDs to ordered list items
    html_body = add_reference_ids(html_body)

    # Convert [N] to clickable reference links
    html_body = convert_reference_links(html_body)

    # Embed local images as base64 data URIs
    md_dir = os.path.dirname(os.path.abspath(md_file))
    html_body = embed_images(html_body, md_dir)

    # Extract title from first h1 or use provided/filename
    if title is None:
        title = os.path.basename(md_file).replace('.md', '')
        if '<h1>' in html_body:
            start = html_body.find('<h1>') + 4
            end = html_body.find('</h1>')
            title = html_body[start:end]

    # Create full HTML document with MathJax
    html_template = '''<!DOCTYPE html>
<html lang="ja">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script>
    MathJax = {{
      tex: {{
        inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
        displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
        processEscapes: true
      }},
      options: {{
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
      }}
    }};
    </script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.8;
            max-width: 900px;
            margin: 0 auto;
            padding: 20px;
            background-color: #fafafa;
            color: #333;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-bottom: 1px solid #bdc3c7;
            padding-bottom: 5px;
            margin-top: 40px;
        }}
        h3 {{
            color: #555;
            margin-top: 30px;
        }}
        h4 {{
            color: #666;
            margin-top: 20px;
        }}
        code {{
            background-color: #ecf0f1;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Consolas', 'Monaco', monospace;
        }}
        pre {{
            background-color: #2c3e50;
            color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
        }}
        pre code {{
            background-color: transparent;
            padding: 0;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #bdc3c7;
            padding: 10px;
            text-align: left;
        }}
        th {{
            background-color: #3498db;
            color: white;
        }}
        tr:nth-child(even) {{
            background-color: #ecf0f1;
        }}
        blockquote {{
            border-left: 4px solid #3498db;
            margin: 20px 0;
            padding: 10px 20px;
            background-color: #ecf0f1;
        }}
        .highlight {{
            background-color: #fffacd;
            padding: 2px 5px;
            border-radius: 3px;
        }}
        ul, ol {{
            padding-left: 30px;
        }}
        li {{
            margin: 8px 0;
        }}
        a {{
            color: #3498db;
            text-decoration: none;
        }}
        a:hover {{
            text-decoration: underline;
        }}
        hr {{
            border: none;
            border-top: 1px solid #ddd;
            margin: 30px 0;
        }}
        /* Reference link styling */
        .ref-link {{
            color: #e74c3c;
            font-weight: bold;
            text-decoration: none;
        }}
        .ref-link:hover {{
            text-decoration: underline;
        }}
        /* References section styling */
        ol.references {{
            background-color: #f5f5f5;
            padding: 20px 20px 20px 45px;
            border-radius: 8px;
        }}
        ol.references li {{
            padding: 10px;
            background-color: #fff;
            border-radius: 5px;
            border-left: 3px solid #3498db;
            margin: 10px 0;
        }}
        ol.references li:target {{
            background-color: #fffacd;
            border-left-color: #e74c3c;
        }}
        /* Math block styling (for display math) */
        .MathJax_Display {{
            background-color: #f8f9fa;
            padding: 10px;
            margin: 15px 0;
            border-radius: 5px;
            overflow-x: auto;
        }}
    </style>
</head>
<body>
{body}
</body>
</html>
'''

    html_content = html_template.format(title=title, body=html_body)

    with open(output_file, 'w', encoding='utf-8-sig') as f:
        f.write(html_content)

    print(f"Converted: {md_file} -> {output_file}")
    print(f"  Math blocks found: {len(math_blocks)}")
    return output_file

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] in ['--help', '-h']:
        print(__doc__)
        sys.exit(0)

    if len(sys.argv) >= 4:
        md_path = sys.argv[1]
        html_path = sys.argv[2]
        title = sys.argv[3]
        md_to_html(md_path, html_path, title)
    elif len(sys.argv) >= 2:
        md_path = sys.argv[1]
        md_to_html(md_path)
    else:
        print("Usage:")
        print("  python md2html.py <input.md> <output.html> <title>")
        print("  python md2html.py <input.md>  # auto-generate output and title")
        print("  md2html <input.md>            # with wrapper command")
        sys.exit(1)
