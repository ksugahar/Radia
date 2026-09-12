"""Small shared TeX lexical boundary, not macro expansion or interpretation."""
import re


_LITERAL_BEGIN = re.compile(r"\\begin\{(verbatim\*?|Verbatim|lstlisting|minted)\}")


def mask_tex_noncode(text: str) -> str:
    """Mask comments and common literal-code forms, preserving source offsets."""
    output = list(text)

    def mask(start, end):
        output[start:end] = [char if char in "\r\n" else " " for char in text[start:end]]

    i = 0
    while i < len(text):
        if text[i] == "%":
            end = text.find("\n", i)
            end = len(text) if end < 0 else end
            mask(i, end)
            i = end
        elif text[i] == "\\":
            literal = _LITERAL_BEGIN.match(text, i)
            if literal:
                closing = r"\end{" + literal[1] + "}"
                end = text.find(closing, literal.end())
                if end < 0:
                    raise ValueError(f"unterminated literal environment: {literal[1]}")
                end += len(closing)
                mask(i, end)
                i = end
                continue
            verb = re.match(r"\\verb\*?(?![A-Za-z])", text[i:])
            if verb:
                start = i + verb.end()
                if start >= len(text) or text[start].isspace():
                    raise ValueError("missing inline verb delimiter")
                end = text.find(text[start], start + 1)
                if end < 0 or "\n" in text[start:end]:
                    raise ValueError("unterminated inline verb")
                mask(i, end + 1)
                i = end + 1
                continue
            # A control symbol consumes the next character: \%, \\, etc.
            if i + 1 < len(text) and not text[i + 1].isalpha():
                i += 2
            else:
                command = re.match(r"\\[A-Za-z]+", text[i:])
                i += command.end() if command else 1
        else:
            i += 1
    return "".join(output)
