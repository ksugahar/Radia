"use strict";

const assert = require("node:assert/strict");
const path = require("node:path");

const editor = require(path.join(__dirname, "..", "web", "equation-editor.js"));

assert.deepEqual(
  editor.mathAlphabets.map((item) => item.snippet),
  ["\\mathrm{}", "\\mathit{}", "\\mathbf{}"]
);

let edit;
for (const snippet of ["\\mathsf{}", "\\mathtt{}", "\\mathcal{}",
                       "\\mathbb{}", "\\mathfrak{}", "\\bm{}",
                       "\\mathnormal{}"]) {
  edit = editor.composeInsertion("x", 0, 1, snippet);
  assert.equal(edit.value, snippet.slice(0, -1) + "x}");
}
assert.equal(editor.mathJaxTex("\\bm{\\alpha}+\\bmatrix"),
             "\\boldsymbol{\\alpha}+\\bmatrix");
assert.equal(
  editor.mathJaxTex("\\begin{aligned}short \\\\ muchlonger\\end{aligned}"),
  "\\begin{aligned}& short \\\\ & muchlonger\\end{aligned}"
);
assert.equal(
  editor.mathJaxTex("\\begin{aligned}F&=ma \\\\ E&=mc^2\\end{aligned}"),
  "\\begin{aligned}F&=ma \\\\ E&=mc^2\\end{aligned}"
);

edit = editor.composeInsertion("", 0, 0, "\\mathrm{}");
assert.deepEqual(edit, { value: "\\mathrm{}", caret: 8 });

edit = editor.composeInsertion("abc", 0, 3, "\\mathit{}");
assert.deepEqual(edit, { value: "\\mathit{abc}", caret: 12 });

edit = editor.composeInsertion("a+b", 2, 3, "\\mathbf{}");
assert.deepEqual(edit, { value: "a+\\mathbf{b}", caret: 12 });

/* Existing multi-hole behavior must remain intact: selection fills the first
 * hole and the caret advances into the second. */
edit = editor.composeInsertion("xy", 0, 2, "\\frac{}{}");
assert.deepEqual(edit, { value: "\\frac{xy}{}", caret: 10 });


/* nextHole reaches a hole that ends the source when walking backwards: the
 * caret sits inside the trailing hole at 5, and Shift+Tab must land in the
 * earlier one instead of reporting "no more holes". */
assert.equal(editor.nextHole("a{}b{}", 5, true), 2);
assert.equal(editor.nextHole("{}", 1, true), null);
assert.equal(editor.nextHole("\\frac{}{}", 0, false), 6);
assert.equal(editor.nextHole("\\frac{}{}", 7, false), 8);

/* prettyTex is the native source layout (GUI_SPEC 3.5) applied to inserted
 * templates: break after \begin, before \end, after a row break, and indent
 * by environment depth.  Only TeX-ignored whitespace moves. */
assert.equal(
  editor.prettyTex("\\begin{pmatrix} {} & {} \\\\ {} & {} \\end{pmatrix}"),
  "\\begin{pmatrix}\n  {} & {} \\\\\n  {} & {}\n\\end{pmatrix}"
);
assert.equal(
  editor.prettyTex("\\begin{cases} {} & {} \\\\ {} & {} \\end{cases}"),
  "\\begin{cases}\n  {} & {} \\\\\n  {} & {}\n\\end{cases}"
);
/* An array keeps its column specification on the \begin line, and a nested
 * environment indents one more level. */
assert.equal(
  editor.prettyTex("\\begin{aligned} a &= \\begin{array}{cc} 1 & 2 \\end{array} \\end{aligned}"),
  "\\begin{aligned}\n  a &= \\begin{array}{cc}\n    1 & 2\n  \\end{array}\n\\end{aligned}"
);
/* Escaped braces are control symbols, not group depth. */
assert.equal(editor.prettyTex("\\begin{cases} \\{ x \\} \\end{cases}"),
             "\\begin{cases}\n  \\{ x \\}\n\\end{cases}");

/* Structural Enter (native GUI_SPEC 3.5).  Plain text is wrapped in aligned
 * and split at the caret. */
assert.deepEqual(
  editor.composeRowBreak("E=mc^2", 6, 6),
  { value: "\\begin{aligned}\n  E=mc^2 \\\\\n  {}\n\\end{aligned}", caret: 31 }
);
let broken = editor.composeRowBreak("a+b", 1, 1);
assert.equal(broken.value, "\\begin{aligned}\n  a \\\\\n  +b\n\\end{aligned}");
assert.equal(broken.value.slice(broken.caret, broken.caret + 2), "+b");

/* A selection is replaced by the row break, exactly as typing would. */
broken = editor.composeRowBreak("a+bXY", 3, 5);
assert.equal(broken.value, "\\begin{aligned}\n  a+b \\\\\n  {}\n\\end{aligned}");

/* Inside a row environment the caret's row is split and the column is kept,
 * so the new row lines up under the same alignment point. */
const aligned = "\\begin{aligned}\n  F &= ma \\\\\n  E &= mc^2\n\\end{aligned}";
broken = editor.composeRowBreak(aligned, aligned.indexOf("\n\\end{aligned}"),
                                aligned.indexOf("\n\\end{aligned}"));
assert.equal(
  broken.value,
  "\\begin{aligned}\n  F &= ma \\\\\n  E &= mc^2 \\\\\n  & {}\n\\end{aligned}"
);
assert.equal(broken.value.slice(broken.caret - 1, broken.caret + 1), "{}");

/* Splitting the middle of a row moves the suffix down with its column. */
const matrix = "\\begin{pmatrix}\n  a & b \\\\\n  c & d\n\\end{pmatrix}";
const beforeB = matrix.indexOf("a & b") + 4;
broken = editor.composeRowBreak(matrix, beforeB, beforeB);
assert.equal(
  broken.value,
  "\\begin{pmatrix}\n  a & \\\\\n  & b \\\\\n  c & d\n\\end{pmatrix}"
);

/* A group is never split: inside \frac the row break is added after the whole
 * expression instead of cutting the numerator in half. */
broken = editor.composeRowBreak("\\frac{ab}{c}", 8, 8);
assert.equal(broken.value,
             "\\begin{aligned}\n  \\frac{ab}{c} \\\\\n  {}\n\\end{aligned}");
broken = editor.composeRowBreak(
  "\\begin{aligned}\n  x &= \\frac{ab}{c}\n\\end{aligned}", 30, 30);
assert.equal(
  broken.value,
  "\\begin{aligned}\n  x &= \\frac{ab}{c} \\\\\n  {}\n\\end{aligned}"
);

/* Enter just outside a source that is one row environment adds a row to it
 * rather than nesting a second aligned. */
broken = editor.composeRowBreak(aligned, aligned.length, aligned.length);
assert.equal(
  broken.value,
  "\\begin{aligned}\n  F &= ma \\\\\n  E &= mc^2 \\\\\n  & {}\n\\end{aligned}"
);
broken = editor.composeRowBreak(aligned, 0, 0);
assert.equal(
  broken.value,
  "\\begin{aligned}\n  {} \\\\\n  F &= ma \\\\\n  E &= mc^2\n\\end{aligned}"
);

/* A prime is itself a superscript, so TeX rejects one placed straight after
 * another superscript ("Prime causes double exponent").  The palette inserts
 * at the caret, so the insertion carries MathJax's own remedy when needed and
 * stays untouched when it is not. */
assert.deepEqual(editor.composeInsertion("a", 1, 1, "'"), { value: "a'", caret: 2 });
assert.deepEqual(editor.composeInsertion("f", 1, 1, "''"), { value: "f''", caret: 3 });
assert.deepEqual(editor.composeInsertion("a^{2}", 5, 5, "'"),
                 { value: "a^{2}{}'", caret: 8 });
assert.deepEqual(editor.composeInsertion("a^2", 3, 3, "'"),
                 { value: "a^2{}'", caret: 6 });
assert.deepEqual(editor.composeInsertion("a^{n+1} ", 8, 8, "'"),
                 { value: "a^{n+1} {}'", caret: 11 });
/* Already primed, subscripted, or plain text needs no group. */
assert.deepEqual(editor.composeInsertion("a'", 2, 2, "'"), { value: "a''", caret: 3 });
assert.deepEqual(editor.composeInsertion("x_{1}", 5, 5, "'"), { value: "x_{1}'", caret: 6 });
assert.deepEqual(editor.composeInsertion("a^{2}b", 6, 6, "'"), { value: "a^{2}b'", caret: 7 });
/* An escaped brace is a symbol, not the end of a superscript group. */
assert.deepEqual(editor.composeInsertion("\\{x\\}", 5, 5, "'"),
                 { value: "\\{x\\}'", caret: 6 });

/* An unclosed group is the beginner's most frequent mistake, and MathJax
 * answers it by refusing the whole expression, which left raw TeX on screen.
 * The balance is checked before MathJax sees the source, and the message
 * names where the problem is. */
assert.equal(editor.braceProblem("\\frac{1}{2}"), null);
assert.equal(editor.braceProblem(""), null);
assert.match(editor.braceProblem("\\frac{1}{2"), /閉じていない \{ が 1 個/);
assert.match(editor.braceProblem("\\frac{1}{2"), /^閉じていない/);
assert.match(editor.braceProblem("\\sqrt{\\frac{1}{2}"), /が 1 個/);
/* The first unmatched closing brace is reported by position, 1-based. */
assert.match(editor.braceProblem("x}"), /^2 文字目の \} /);
/* Escaped braces are symbols, and environments balance on their own. */
assert.equal(editor.braceProblem("\\left\\{ x \\right\\}"), null);
assert.equal(editor.braceProblem("\\begin{aligned} a &= b \\end{aligned}"), null);
assert.equal(editor.braceProblem("\\{ \\} \\{"), null);

console.log("PASS: Web style selection wrapping and caret placement");
console.log("PASS: Web hole traversal, native source layout, structural Enter");
