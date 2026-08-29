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

console.log("PASS: Web style selection wrapping and caret placement");
