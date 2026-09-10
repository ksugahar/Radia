"use strict";
// Independent expected TeX is shared with the native test, not read from
// production PALETTES. Test the real insertion function with every real key.
const assert = require("node:assert/strict");
const spec = require("../docs/palette_intent.json");
const editor = require("../web/equation-editor.js");
const symbols = new Map(spec.symbols.split(/\s+/).map(entry => {
  const [face, names] = entry.split(":");
  return [face, names.split(",").map(name => "\\" + name + " ")];
}));

function expectedFor(category, face) {
  const special = (spec.web_overrides[category] || {})[face];
  if (special !== undefined) return [special];
  assert(symbols.has(face), `Unreviewed Web key: ${category}/${face}`);
  return symbols.get(face);
}
function check(category, face, snippet) {
  const expected = expectedFor(category, face);
  assert(expected.includes(snippet), `${category}/${face}: ${JSON.stringify(snippet)} != ${JSON.stringify(expected)}`);
  // Both insertion at a caret and replacement/wrapping of a selection.
  // Expected text composition is deliberately independent of composeInsertion.
  for (const end of [1, 2]) {
    const actual = editor.composeInsertion("LxR", 1, end, snippet);
    const selected = end === 2;
    let wanted = snippet;
    if (selected && wanted.includes("{}")) wanted = wanted.replace("{}", "{x}");
    assert.equal(actual.value, "L" + wanted + (selected ? "R" : "xR"), `${category}/${face}: insertion`);
    assert(actual.caret >= 1 && actual.caret <= actual.value.length - 1);
  }
}
let count = 0;
const seen = new Set();
for (const palette of editor.palettes) {
  for (const [face, snippet] of palette.items) {
    const key = palette.label + "/" + face;
    assert(!seen.has(key), `Ambiguous key: ${key}`);
    seen.add(key);
    check(palette.label, face, snippet);
    count++;
  }
}
for (const [category, entries] of Object.entries(spec.web_overrides)) {
  for (const face of Object.keys(entries)) assert(seen.has(category + "/" + face), `Stale intent: ${category}/${face}`);
}
assert.deepEqual(editor.mathAlphabets, [
  {label: "R x", snippet: "\\mathrm{}"},
  {label: "I x", snippet: "\\mathit{}"},
  {label: "B x", snippet: "\\mathbf{}"}
]);
// Stable wrong equations must fail, even though their own round trip succeeds.
assert.throws(() => check("矢印・集合", "←", "\\rightarrow "));
assert.throws(() => check("装飾", "¯□", "\\underline{}"));
assert.throws(() => check("構造", "x²", "_{}"));
assert.throws(() => check("装飾", "sf", "\\mathtt{}"));
assert.throws(() => check("構造", "new unreviewed key", "x"));
console.log(`PASS: ${count} Web keys + 3 persistent styles match independent intent; swaps rejected`);
