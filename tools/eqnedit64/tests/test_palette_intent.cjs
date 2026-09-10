"use strict";
// Independent expected TeX is shared with the native test, not read from
// production PALETTES. Test the real insertion function with every real key.
const assert = require("node:assert/strict");
const spec = require("../docs/palette_intent.json");
const editor = require("../web/equation-editor.js");
const shared = require("../shared/palettes.json");
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
// This lexical comparison ignores optional single-character grouping only.
// Expected mathematics remains the independently reviewed native oracle.
function spelling(tex) { return tex.replace(/\s/g, "").replace(/\{([a-z0-9])\}/g, "$1"); }
const commands = new Set();
for (let pi = 0; pi < shared.palettes.length; ++pi) {
  const palette = editor.palettes[pi];
  assert.equal(palette.label, shared.palettes[pi].label);
  for (const item of palette.items) {
    const [face, snippet, , preview, command, bodySlot, unavailable] = item;
    assert(!commands.has(command)); commands.add(command);
    assert.equal(typeof preview, "string", command); // explicit empty = literal face
    let expected, values = "", order;
    if (command.startsWith("template.matrix")) {
      const size = command.slice("template.matrix".length);
      assert(spec.native_matrix_sizes.includes(size));
      const [rows, cols] = size.split("x").map(Number);
      values = "abcdefghijklmnopqrstuvwxyz0123456789".slice(0, rows * cols);
      expected = "\\begin{matrix}" + Array.from({length: rows}, (_, r) =>
        values.slice(r * cols, (r + 1) * cols).split("").join("&")).join("\\\\") + "\\end{matrix}";
    } else if (command.startsWith("template.")) {
      const name = command.slice(9), contract = spec.native_templates[name];
      assert(contract, command); assert.equal(face, contract[0]);
      values = contract[1]; expected = contract[2];
      if (["over", "under", "nthroot"].includes(name)) order = [1, 0];
      if (name === "slashfrac") order = [1, 3];
      if (["prime", "dprime", "tprime"].includes(name)) {
        assert.equal("a" + snippet, expected); values = ""; expected = snippet;
      }
    } else if (command.startsWith("style.")) {
      const contract = spec.native_styles[command.slice(6)];
      assert(contract); assert.equal(face, contract[0]); expected = contract[1]; values = "a";
    } else if (command.startsWith("symbol.")) {
      assert(symbols.get(face).some(t => t.trim() === command.slice(7)));
      expected = command.slice(7);
    } else if (command.startsWith("latex.")) {
      assert.equal(spec.native_raw[face], command.slice(6)); expected = command.slice(6);
    } else {
      assert(command.startsWith("matrix."));
      assert(spec.matrix_actions[command.slice(7)]);
      assert(unavailable, "Native-only matrix editing must explain its boundary");
      assert.throws(() => editor.composePaletteInsertion("", 0, 0, item));
      count++; continue;
    }
    let holes = [...snippet.matchAll(/\{\}|\[\]/g)].map(m => m.index + 1);
    let result = snippet;
    const insertions = [...values].map((v, i) => [holes[order ? order[i] : i], v]);
    for (const [position, value] of insertions.sort((a,b) => b[0]-a[0])) {
      assert.notEqual(position, undefined, command);
      result = result.slice(0,position) + value + result.slice(position);
    }
    assert.equal(spelling(result), spelling(expected), command);
    for (const end of [1, 2]) {
      const actual = editor.composePaletteInsertion("LxR", 1, end, item);
      if (holes.length) {
        const pos = holes[bodySlot], selected = end === 2 ? "x" : "";
        assert.equal(actual.value, "L" + snippet.slice(0,pos) + selected + snippet.slice(pos) + (end === 2 ? "R" : "xR"));
      }
    }
    count++;
  }
}
for (let pi = 0; pi < shared.web_extras.length; ++pi) {
  const source = shared.web_extras[pi], actual = editor.palettes[shared.palettes.length + pi];
  source.items.forEach((entry, i) => {
    const [face, snippet] = actual.items[i];
    check(entry.legacy_palette, entry.legacy_face, snippet);
    assert.equal(face, entry.face); count++;
  });
}
assert.equal(commands.size, 245);
const shippingSnippets = new Set(editor.palettes.flatMap(p => p.items.map(i => i[1].replace(/\s/g, ""))));
for (const [category, entries] of Object.entries(spec.web_overrides))
  for (const [face, snippet] of Object.entries(entries))
    assert(shippingSnippets.has(snippet.replace(/\s/g, "")), `Lost legacy Web action: ${category}/${face}`);
assert.equal(editor.palettes.length, shared.palettes.length + shared.web_extras.length);
assert.deepEqual(editor.paletteTabs.slice(0,5).map(g => g.label), ["基本","解析","集合・記号","幾何","ギリシャ"]);
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
