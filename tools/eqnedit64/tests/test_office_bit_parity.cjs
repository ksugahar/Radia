"use strict";
const fs = require("node:fs"), path = require("node:path"),
      http = require("node:http"), os = require("node:os"),
      crypto = require("node:crypto");
const {execFileSync} = require("node:child_process");
const assert = require("node:assert/strict");
const {chromium} = require("playwright");
const web = path.resolve(__dirname, "../web");
const mathjax = path.dirname(require.resolve("mathjax/es5/tex-chtml.js"));
const probe = process.argv[2];
if (!probe)
  throw Error("Usage: node test_office_bit_parity.cjs office_payload_probe");
const wire = process.argv.includes("--wire");
const app = process.argv[3];
if (wire && (process.platform !== "win32" ||
             process.env.EQNEDIT64_ISOLATED_TEST_SESSION !== "1" ||
             [ "INTEL11", "LAB" ].includes(process.env.COMPUTERNAME)))
  throw Error("Wire gate requires disposable Windows isolated host");
if (wire && !app)
  throw Error("Wire gate requires Eqnedit64 app path");
const corpus = JSON.parse(
    fs.readFileSync(path.join(__dirname, "office_bit_corpus.json"), "utf8"));
const fixtureIds = new Set(corpus.map(fixture => fixture.id));
assert.equal(fixtureIds.size, corpus.length, "Duplicate corpus fixture ID");
for (const id of ["H1-H5", "H6-anchored", "H6-unanchored", "cases", "pmatrix",
                  "red", "mathbb", "mathcal", "mathfrak", "mathsf"])
  assert(fixtureIds.has(id), "Missing mandatory fixture: " + id);
const scratch = fs.mkdtempSync(
    path.join(process.env.RUNNER_TEMP ||
                  (process.platform === "win32" ? "C:\\temp" : os.tmpdir()),
              "eqnedit64-office-bits-"));
const server = http.createServer((req, res) => {
  if (req.url.startsWith("/mathjax/")) {
    const file = path.resolve(mathjax, "." + req.url.slice(8));
    if (!file.startsWith(mathjax + path.sep)) {
      res.writeHead(404);
      res.end();
      return;
    }
    res.setHeader("Content-Type", "application/javascript");
    res.end(fs.readFileSync(file));
  } else if (req.url === "/equation-editor.js") {
    res.setHeader("Content-Type", "application/javascript");
    res.end(fs.readFileSync(path.join(web, "equation-editor.js")));
  } else {
    res.setHeader("Content-Type", "text/html; charset=utf-8");
    res.end(
        '<!doctype html><meta charset="utf-8"><script>MathJax={startup:{typeset:false}};SugaharaMath={typeset:()=>Promise.resolve()};</script><script src="/mathjax/tex-chtml.js"></script>' +
        fs.readFileSync(path.join(web, "equation-editor.fragment.html"),
                        "utf8"));
  }
});
// Only transport envelopes/fallback are excluded. Compare raw UTF-8 conditional
// OMML verbatim: namespace, text, scripts, styles, colours, rows all remain.
function officeBytes(html) {
  if (/^Version:/.test(html)) {
    const bytes = Buffer.from(html, "utf8");
    const offset = name => {
      const match = html.match(new RegExp("(?:^|\\r?\\n)" + name + ":(\\d+)"));
      assert(match, "Missing CF_HTML offset: " + name);
      return Number(match[1]);
    };
    const startHtml = offset("StartHTML"), endHtml = offset("EndHTML");
    const start = offset("StartFragment"), end = offset("EndFragment");
    assert(startHtml <= start && start < end && end <= endHtml && endHtml <= bytes.length,
           "Invalid CF_HTML fragment byte range");
    html = bytes.subarray(start, end).toString("utf8");
  }
  const matches = [...html.matchAll(/<!--\[if gte msEquation 12\]>([\s\S]*?)<!\[endif\]-->/g)];
  assert.equal(matches.length, 1, "Expected exactly one primary OMML branch");
  return Buffer.from(matches[0][1], "utf8");
}
(async () => {
  let browser;
  await new Promise(resolve => server.listen(0, "127.0.0.1", resolve));
  try {
    browser = await chromium.launch({
      headless : true,
      ...(process.env.EQNEDIT64_BROWSER_PATH
              ? {executablePath : process.env.EQNEDIT64_BROWSER_PATH}
              : {})
    });
    const page = await browser.newPage();
    if (wire)
      await page.context().grantPermissions(
          [ "clipboard-read", "clipboard-write" ]);
    await page.goto("http://127.0.0.1:" + server.address().port);
    await page.waitForFunction(
        () => !document.querySelector(".eqed-copy-office").disabled);
    if (!wire)
      await page.evaluate(() => {
        window.lastOfficeCopy = null;
        document.execCommand = command => {
          if (command !== "copy")
            return false;
          const data = new DataTransfer();
          const e = new ClipboardEvent(
              "copy", {clipboardData : data, cancelable : true});
          document.dispatchEvent(e);
          window.lastOfficeCopy = {
            html : data.getData("text/html"),
            tex : data.getData("text/plain")
          };
          return e.defaultPrevented;
        };
      });
    let results = [];
    for (const fixture of corpus) {
      const input = path.join(scratch, fixture.id + ".tex");
      fs.writeFileSync(input, fixture.tex);
      const produced = execFileSync(probe, [ input ], {encoding : "utf8"});
      let native = produced;
      if (wire) {
        execFileSync(app, [ input, "office" ], {windowsHide : true});
        const output = path.join(scratch, fixture.id + ".native.cfhtml");
        execFileSync("pwsh", [
          "-NoProfile", "-STA", "-File",
          path.join(__dirname, "read_office_wire.ps1"), "-OutputPath", output
        ]);
        native = fs.readFileSync(output, "utf8");
        assert(officeBytes(native).equals(officeBytes(produced)),
               fixture.id + ": actual EXE clipboard differs from its producer");
      }
      await page.locator(".eqed-source").fill(fixture.tex);
      await page.evaluate(() => window.lastOfficeCopy = null);
      await page.locator(".eqed-copy-office").click();
      const copied = wire ? await page.evaluate(async () => {
        const items = await navigator.clipboard.read();
        for (const item of items)
          if (item.types.includes("text/html"))
            return {html : await (await item.getType("text/html")).text()};
        return null;
      })
                          : await page.evaluate(() => window.lastOfficeCopy);
      assert(copied, fixture.id + ": no Web Office copy");
      if (!wire)
        assert.equal(copied.tex, fixture.tex);
      if (wire) {
        const output = path.join(scratch, fixture.id + ".web.cfhtml");
        execFileSync("pwsh", [
          "-NoProfile", "-STA", "-File",
          path.join(__dirname, "read_office_wire.ps1"), "-OutputPath", output
        ]);
        assert(officeBytes(fs.readFileSync(output, "utf8"))
                   .equals(officeBytes(copied.html)),
               fixture.id +
                   ": browser read differs from raw Windows clipboard");
      }
      const n = officeBytes(native), w = officeBytes(copied.html);
      fs.writeFileSync(path.join(scratch, fixture.id + ".native.html"), native);
      fs.writeFileSync(path.join(scratch, fixture.id + ".web.html"),
                       copied.html);
      const equal = n.equals(w);
      results.push({
        id : fixture.id,
        equal,
        native_bytes : n.length,
        web_bytes : w.length,
        native_sha256 : crypto.createHash("sha256").update(n).digest("hex"),
        web_sha256 : crypto.createHash("sha256").update(w).digest("hex")
      });
    }
    fs.writeFileSync(path.join(scratch, "results.json"), JSON.stringify({
      contract : "eqnedit64.office-omml-bits.v1",
      source_sha : process.env.GITHUB_SHA || null,
      web_sha256 :
          crypto.createHash("sha256")
              .update(fs.readFileSync(path.join(web, "equation-editor.js")))
              .digest("hex"),
      corpus_sha256 : crypto.createHash("sha256")
                          .update(fs.readFileSync(
                              path.join(__dirname, "office_bit_corpus.json")))
                          .digest("hex"),
      mode : wire ? "windows-clipboard-wire" : "producer-bytes",
      results
    },
                                                                        null,
                                                                        2));
    const failures = results.filter(r => !r.equal);
    console.log("Office byte evidence: " + scratch);
    assert.deepEqual(failures, [], "Office payload bit mismatch");
    console.log(
        "PASS: " + results.length +
        " native/Web primary conditional OMML branches are byte-identical (" +
        (wire ? "actual Windows clipboard" : "producer only") +
        "). No Office result inferred.");
  } finally {
    if (browser)
      await browser.close();
    server.close();
  }
})().catch(e => {
  console.error(e);
  process.exitCode = 1;
});
