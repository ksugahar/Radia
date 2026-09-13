"use strict";
// Real browser gate: successful conversion alone does not prove that CHTML
// installed its font CSS, or that noundefined did not print red TeX text.
const fs = require("node:fs");
const path = require("node:path");
const http = require("node:http");
const assert = require("node:assert/strict");
const {chromium} = require("playwright");
const web = path.resolve(__dirname, "../web");
const mathjax = path.dirname(require.resolve("mathjax/es5/tex-chtml.js"));
const server = http.createServer((req, res) => {
  res.setHeader("Content-Type", "application/javascript; charset=utf-8");
  if (req.url.startsWith("/mathjax/")) {
    const file = path.resolve(mathjax, "." + req.url.slice(8));
    if (!file.startsWith(mathjax + path.sep) || !fs.existsSync(file)) {
      res.writeHead(404); res.end(); return;
    }
    res.end(fs.readFileSync(file));
  } else if (req.url === "/equation-editor.js") {
    res.end(fs.readFileSync(path.join(web, "equation-editor.js")));
  } else {
    res.setHeader("Content-Type", "text/html; charset=utf-8");
    res.end('<!doctype html><html lang="ja"><meta charset="utf-8">' +
      '<script>window.MathJax={startup:{typeset:false}};</script>' +
      '<script src="/mathjax/tex-chtml.js"></script>' +
      // The homepage supplies this adapter; exercise the real preview engine.
      '<script>window.SugaharaMath={typeset:function(node){' +
      'return MathJax.startup.promise.then(function(){return MathJax.typesetPromise([node]);}).catch(function(){});}};</script>' +
      fs.readFileSync(path.join(web, "equation-editor.fragment.html"), "utf8")
        .replace(req.url === "/?cold" ? '<script src="./equation-editor.js"></script>' : "__unused__", ""));
  }
});
(async () => {
  let browser;
  await new Promise(resolve => server.listen(0, "127.0.0.1", resolve));
  try {
    browser = await chromium.launch({headless:true,
      ...(process.env.EQNEDIT64_BROWSER_PATH ? {executablePath:process.env.EQNEDIT64_BROWSER_PATH} : {})});
    const cold = await browser.newPage();
    let releaseCancel;
    const cancelGate = new Promise(resolve => {releaseCancel = resolve;});
    await cold.route("**/cancel.js", async route => {await cancelGate; await route.continue();});
    await cold.goto("http://127.0.0.1:" + server.address().port + "/?cold");
    await cold.evaluate(async () => {
      await MathJax.startup.promise;
      const typeset = MathJax.typesetPromise;
      const gate = new Promise(resolve => {window.releasePalettePreviews = resolve;});
      MathJax.typesetPromise = function (nodes) {
        if (nodes && nodes[0].classList.contains("eqed-math-face")) {
          window.palettePreviewStarted = true;
          return gate.then(() => typeset(nodes));
        }
        return typeset(nodes);
      };
      // Capture the application's real copy event without touching the host
      // clipboard or granting browser clipboard permissions.
      window.capturedCopies = [];
      const exec = document.execCommand.bind(document);
      document.execCommand = function (command, ...args) {
        if (command !== "copy") return exec(command, ...args);
        const data = new DataTransfer();
        const event = new ClipboardEvent("copy", {clipboardData:data,cancelable:true});
        document.dispatchEvent(event);
        window.capturedCopies.push({html:data.getData("text/html"),tex:data.getData("text/plain")});
        return event.defaultPrevented;
      };
    });
    await cold.addScriptTag({url:"http://127.0.0.1:" + server.address().port + "/equation-editor.js"});
    assert(await cold.locator(".eqed-copy-office").isDisabled(), "Copy is unavailable while macros load");
    await cold.locator(".eqed-source").fill("\\bm{x}+\\cancel{x}");
    releaseCancel();
    await cold.waitForFunction(() => !document.querySelector(".eqed-copy-office").disabled && window.palettePreviewStarted);
    assert.equal(await cold.locator(".eqed-math-face mjx-container").count(),0);
    await cold.locator(".eqed-copy-office").click();
    const copies = await cold.evaluate(() => window.capturedCopies);
    assert.equal(copies.length,1);
    assert.equal(copies[0].tex,"\\bm{x}+\\cancel{x}");
    assert.match(copies[0].html,/menclose/);
    assert.match(copies[0].html,/mathvariant="bold-italic"/);
    assert.doesNotMatch(copies[0].html,/<merror|mathcolor="red"/);
    await cold.evaluate(() => window.releasePalettePreviews());
    await cold.waitForFunction(() => document.querySelectorAll(".eqed-math-face").length ===
      document.querySelectorAll(".eqed-math-face mjx-container, .eqed-preview-error .eqed-math-face").length);
    assert.equal(await cold.locator(".eqed-preview-error").count(),0);
    await cold.close();
    const failed = await browser.newPage();
    await failed.route("**/cancel.js", route => route.abort());
    await failed.goto("http://127.0.0.1:" + server.address().port + "/?cold");
    await failed.evaluate(() => MathJax.startup.promise);
    await failed.addScriptTag({url:"http://127.0.0.1:" + server.address().port + "/equation-editor.js"});
    await failed.locator(".eqed-source").fill("a^{2}");
    await failed.waitForFunction(() => document.querySelector(".eqed-copy-office").title.includes("再読み込み"));
    assert(await failed.locator(".eqed-copy-office").isDisabled());
    await failed.locator(".eqed-preview mjx-container").waitFor({timeout:10000});
    assert.doesNotMatch(await failed.locator(".eqed-preview").innerText(), /準備しています/);
    await failed.locator(".eqed-source").fill("b+1");
    await failed.waitForFunction(() => document.querySelector(".eqed-preview").textContent.includes("b"));
    await failed.close();
    const page = await browser.newPage({viewport:{width:1400,height:1000}});
    const errors = [];
    page.on("pageerror", error => errors.push(error.message));
    await page.goto("http://127.0.0.1:" + server.address().port);
    await page.waitForFunction(() => document.querySelectorAll(".eqed-math-face").length > 0 &&
      document.querySelectorAll(".eqed-math-face").length ===
      document.querySelectorAll(".eqed-math-face mjx-container, .eqed-preview-error .eqed-math-face").length,
      {}, {timeout:60000});
    await page.evaluate(() => document.fonts.ready);
    const result = await page.evaluate(() => ({
      keys:document.querySelectorAll(".eqed-key").length,
      disabled:document.querySelectorAll(".eqed-key:disabled").length,
      errors:[...document.querySelectorAll(".eqed-preview-error")].map(b=>b.title),
      sans:getComputedStyle(document.querySelector(".TEX-SS")).fontFamily,
      fraktur:getComputedStyle(document.querySelector(".TEX-FR")).fontFamily,
      strike:document.querySelectorAll(".eqed-math-face mjx-ustrike").length
    }));
    assert.equal(result.keys,291); assert.equal(result.disabled,4);
    assert.deepEqual(result.errors,[]); assert.deepEqual(errors,[]);
    // Independent MathJax semantics: an explicit prime glyph is in the outer
    // exponent, not an apostrophe creating another exponent inside it.
    const primeChecks = await page.evaluate(() => [1,2,3].map(count => {
      const tex = "{a^{2}}^{" + "\\prime ".repeat(count) + "}";
      const xml = new DOMParser().parseFromString(MathJax.tex2mml(tex), "text/xml");
      return {scripts:xml.getElementsByTagName("msup").length,
              errors:xml.getElementsByTagName("merror").length,
              primes:(xml.documentElement.textContent.match(/′/g)||[]).length};
    }));
    assert.deepEqual(primeChecks,[1,2,3].map(primes => ({scripts:2,errors:0,primes})));
    assert.match(result.sans,/MJXTEX-SS/); assert.match(result.fraktur,/MJXTEX-FR/);
    assert(result.strike > 0);
    // The actual click must put selected text in the radicand, not its index.
    await page.locator(".eqed-source").fill("x");
    await page.locator(".eqed-source").evaluate(e=>e.setSelectionRange(0,1));
    await page.locator('.eqed-key').filter({hasText:"n√□"}).click();
    assert.equal(await page.locator(".eqed-source").inputValue(),"\\sqrt[]{x}");
    for (const count of [1,2,3]) {
      await page.locator(".eqed-source").fill("x");
      await page.locator(".eqed-source").evaluate(e=>e.setSelectionRange(1,1));
      const face = ["x′", "x″", "x‴"][count-1];
      await page.getByRole('button', {name:face, exact:true}).click();
      assert.equal(await page.locator(".eqed-source").inputValue(),
                   "x^{" + "\\prime ".repeat(count) + "}");
    }
    await page.getByRole('tab', {name:"Web追加", exact:true}).click();
    for (const [base, selected, expected] of [
      ["x", true, "{x}^{\\circ}"],
      ["a+b", true, "{a+b}^{\\circ}"],
      ["a^{2}", false, "a^{2}{}^{\\circ}"],
      ["x^n_i", false, "x^n_i{}^{\\circ}"]
    ]) {
      await page.locator(".eqed-source").fill(base);
      await page.locator(".eqed-source").evaluate((e, selected) =>
        e.setSelectionRange(selected ? 0 : e.value.length, e.value.length), selected);
      await page.getByRole('button', {name:"°", exact:true}).click();
      const tex = await page.locator(".eqed-source").inputValue();
      assert.equal(tex.trim(), expected);
      const errorCount = await page.evaluate(tex => {
        const xml = new DOMParser().parseFromString(MathJax.tex2mml(tex), "text/xml");
        return xml.getElementsByTagName("merror").length;
      }, tex);
      assert.equal(errorCount, 0, tex);
    }
    console.log("PASS: cold copy, 291 keys, fonts, root body, and degree selection/script attachment");
  } finally {
    if (browser) await browser.close();
    server.close();
  }
})().catch(error => {console.error(error); process.exitCode = 1;});
