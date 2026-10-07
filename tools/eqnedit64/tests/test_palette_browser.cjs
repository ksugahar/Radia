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
      (req.url === '/?delayed' ? '' : '<script src="/mathjax/tex-chtml.js"></script>') +
      // The homepage supplies this adapter; exercise the real preview engine.
      '<script>window.SugaharaMath={typeset:function(node){' +
      'return MathJax.startup.promise.then(function(){return MathJax.typesetPromise([node]);}).catch(function(){});}};</script>' +
      fs.readFileSync(path.join(web, "equation-editor.fragment.html"), "utf8")
        .replace((req.url === "/?cold" || req.url === "/?delayed") ? '<script src="./equation-editor.js"></script>' : "__unused__", ""));
  }
});
(async () => {
  let browser;
  await new Promise(resolve => server.listen(0, "127.0.0.1", resolve));
  try {
    browser = await chromium.launch({headless:true,
      ...(process.env.EQNEDIT64_BROWSER_PATH ? {executablePath:process.env.EQNEDIT64_BROWSER_PATH} : {})});
    for (const config of [true, false]) {
      const delayed = await browser.newPage();
      const delayedErrors = [];
      delayed.on("pageerror", error => delayedErrors.push(error.message));
      await delayed.goto("http://127.0.0.1:" + server.address().port + "/?delayed");
      await delayed.evaluate(config => {
        window.MathJax = config ? {startup:{typeset:false}} : undefined;
      }, config);
      await delayed.addScriptTag({url:"http://127.0.0.1:" + server.address().port + "/equation-editor.js"});
      assert.equal(await delayed.locator(".eqed-key").count(), 291);
      assert(await delayed.locator(".eqed-copy-office").isDisabled());
      await delayed.evaluate(() => {
        const originalNow = Date.now;
        Date.now = () => originalNow() + 61000;
      });
      await delayed.waitForFunction(() => document.querySelector(".eqed-copy-office").title.includes("読み込みを待っています"));
      await delayed.addScriptTag({url:"http://127.0.0.1:" + server.address().port + "/mathjax/tex-chtml.js"});
      await delayed.waitForFunction(() => !document.querySelector(".eqed-copy-office").disabled);
      assert.deepEqual(delayedErrors, []);
      await delayed.locator(".eqed-source").fill("\\foo x");
      await delayed.locator(".eqed-copy-office").click();
      await delayed.waitForFunction(() => document.querySelector(".eqed-status").textContent.includes("MathMLに変換できません"));
      await delayed.evaluate(() => {
        window.copyAttempts = 0;
        document.execCommand = function (command) {
          if (command !== "copy") return false;
          window.copyAttempts++;
          const data = new DataTransfer();
          const event = new ClipboardEvent("copy", {clipboardData:data, cancelable:true});
          document.dispatchEvent(event);
          return event.defaultPrevented;
        };
      });
      await delayed.locator(".eqed-source").fill("\\color{red}{x}");
      await delayed.locator(".eqed-copy-office").click();
      await delayed.waitForFunction(() => window.copyAttempts === 1);
      await delayed.waitForFunction(() => document.querySelector(".eqed-status").textContent.includes("コピーしました"));
      await delayed.close();
    }
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
    // Exercise the production copy handler. This verifies generated payloads,
    // not PowerPoint import; host clipboard stays untouched.
    for (const [tex, expected] of [
      ["y=\\int_a^b x^2\\,dx", {limits:true,naryBody:"x2dx"}],
      ["\\begin{aligned}y&=\\int_a^b x^2\\,dx\\\\z&=1\\end{aligned}", {limits:true,rows:2}],
      ["\\begin{aligned}y=\\int_a^b x^2\\,dx\\\\z=1\\end{aligned}", {limits:true,rows:2}],
      ["f=\\begin{cases}\\int_a^b x\\,dx & x>0\\end{cases}", {limits:true,matrix:true}],
      ["\\begin{pmatrix}\\int_a^b x\\,dx & 0\\\\0 & 1\\end{pmatrix}", {limits:true,matrix:true}],
      ["\\color{red}{x}", {red:true}],
      ["\\textcolor{red}{x}+\\textcolor{blue}{y}", {red:true,blue:true}],
      ["\\frac{x_1^2}{\\sqrt[3]{y}}", {fraction:true,radical:true}],
      ["\\sqrt{x+y}", {radical:true,radicalBody:"x+y"}],
      ["\\binom{n}{k}", {fraction:true,noBar:true}],
      ["\\sum_{i=1}^{n}x_i", {nary:true}]
    ]) {
      await cold.locator(".eqed-source").fill(tex);
      const before = await cold.evaluate(() => window.capturedCopies.length);
      await cold.locator(".eqed-copy-office").click();
      const payload = await cold.evaluate(() => window.capturedCopies.at(-1));
      assert.equal(await cold.evaluate(() => window.capturedCopies.length), before+1,tex);
      assert.equal(payload.tex,tex);
      assert.match(payload.html,/<!--\[if gte msEquation 12\]>/);
      const structure = await cold.evaluate(html => {
        const omml = html.match(/<!--\[if gte msEquation 12\]>([\s\S]*?)<!\[endif\]-->/)[1];
        const xml = new DOMParser().parseFromString(omml,"application/xml");
        const nodes = name => [...xml.getElementsByTagNameNS("*",name)];
        return {errors:nodes("parsererror").length,
          sub:nodes("sub").map(n=>n.textContent),sup:nodes("sup").map(n=>n.textContent),
          rows:nodes("eqArr").map(n=>[...n.children].filter(c=>c.localName==="e").length),
          matrix:nodes("m").length,fraction:nodes("f").length,radical:nodes("rad").length,
          nary:nodes("nary").length,
          naryBodies:nodes("nary").map(n=>[...n.children].find(c=>c.localName==="e").textContent.replace(/\s/g,"")),
          radicalBodies:nodes("rad").map(n=>[...n.children].find(c=>c.localName==="e").textContent),
          noBar:nodes("type").some(n=>n.getAttribute("m:val")==="noBar"),
          styles:nodes("span").map(n=>n.getAttribute("style")||"").join(";")};
      },payload.html);
      assert.equal(structure.errors,0,tex);
      if(expected.limits){assert(structure.sub.includes("a"),tex);assert(structure.sup.includes("b"),tex);}
      if(expected.rows)assert.deepEqual(structure.rows,[expected.rows],tex);
      if(expected.matrix)assert(structure.matrix>0,tex);
      if(expected.fraction)assert(structure.fraction>0,tex);
      if(expected.radical)assert(structure.radical>0,tex);
      if(expected.radicalBody)assert(structure.radicalBodies.includes(expected.radicalBody),tex);
      if(expected.noBar)assert(structure.noBar,tex);
      if(expected.nary)assert(structure.nary>0,tex);
      if(expected.naryBody)assert(structure.naryBodies.includes(expected.naryBody),tex);
      if(expected.red)assert.match(structure.styles,/color:red/);
      if(expected.blue)assert.match(structure.styles,/color:blue/);
    }
    await cold.locator(".eqed-source").fill("\\smash{x}");
    const beforeRejectedCopy = await cold.evaluate(() => window.capturedCopies.length);
    await cold.locator(".eqed-copy-office").click();
    assert.equal(await cold.evaluate(() => window.capturedCopies.length),beforeRejectedCopy);
    assert.match(await cold.locator(".eqed-status").textContent(),/Office形式に変換できません/);
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
    // Compare native RGB constants with the actual independent MathJax package.
    const colorHeader = fs.readFileSync(path.join(web, "../src/named_colors.h"), "utf8");
    const colorPairs = [...colorHeader.matchAll(/\{"([A-Za-z]+)", "(#[0-9A-F]{6})"\}/g)]
      .map(match => [match[1], match[2]]);
    assert.equal(colorPairs.length, 216);
    colorPairs.push([" red ", "#FF0000"], [" Red ", "#FF0000"]);
    const colors = await page.evaluate(pairs => pairs.map(([name, expected]) => {
      const xml = new DOMParser().parseFromString(MathJax.tex2mml("\\textcolor{"+name+"}{x}"), "text/xml");
      const styled = xml.querySelector("[mathcolor]");
      const canvas = document.createElement("canvas");
      const context = canvas.getContext("2d");
      context.fillStyle = styled ? styled.getAttribute("mathcolor") : "#000000";
      return {name, expected:expected.toLowerCase(), actual:context.fillStyle.toLowerCase(),
              errors:xml.getElementsByTagName("merror").length};
    }), colorPairs);
    for (const color of colors) {
      assert.equal(color.errors, 0, color.name);
      assert.equal(color.actual, color.expected, color.name);
    }
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
