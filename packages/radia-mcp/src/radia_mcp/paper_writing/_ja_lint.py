"""Inline copy of 8 JA-lint helpers from mcp-server-document.grant_writing.

Copied 2026-05-26 so radia_mcp.paper_writing (public PyPI) is
self-contained.  Canonical source: mcp-server-document.grant_writing.tools.
Do NOT add new dependencies on grant_writing-specific plans modules —
the 8 functions below are intentionally standalone.
"""
from __future__ import annotations

import pathlib
import re

# Local _shared.hedges -- needed by grant_writing_lint_bedrock which
# reports hedge counts as a sub-axis of bedrock analysis.  The original
# grant_writing/tools.py imported scan_hedges from ..\\_shared.hedges;
# after the 2026-05-26 migration the cross-package _shared.hedges was
# inline-copied into paper_writing/_shared/, so the import path becomes
# the LOCAL ._shared.hedges below.  Without this restore lint_bedrock
# raises NameError: name '_scan_hedges' is not defined.
from ._shared.hedges import scan_hedges as _scan_hedges

# Real content starts below.
_HERE = pathlib.Path(__file__).resolve().parent
KNOWLEDGE = _HERE
_REDUNDANCY = [
    (r"することができ(る|ない)", r"でき\1"),
    (r"ということができる", "といえる"),
    (r"に関する", "の"),
    (r"に関して(は)?", ""),
    (r"においては", "では"),
    (r"における", "の"),
    (r"ではないだろうか", "であろう"),
    (r"にほかならない", "である"),
    (r"することによって", "して"),
    (r"することとなる", "する"),
    (r"のような形で", "で"),
    (r"の中で", "で"),
    (r"にあたり", "するとき"),
    (r"ということである", "である"),
    (r"高い精度を持つ", "精度が高い"),
    (r"大きな影響を与える", "大きく影響する"),
    (r"の形で存在する", "として存在する"),
    (r"を目的として行う", "のために"),
    (r"を行う", "する"),
    (r"を行なう", "する"),
    (r"増加することになる", "増加する"),
    (r"〜という[^、。]", ""),
    (r"することとする", "する"),
    (r"〜に対しての", "への"),
    (r"の場合において", "で"),
]
_MISUSE = [
    (r"よろしかったでしょうか", "「よろしいでしょうか」に"),
    (r"のほう(?=を|が|は)", "「のほう」はぼかし。削除するか具体語に"),
    (r"じゃないですか(?=[、。])", "「のですが」に (同意誘導は弱い)"),
    (r"^なので、", "文頭「なので」は口語。「したがって」「そのため」に"),
    (r"(?<=[^で])ないです(?=[。、])", "口頭語。書面では「ありません」"),
    (r"やむ[お]えない", "正: 「やむを得ない」"),
    (r"(わたし|僕|自分|個人)的に(は)?", "「私としては」「個人としては」"),
    (r"こんにち[わ]", "正: 「こんにちは」"),
    (r"お(連絡|報告|通知|返事|案内|質問|説明)", "漢語には「ご」"),
    (r"違(く[てかっ]|かった)", "正: 「違って」「違った」"),
    (r"(きもい|うざい|きしょい|ぶっちゃけ)", "俗語。書面では使用禁止"),
    (r"みたいな[。、]", "「という感じだ」「と言った」"),
    (r"なにげに", "正: 「何気なく」「それとなく」"),
    (r"とんでもあ(りません|ございません)", "正: 「とんでもない (ことです)」"),
    (r"(咲い|降っ|流れ|吹い)ております", "主体が物なら「〜ています」"),
]
_DEFAULT_ACRONYM_WHITELIST = {
    # 機関・組織
    "IEEE", "ACM", "JSPS", "JST", "NEDO", "MEXT", "AMED", "RIKEN", "NIMS",
    "NIH", "NSF", "DARPA", "DOE", "CERN", "ESRF", "MIT", "TU",
    # ファイル・標準
    "PDF", "HTML", "XML", "JSON", "CSV", "URL", "URI", "API", "SDK",
    "README",
    "GUI", "CLI", "IDE", "OS",
    # ライセンス
    "BSD", "GPL", "LGPL", "MPL",
    # 一般技術
    "CPU", "GPU", "RAM", "ROM", "SSD", "HDD", "IO", "UI", "UX",
    "AI", "ML", "DL", "NLP", "CV", "CAD", "CAM", "CAE",
    # Unicode / encoding
    "UTF", "ASCII", "ISO",
    # 記法・場所表記 (スケジュール)
    # D1, D2 ... D31 (day-1..day-31) は schedule markers — whitelist
    # I, II, III, IV 等 ローマ数字 は section markers — whitelist
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
}
_KANJI = r"一-鿿"
_NOTATION_PAIRS = [
    # 形式名詞 (ひらがな推奨)。前に漢字が無く、直後に助詞が来るときのみ候補化。
    # これで 仕事/事業/物事/物質/研究所/場所/行為/開発時 等の熟語を誤検出しない。
    (rf"(?<![{_KANJI}])事(?=[はがもをにで、。])", "こと"),
    (rf"(?<![{_KANJI}])時(?=[はがもをにで、。])", "とき"),
    (rf"(?<![{_KANJI}])物(?=[はがもをにで、。])", "もの"),
    (rf"(?<![{_KANJI}])所(?=[はがもをにで、。])", "ところ"),
    (rf"(?<![{_KANJI}])為(?=[に、])", "ため"),
    # 形式動詞 (補助動詞はひらがな推奨)
    (r"(?<=[てでにを])居る", "いる"),
    (r"(?<=[てでにを])在る", "ある"),
    (r"(?<=[てでにを])頂く", "いただく"),
    (r"(?<=て)下さい", "ください"),
    # 送り仮名の揺れ
    (r"行なう", "行う"),   # 行う に統一
    (r"表わす", "表す"),   # 表す に統一
    (r"現わす", "現す"),
    (r"表われる", "表れる"),
    # カタカナ語末長音 (JIS Z 8301:2005 / MS スタイル 2008 改訂以降)
    # - 2 音節以下は長音なし (データ が慣例、データー は誤り)
    # - 3 音節以上は長音あり (コンピューター, サーバー, ユーザー)
    (r"データー", "データ"),
    (r"コンピュータ(?!ー)", "コンピューター"),
    (r"サーバ(?!ー)", "サーバー"),
    (r"ユーザ(?!ー)", "ユーザー"),
    # 英字括弧スタイル (半角統一推奨)
    # 数値 + 単位の間隔 (半角スペース + \, 推奨)
]

def _load_skill() -> str:
    return (KNOWLEDGE / "skill.md").read_text(encoding="utf-8")

def grant_writing_lint_bedrock(text: str) -> dict:
    """木下 10 原則 + 本多 + 知的 による和文技術文章 bedrock 診断。

    1. 緩衝表現 (「と思われる」「と考えられる」)
    2. 受身連用 (3 文連続で受身末尾)
    3. 逆茂木型文 (「、」4 個以上 + 80 字以上)
    4. 「の」連打 3 個以上
    5. 一文中「は」2 個以上
    6. 文末モノトーン (3 文連続同じ文末)
    7. 意見の事実断定 (主観修飾語 + 断定)
    8. 二重否定 (「ないではない」「なくはない」等)

    すべての writing-family MCP で共通利用可能。
    """
    issues: list[dict] = []

    # Split sentences for several checks
    sentences = [s.strip() for s in re.split(r"[。．]", text) if s.strip()]

    # 1. 緩衝表現 — _HEDGE_PATTERNS で集中管理 (count_weak_expressions と共通)
    hedge_total, hedge_by_pat = _scan_hedges(text)
    if hedge_total:
        issues.append({
            "rule": "hedging_expressions",
            "severity": "HIGH",
            "count": hedge_total,
            "examples": list(hedge_by_pat.keys())[:5],
            "fix": "断定表現に置換 (「と考える」等の基本形)",
            "source": "木下 p.95",
        })

    # 2. 受身連用 (windowed)
    passive_ends = [bool(re.search(r"(される|られる|された|られた)$", s))
                    for s in sentences]
    consecutive_passive = 0
    max_consec = 0
    for p in passive_ends:
        if p:
            consecutive_passive += 1
            max_consec = max(max_consec, consecutive_passive)
        else:
            consecutive_passive = 0
    if max_consec >= 3:
        issues.append({
            "rule": "passive_voice_chain",
            "severity": "MODERATE",
            "max_consecutive": max_consec,
            "fix": "能動態に書き換え (主語を明示)",
            "source": "木下 p.140",
        })

    # 3. 逆茂木型
    long_sentences = [s for s in sentences
                      if len(s) >= 80 and s.count("、") >= 4]
    if long_sentences:
        issues.append({
            "rule": "gyakumogi_sentence",
            "severity": "MODERATE",
            "count": len(long_sentences),
            "examples": [s[:60] + "..." for s in long_sentences[:3]],
            "fix": "長い前置修飾節は 2 個まで。文を分割。",
            "source": "木下 p.87",
        })

    # 4. 「の」連打 3 個
    no_runs = re.findall(r"の[^の、。\s]{0,8}の[^の、。\s]{0,8}の", text)
    if no_runs:
        issues.append({
            "rule": "no_particle_triple",
            "severity": "LOW",
            "count": len(no_runs),
            "examples": no_runs[:5],
            "fix": "「の」は連続 2 個まで。別助詞へ。",
            "source": "知的 p.46",
        })

    # 5. 一文中「は」2 個
    wa_heavy = []
    for s in sentences:
        # Remove quoted content
        s2 = re.sub(r"「[^」]*」", "", s)
        if s2.count("は") >= 2:
            wa_heavy.append(s)
    if wa_heavy:
        issues.append({
            "rule": "double_wa_in_sentence",
            "severity": "LOW",
            "count": len(wa_heavy),
            "examples": [s[:60] + "..." for s in wa_heavy[:3]],
            "fix": "1 文に「は」は 1 つまで。主題を 1 つに絞る。",
            "source": "知的 p.48",
        })

    # 6. 文末モノトーン
    if len(sentences) >= 3:
        endings = [s[-3:] if len(s) >= 3 else s for s in sentences]
        mono_count = 0
        for i in range(len(endings) - 2):
            if endings[i] == endings[i+1] == endings[i+2]:
                mono_count += 1
        if mono_count >= 1:
            issues.append({
                "rule": "monotone_sentence_ending",
                "severity": "LOW",
                "mono_3_in_a_row_count": mono_count,
                "fix": "文末多様化 (「である」「挙げられる」「報告されている」等を混在)",
                "source": "知的 p.36",
            })

    # 7. 意見の事実断定
    opinion_fact = re.findall(
        r"(便利|すぐれ|優秀|重要|有用|効果的|画期的)[な]?[^。]*である[。\.]",
        text
    )
    if opinion_fact:
        issues.append({
            "rule": "opinion_stated_as_fact",
            "severity": "HIGH",
            "count": len(opinion_fact),
            "examples": [s[:60] + "..." for s in opinion_fact[:3]],
            "fix": "意見と事実を峻別。数値で裏付けるか「〜と考える」の形へ。",
            "source": "木下 p.117",
        })

    # 8. 二重否定
    double_neg = re.findall(r"ない[^。]{0,12}(こと|わけ|の)(は|が)?ない",
                            text)
    if double_neg:
        issues.append({
            "rule": "double_negative",
            "severity": "LOW",
            "count": len(double_neg),
            "fix": "肯定文で書き直す (「有効でないとは言えない」→「有効である」)",
            "source": "中学生 p.190",
        })

    return {
        "total_sentences": len(sentences),
        "issue_count": len(issues),
        "issues": issues,
        "bedrock_principles": (
            "木下 10 原則 (目標規定文 / 重点先行 / パラグラフ 1 主張 / "
            "事実意見峻別 / 逆茂木排除 / 修飾近接 / 受身征伐 / "
            "緩衝表現削除 / 誤解排除 / 主述直結) + 本多テン二大原則"
        ),
    }

def grant_writing_suggest_redundancy_fixes(text: str) -> dict:
    """和文の典型的冗長表現 25 パターンを検出し置換候補を示す。

    『日本語の作文技術』+ 『知的な科学・技術文章の書き方』+ 『問題な日本語』
    の定番 lint。
    """
    suggestions: list[dict] = []
    for pat, replacement in _REDUNDANCY:
        for m in re.finditer(pat, text):
            is_delete = (replacement == "")
            suggestions.append({
                "pattern": pat,
                "match": m.group(0),
                "position": m.start(),
                "action": "delete" if is_delete else "replace",
                "suggested": replacement if not is_delete else None,
                "note": ("削除 (文脈確認必須)" if is_delete
                         else f"'{replacement}' に置換"),
            })
    # De-duplicate patterns (first occurrence only per pattern)
    seen = set()
    dedup = []
    for s in suggestions:
        if s["pattern"] not in seen:
            seen.add(s["pattern"])
            dedup.append(s)
    return {
        "total_matches": len(suggestions),
        "unique_patterns_hit": len(dedup),
        "top_fixes": dedup[:20],
        "note": (
            "置換候補は文脈依存。機械的に全置換せず、1 件ずつ判定。"
            "action='delete' は特に文脈確認必須 (削るとねじれる場合がある)。"
        ),
        "source": "本多『日本語の作文技術』+ 『知的な科学・技術文章の書き方』",
    }

def grant_writing_check_misuse_japanese(text: str) -> dict:
    """『問題な日本語』由来の現代誤用 15 パターン検出。

    申請書で低頻度だが、混入すると「校正が甘い」印象で減点対象。
    """
    hits: list[dict] = []
    for pat, hint in _MISUSE:
        for m in re.finditer(pat, text):
            hits.append({
                "pattern": pat,
                "match": m.group(0),
                "position": m.start(),
                "fix": hint,
            })
    return {
        "total_matches": len(hits),
        "issues": hits[:30],
        "source": "北原保雄『問題な日本語』",
    }

def grant_writing_check_kanji_ratio(text: str,
                                     min_ratio: float = 0.18,
                                     max_ratio: float = 0.40,
                                     ) -> dict:
    """漢字比率の偏りを検出。本多『日本語の作文技術』第四章に基づく。

    推奨 range: 0.25-0.35。< 0.18 = 幼稚/冗長、> 0.40 = 「官庁漢語」感。

    漢字・ひらがな・カタカナの unicode 範囲で単純分類。
    """
    kanji = sum(1 for c in text if "一" <= c <= "鿿")
    hira = sum(1 for c in text if "぀" <= c <= "ゟ")
    kata = sum(1 for c in text if "゠" <= c <= "ヿ")
    total = kanji + hira + kata
    if total == 0:
        return {"error": "no Japanese characters found"}
    ratio = kanji / total
    if ratio < min_ratio:
        status = "hiragana_heavy"
        hint = "漢字を増やして締まりを持たせる"
    elif ratio > max_ratio:
        status = "kanji_heavy"
        hint = "形式名詞 (こと/とき/もの) を漢字にしない。官庁漢語を削減。"
    else:
        status = "balanced"
        hint = "良好 (推奨範囲)"
    return {
        "total_jp_chars": total,
        "kanji": kanji,
        "hiragana": hira,
        "katakana": kata,
        "kanji_ratio": round(ratio, 3),
        "target_range": f"{min_ratio:.2f} - {max_ratio:.2f}",
        "status": status,
        "hint": hint,
        "source": "本多『日本語の作文技術』第四章",
    }

def grant_writing_check_subject_predicate_distance(text: str,
                                                     max_chars: int = 40,
                                                     ) -> dict:
    """主述の直結原則 (本多 p.22): 主語と述語の間の距離が遠い文を検出。

    判定法: 「は、」または「が、」 のパターンで **明示的な主題マーカー** を
    特定し、その位置から文末までの距離を計測。40 字超で warning、60 字超
    で fatal。

    設計上の制約 (2026-04-23):
      - 読点の伴わない「は/が」(例:「Xは重要である」) は分析対象外。
        主述距離が問題になるのは主題の後に長い修飾節が挟まる文のみ。
      - この制約で「については」「にとっては」等の副助詞「は」を主題と
        誤認する false positive を大きく削減している。
      - 結果として再現率は下がるが hint tool として精度を優先。
    """
    sentences = [s.strip() for s in re.split(r"[。．]", text) if s.strip()]
    violations: list[dict] = []
    skipped_no_topic_comma = 0
    for s in sentences:
        m = re.search(r"(は|が)、", s)
        if not m:
            skipped_no_topic_comma += 1
            continue
        subj_pos = m.start()
        distance = len(s) - subj_pos
        if distance > max_chars:
            violations.append({
                "sentence_head": s[:50] + ("..." if len(s) > 50 else ""),
                "subject_marker_pos": subj_pos,
                "distance_to_end": distance,
            })
    return {
        "total_sentences": len(sentences),
        "analyzed_sentences": len(sentences) - skipped_no_topic_comma,
        "skipped_no_topic_comma": skipped_no_topic_comma,
        "violation_count": len(violations),
        "violations": violations[:10],
        "target": f"距離 <= {max_chars} 字",
        "fatal": "距離 > 60 字 — 読解困難",
        "fix": "長い修飾節を前置せず、主語と述語を近接させる。文を分割。",
        "hint": (
            "hint として使う。「は、」「が、」を含まない文 (単純文・"
            "副助詞「は」含む文) は解析対象外。全件を 0 に潰そうとしない。"
        ),
        "source": "本多『日本語の作文技術』第一章",
    }

def grant_writing_find_undefined_acronyms(
    text: str,
    whitelist: str = "",
    min_len: int = 2,
    max_len: int = 6,
    context_window: int = 200,
) -> dict:
    """Latin 略語の初出で定義 (〜 or 〜の略) が近くに無いものを検出。

    検出ロジック:
      1. 2-6 字の大文字 Latin + 数字 acronym を抽出 (単語境界)
      2. 初出位置周辺 (±context_window 文字) で
         `ACRONYM (expansion)` or `expansion (ACRONYM)` パターンを検索
      3. どちらも見つからない場合 = 未定義
      4. whitelist (既知の略語: IEEE / JSPS / PDF 等) は除外

    Args:
        text: 検査対象テキスト (tex コマンドは事前に除去推奨)
        whitelist: 追加で除外する acronym をカンマ区切りで指定
            (例: "ESIM,MSFEM,POD" のようにドメイン固有用語の後回しも可)
        min_len, max_len: 検出対象の文字長 (既定 2-6)
        context_window: 初出周辺の検索範囲 (既定 200 文字)

    Returns:
        dict with:
          - total_acronyms: 検出された unique acronym 数
          - undefined: 未定義 acronym のリスト [{acronym, first_pos, context}]
          - defined: 定義済み acronym のリスト
          - whitelist_skipped: whitelist 除外された件数
    """
    user_whitelist = {w.strip().upper() for w in whitelist.split(",") if w.strip()}
    combined_whitelist = _DEFAULT_ACRONYM_WHITELIST | user_whitelist

    # Also whitelist D1-D99 (day-N schedule markers) and ローマ数字 I-XX
    schedule_marker_re = re.compile(r"^D\d{1,2}$")
    roman_re = re.compile(r"^[IVX]{1,5}$")

    # Python's ``\b`` treats Japanese letters as word characters, so ``MCPで``
    # has no boundary and was not counted. Use an ASCII identifier boundary.
    acronym_re = re.compile(
        r"(?<![A-Za-z0-9])([A-Z][A-Z0-9]{"
        + str(min_len - 1)
        + ","
        + str(max_len - 1)
        + r"})(?![A-Za-z0-9])"
    )
    first_seen: dict[str, int] = {}
    for m in acronym_re.finditer(text):
        a = m.group(1)
        # JP-MARs and similarly hyphenated project names are one proper noun,
        # not an unexplained acronym followed by an unrelated word.
        if m.end() + 1 < len(text) and text[m.end()] in "-–—" and text[m.end() + 1].isalpha():
            continue
        if a not in first_seen:
            first_seen[a] = m.start()

    undefined: list[dict] = []
    defined: list[str] = []
    skipped = 0
    for acronym, pos in sorted(first_seen.items(), key=lambda kv: kv[1]):
        if acronym in combined_whitelist:
            skipped += 1
            continue
        if schedule_marker_re.match(acronym) or roman_re.match(acronym):
            skipped += 1
            continue

        # Check (expansion) right after: "ACRONYM (expansion ...)"
        after_window = text[pos:pos + context_window]
        after_re = re.compile(
            r"\b" + re.escape(acronym) + r"\s*[\(（]([^\)）]+)[\)）]"
        )
        has_after = bool(after_re.search(after_window))

        # Check "expansion (ACRONYM)" appearing anywhere up to first-use
        before_all = text[:pos + len(acronym) + 2]
        before_re = re.compile(
            r"[\(（]\s*" + re.escape(acronym) + r"\s*[\)）]"
        )
        gloss_re = re.compile(
            r"[\(（][^()（）]{3,100}[;；,，]\s*"
            + re.escape(acronym)
            + r"\s*[\)）]"
        )
        has_before = bool(before_re.search(before_all) or gloss_re.search(before_all))

        if has_after or has_before:
            defined.append(acronym)
        else:
            context_start = max(0, pos - 30)
            context_end = min(len(text), pos + len(acronym) + 40)
            snippet = text[context_start:context_end].replace("\n", " ")
            undefined.append({
                "acronym": acronym,
                "first_pos": pos,
                "context": snippet,
            })

    return {
        "total_acronyms": len(first_seen),
        "undefined_count": len(undefined),
        "defined_count": len(defined),
        "whitelist_skipped": skipped,
        "undefined": undefined[:40],
        "defined_sample": defined[:20],
        "hint": (
            "未定義 acronym は **選択的に** 初出で `ACRONYM (展開)` を足す。"
            "ドメイン固有語 (ESIM, MSFEM 等) を reviewer が知っている前提なら "
            "定義不要。全件を 0 に潰そうとすると冗長で逆に読みにくくなる。"
            "whitelist 引数で skip 対象をカンマ区切り指定可。"
            "**使用頻度を含めた総合判定は grant_writing_acronym_usage_audit** を使う。"
        ),
        "source": "木下『理科系の作文技術』原則 9「誤解できないように書け」"
                  " (p.131)。 **hint として使用** — violation count ではない。",
    }

def grant_writing_acronym_usage_audit(
    text: str,
    whitelist: str = "",
    min_uses_for_abbrev: int = 3,
    min_len: int = 2,
    max_len: int = 6,
) -> dict:
    """略語の使用頻度と初出形式を監査し、3 段階の推奨を返す。

    **ルール** (grant-writing 独自、2026-04-22 採用):

    1. **初出は略語でなく full form で書く**。`(ACRONYM)` を括弧で添える
       形が標準。例: ✅「Effective Surface Impedance Method (ESIM)」、
       ❌「ESIM (Effective Surface Impedance Method)」。
    2. **2 回目以降は略語 OK**。
    3. **登場回数が少ない (既定 3 回未満) 時は そもそも略語を導入しない**。
       full form で通した方が読みやすい。

    Args:
        text: 検査対象テキスト。
        whitelist: skip する略語 (IEEE / JSPS 等は組み込み whitelist で除外済)。
        min_uses_for_abbrev: この回数未満なら略語を導入しない推奨 (既定 3)。
        min_len, max_len: 略語の文字長範囲 (既定 2-6)。

    Returns:
        dict:
          - total_acronyms
          - findings: list of {acronym, count, first_form ("expansion"/"acronym"),
                                verdict, recommendation}
          - summary_counts: {drop_abbrev, fix_first_use, ok}
          - hint / source
    """
    user_whitelist = {w.strip().upper() for w in whitelist.split(",") if w.strip()}
    combined = _DEFAULT_ACRONYM_WHITELIST | user_whitelist
    schedule_re = re.compile(r"^D\d{1,2}$")
    roman_re = re.compile(r"^[IVX]{1,5}$")
    acronym_re = re.compile(
        r"(?<![A-Za-z0-9])([A-Z][A-Z0-9]{"
        + str(min_len - 1)
        + ","
        + str(max_len - 1)
        + r"})(?![A-Za-z0-9])"
    )
    # Count occurrences and first position per acronym
    occurrences: dict[str, list[int]] = {}
    for m in acronym_re.finditer(text):
        a = m.group(1)
        if a in combined or schedule_re.match(a) or roman_re.match(a):
            continue
        if m.end() + 1 < len(text) and text[m.end()] in "-–—" and text[m.end() + 1].isalpha():
            continue
        occurrences.setdefault(a, []).append(m.start())

    findings = []
    counts_by_verdict = {"drop_abbrev": 0, "fix_first_use": 0, "ok": 0}

    for acronym, positions in sorted(occurrences.items(),
                                       key=lambda kv: kv[1][0]):
        count = len(positions)
        first_pos = positions[0]

        # Determine first-use form: is `expansion (ACRONYM)` present BEFORE
        # or AT the first standalone acronym occurrence?
        before_slice = text[:first_pos + len(acronym) + 2]
        # Match `expansion (ACRONYM)` — any parenthesised ACRONYM before first_pos
        before_paren_re = re.compile(
            r"[\(（]\s*" + re.escape(acronym) + r"\s*[\)）]"
        )
        gloss_re = re.compile(
            r"[\(（][^()（）]{3,100}[;；,，]\s*"
            + re.escape(acronym)
            + r"\s*[\)）]"
        )
        has_paren_before = any(
            m.start() < first_pos for m in before_paren_re.finditer(before_slice)
        ) or bool(gloss_re.search(before_slice))

        first_form = "expansion" if has_paren_before else "acronym"

        # Apply the three-tier rule
        if count < min_uses_for_abbrev:
            verdict = "drop_abbrev"
            rec = (f"登場 {count} 回のみ。略語を導入せず full form で通す方が "
                   f"読みやすい。")
        elif first_form == "acronym":
            verdict = "fix_first_use"
            rec = (f"登場 {count} 回。初出は full form「... ({acronym})」の形に "
                   f"書き換え、2 回目以降を略語に。")
        else:
            verdict = "ok"
            rec = (f"登場 {count} 回、初出 full form。適切。")

        counts_by_verdict[verdict] += 1

        # Small context around first use
        cs = max(0, first_pos - 30)
        ce = min(len(text), first_pos + len(acronym) + 40)
        findings.append({
            "acronym": acronym,
            "count": count,
            "first_pos": first_pos,
            "first_form": first_form,
            "verdict": verdict,
            "recommendation": rec,
            "context": text[cs:ce].replace("\n", " "),
        })

    return {
        "total_acronyms": len(occurrences),
        "min_uses_for_abbrev": min_uses_for_abbrev,
        "findings": findings[:60],
        "summary_counts": counts_by_verdict,
        "hint": (
            "grant-writing 独自ルール: **初出 full form → 2 回目以降 略語**。"
            f"**{min_uses_for_abbrev} 回未満の登場なら略語を導入しない**。"
            "hint として使う — 審査員の既知略語 (IEEE / JSPS 等) は "
            "組み込み whitelist + user whitelist で skip 可。"
        ),
        "source": "grant-writing rule (2026-04-22): full-form-first, "
                  "abbreviate only if it earns its keep (>=3 occurrences)。",
    }

def grant_writing_check_notation_variants(text: str) -> dict:
    """同一テキスト内で同じ概念が複数の表記で書かれていないかを検出。

    検出対象:
      1. 形式名詞の漢字/ひらがな混在 (事/こと, 時/とき, 物/もの, 所/ところ, 為/ため)
      2. 補助動詞の漢字/ひらがな混在 (居る/いる, 頂く/いただく, 下さい/ください)
      3. 送り仮名の揺れ (行う/行なう, 表す/表わす)
      4. カタカナ末尾長音 (データ/データー)
      5. 半角・全角の括弧混在
      6. 数値+単位の空白揺れ (50 kHz / 50kHz / 50\\,kHz)
      7. 大小文字揺れ (固有名詞、Apache-stylen fuzzy match)
      8. ハイフン揺れ (Go-Tech / GoTech / Go Tech)

    Args:
        text: 検査対象テキスト

    Returns:
        dict: 検出された表記ゆれのペアと出現回数。多数派を推奨として示す。
    """
    findings: list[dict] = []

    # (1)-(4): 既知表記ペアの検出
    for pat_src, preferred in _NOTATION_PAIRS:
        try:
            matches = list(re.finditer(pat_src, text))
            if matches:
                # Count preferred form too
                pref_count = len(re.findall(re.escape(preferred), text))
                findings.append({
                    "type": "kana_kanji_variant",
                    "pattern": pat_src,
                    "preferred": preferred,
                    "variant_count": len(matches),
                    "preferred_count": pref_count,
                    "action": (f"'{preferred}' に統一推奨"
                                if pref_count >= len(matches)
                                else f"'{preferred}' に統一 (現在少数派)"),
                    "examples": [m.group(0) for m in matches[:3]],
                })
        except re.error:
            continue

    # (5): 半角/全角括弧混在
    han_open = text.count("(")
    zen_open = text.count("（")
    if han_open > 0 and zen_open > 0:
        findings.append({
            "type": "paren_mixed",
            "half_width_open": han_open,
            "full_width_open": zen_open,
            "action": (
                f"半角 '(', ')' {han_open} 個 + 全角 '（', '）' {zen_open} 個 混在。"
                "JSPS は本文半角・法令表記のみ全角 が標準。統一推奨。"
            ),
        })

    # (6): 数値 + 単位の空白ゆれ
    unit_variants: dict[str, dict] = {}
    for unit in ("kHz", "mm", "km", "MHz", "GHz", "Hz", "V", "A"):
        no_space = len(re.findall(r"\d+(?:\.\d+)?" + unit + r"\b", text))
        half_space = len(re.findall(r"\d+(?:\.\d+)? " + unit + r"\b", text))
        thinspace = len(re.findall(r"\d+(?:\.\d+)?\\,\s*" + unit + r"\b", text))
        total = no_space + half_space + thinspace
        if total >= 2 and len(
            [c for c in (no_space, half_space, thinspace) if c > 0]
        ) >= 2:
            unit_variants[unit] = {
                "no_space": no_space,
                "half_space": half_space,
                "thinspace_backslash_comma": thinspace,
            }
    if unit_variants:
        findings.append({
            "type": "unit_spacing_variants",
            "variants": unit_variants,
            "action": (
                "LaTeX 組版では `\\,` (thin space) で数値と単位を区切るのが標準。"
                " 統一推奨。"
            ),
        })

    # (7): 大小文字揺れ (固有名詞・識別子)
    # 例: "NGSolve" vs "ngsolve" vs "NGsolve"
    # 普通の英単語 (This / from / when 等) を拾わないため、
    # group 内に「大文字 2 個以上 or 数字を含む」識別子様 token が
    # 1 つでもある場合のみ flag する。
    def _is_identifier_like(tok: str) -> bool:
        upper = sum(1 for c in tok if c.isupper())
        has_digit = any(c.isdigit() for c in tok)
        return upper >= 2 or has_digit

    token_re = re.compile(r"\b([A-Za-z][A-Za-z0-9]{3,})\b")
    tokens: dict[str, set[str]] = {}
    for m in token_re.finditer(text):
        t = m.group(1)
        tokens.setdefault(t.lower(), set()).add(t)
    case_variants = [
        {
            "lowercase": k,
            "variants": sorted(v),
            "counts": {form: len(re.findall(r"\b" + re.escape(form) + r"\b", text))
                       for form in v},
        }
        for k, v in tokens.items()
        if len(v) >= 2 and any(_is_identifier_like(f) for f in v)
    ]
    if case_variants:
        findings.append({
            "type": "case_variants",
            "count": len(case_variants),
            "examples": case_variants[:10],
            "action": "大小文字の統一。多数派 or 公式表記に揃える。",
        })

    # (8): ハイフン揺れ (Go-Tech / GoTech / Go Tech)
    # ハイフン形の各 token について、同じ文書内に
    #   - 連結形 (GoTech)
    #   - スペース形 (Go Tech)
    # のいずれかが存在すれば表記ゆれとして flag する。
    hyphen_tokens: dict[str, set[tuple[str, str]]] = {}
    hyph_re = re.compile(r"\b([A-Za-z][A-Za-z0-9]{1,10})-([A-Za-z][A-Za-z0-9]{1,10})\b")
    for m in hyph_re.finditer(text):
        a, b = m.group(1), m.group(2)
        key = (a + b).lower()
        hyphen_tokens.setdefault(key, set()).add((a, b))
    hyphen_variants = []
    for key, parts_set in hyphen_tokens.items():
        a, b = min(parts_set)
        hyphen_forms = sorted({f"{pa}-{pb}" for pa, pb in parts_set})
        concat = a + b
        spaced = f"{a} {b}"
        concat_count = len(re.findall(r"\b" + re.escape(concat) + r"\b", text))
        spaced_count = len(re.findall(
            r"\b" + re.escape(a) + r" " + re.escape(b) + r"\b", text
        ))
        if concat_count == 0 and spaced_count == 0:
            continue
        hyphen_variants.append({
            "base": key,
            "hyphen_forms": hyphen_forms,
            "concat_form": concat,
            "spaced_form": spaced,
            "hyphen_count": sum(
                len(re.findall(r"\b" + re.escape(f) + r"\b", text))
                for f in hyphen_forms
            ),
            "concat_count": concat_count,
            "spaced_count": spaced_count,
        })
    if hyphen_variants:
        findings.append({
            "type": "hyphen_variants",
            "examples": hyphen_variants[:5],
            "action": "ハイフン・スペース・連結形を統一する。",
        })

    return {
        "total_findings": len(findings),
        "findings": findings,
        "hint": (
            "表記ゆれは審査員の集中力を削ぐ。1 文書内で同一概念は 1 表記に統一する。"
            "形式名詞はひらがな、固有名詞は公式表記、数値+単位は LaTeX の `\\,` "
            "(thin space) を標準とする。"
        ),
        "source": "本多『日本語の作文技術』 + 共通工学日本語スタイルガイド",
    }
