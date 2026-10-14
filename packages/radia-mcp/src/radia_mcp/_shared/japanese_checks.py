"""Shared Japanese checks with caller-owned policy and preprocessing.

These private pure functions are not MCP tools. Wrappers retain public
signatures/docstrings and supply the appropriate rules or hedge scanner.
"""

import re


def _lint_bedrock(text: str, scan_hedges) -> dict:
    issues: list[dict] = []

    # Split sentences for several checks
    sentences = [s.strip() for s in re.split(r"[。．]", text) if s.strip()]

    # 1. 緩衝表現 — _HEDGE_PATTERNS で集中管理 (count_weak_expressions と共通)
    hedge_total, hedge_by_pat = scan_hedges(text)
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

def _suggest_redundancy_fixes(text: str, rules) -> dict:
    suggestions: list[dict] = []
    for pat, replacement in rules:
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

def _check_misuse_japanese(text: str, rules) -> dict:
    hits: list[dict] = []
    for pat, hint in rules:
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
