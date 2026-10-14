"""Plan B T11: Reproducibility / Open Science 6 軸チェック
(paper_writing_reproducibility_open_science_check).

近年 IEEE / Nature / Science / PNAS / PLoS 等の主要 journal は
**Reproducibility & Open Science** の宣言文を必須化または強く推奨している。
書き忘れによる desk reject や camera-ready 段階での add-back 要求が頻発。

検出する 6 軸:
  1. **code_availability** — コード公開声明 (GitHub/Zenodo URL等)
  2. **data_availability** — データ公開声明 (Zenodo/Figshare/repo URL等)
  3. **preregistration** — 事前登録 (OSF/AsPredicted/clinicaltrials.gov)
  4. **methods_replicability** — version / dependency / hardware spec
  5. **random_seed_disclosure** — seed 値の明示 (機械学習 / 確率実験)
  6. **competing_interests** — COI 宣言 (Conflicts / Funding 等)

設計方針:
- score は 6 軸の coverage 率 (0-10)。
- PRIMARY OUTPUT は per_axis findings + missing_required + journal hints。
- LLM-based ではなく URL pattern + keyword detection の高速 lint。
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# 軸別 patterns
# ------------------------------------------------------------------
# 1. Code availability
_CODE_AVAIL_PATTERNS = [
    r"(?:code|source\s*code|software|implementation)\s+"
    r"(?:is\s+)?(?:available|released|public|open)",
    r"(?:our|the)\s+(?:code|implementation)\s+(?:is|will be)",
    r"github\.com/[\w\-]+/[\w\-\.]+",
    r"gitlab\.com/[\w\-]+",
    r"zenodo\.org/record/\d+",
    r"コード(?:は)?(?:公開|提供|配布)",
    r"ソースコード(?:を)?(?:公開|提供)",
    r"(?:code|software)\s+availability",
]

# 2. Data availability
_DATA_AVAIL_PATTERNS = [
    r"data\s+(?:is|are)?\s*(?:available|deposited|released)",
    r"data\s+availability\s+statement",
    r"figshare\.com/(?:articles|s)/",
    r"zenodo\.org/record/\d+",
    r"dryad\.[a-z]+\.[a-z]+",
    r"openneuro\.org",
    r"データ(?:は)?(?:公開|提供|オープン化)",
    r"データセット(?:を)?(?:公開|提供)",
    r"raw\s+data",
    r"supplementary\s+data",
]

# 3. Pre-registration
_PREREG_PATTERNS = [
    r"pre[\-\s]?regist",
    r"事前登録",
    r"osf\.io/(?:[\w\-]+)",
    r"aspredicted\.org/",
    r"clinicaltrials\.gov",
    r"isrctn\.com",
    r"protocol\s+(?:was|is)\s+(?:registered|preregistered)",
]

# 4. Methods replicability (versions, dependencies, hardware)
_METHODS_REPL_PATTERNS = [
    # software versions
    r"\b(?:Python|MATLAB|Simulink|PyTorch|TensorFlow|NumPy|SciPy|NGSolve|"
    r"Netgen|Radia|COMSOL|ANSYS|JMAG|FEMM)\s+v?\d+\.\d+(?:\.\d+)?",
    r"\b[A-Za-z][A-Za-z0-9_-]*\s+(?:version\s+|v)\d+\.\d+(?:\.\d+)?",
    r"version\s+\d+\.\d+",
    r"バージョン\s*\d+\.\d+",
    # OS / hardware
    r"(?:Ubuntu|Debian|CentOS|Windows|macOS|Linux)\s+\d+",
    r"(?:NVIDIA|AMD|Intel)\s+[A-Z]\d+",
    r"GPU\s*[:：]?\s*[A-Z]\d+|GeForce|Tesla|A100|V100|H100|RTX",
    # dependencies
    r"requirements\.txt",
    r"environment\.ya?ml",
    r"Dockerfile",
    r"docker\s+image",
]

# 5. Random seed
_SEED_PATTERNS = [
    r"random[_\s]?seed[s]?",
    r"seed\s*[=:：]\s*\d+",
    r"random[_\s]?state\s*[=:：]\s*\d+",
    r"乱数(?:の)?(?:種|シード)",
    r"再現性(?:の|のため)",
    r"reproducibility[\-\s]?seed",
]

# 6. Competing interests / COI / Funding
_COI_PATTERNS = [
    r"(?:competing|conflict[s]?\s+of)\s+interest",
    r"declaration\s+of\s+(?:competing|conflict)",
    r"the\s+authors?\s+declare",
    r"利益相反",
    r"\bCOI\b",
    r"funding[\s:：]",
    r"acknowledg(?:ment|e?ment)s?",
    r"謝辞",
    r"助成金",
    r"科研費",
    r"grant\s+number",
]


def _check_axis(text: str, patterns: list[str]) -> dict:
    hits: list[str] = []
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            hits.append(m.group(0))
            if len(hits) >= 8:
                break
    return {
        "hit_count": len(hits),
        "examples": hits[:5],
        "covered": len(hits) >= 1,
    }


def paper_writing_reproducibility_open_science_check(
    text: str, journal_required: str = ""
) -> dict:
    """Reproducibility & Open Science の 6 軸を一括診断。

    Args:
        text: 論文全文 (LaTeX or markdown)
        journal_required: カンマ区切りで journal 必須項目を指定
                          (例: "code,data,coi" — 不在ならば critical 扱い)

    Returns:
        dict: score / per_axis / missing_required / comments / source
    """
    if text is None:
        text = ""

    axes = {
        "code_availability": _check_axis(text, _CODE_AVAIL_PATTERNS),
        "data_availability": _check_axis(text, _DATA_AVAIL_PATTERNS),
        "preregistration": _check_axis(text, _PREREG_PATTERNS),
        "methods_replicability": _check_axis(text, _METHODS_REPL_PATTERNS),
        "random_seed_disclosure": _check_axis(text, _SEED_PATTERNS),
        "competing_interests": _check_axis(text, _COI_PATTERNS),
    }

    n_covered = sum(1 for a in axes.values() if a["covered"])
    score = round(n_covered / 6 * 10, 1)

    # Journal-required check
    required_set = {
        s.strip().lower() for s in journal_required.split(",") if s.strip()
    }
    axis_alias = {
        "code": "code_availability",
        "data": "data_availability",
        "prereg": "preregistration",
        "preregistration": "preregistration",
        "methods": "methods_replicability",
        "seed": "random_seed_disclosure",
        "coi": "competing_interests",
        "funding": "competing_interests",
    }
    missing_required: list[str] = []
    for req in required_set:
        axis_name = axis_alias.get(req)
        if axis_name and not axes[axis_name]["covered"]:
            missing_required.append(axis_name)

    if missing_required:
        score = max(0.0, score - 2 * len(missing_required))

    # ---- comments ----
    comments: list[str] = []

    if n_covered == 0:
        comments.append(
            "Reproducibility & Open Science 関連の記述が一切検出されず。"
            "近年の主要 journal (Nature/Science/IEEE/PLoS) は最低 1 軸 "
            "(典型的には COI と code/data availability) を必須化。"
            "desk reject 回避のため最低限の宣言文を追加推奨。"
        )
    else:
        comments.append(
            f"6 軸中 {n_covered} 軸で記述検出。"
            f"covered: {[k for k, v in axes.items() if v['covered']]}。"
        )

    not_covered = [k for k, v in axes.items() if not v["covered"]]

    if "code_availability" in not_covered:
        comments.append(
            "Code availability 宣言なし。GitHub/Zenodo URL を 1 文 + "
            "「Code is available at https://...」のような宣言を Methods 末尾 "
            "or 専用セクションに配置。実装系論文では事実上必須。"
        )

    if "data_availability" in not_covered:
        comments.append(
            "Data availability 宣言なし。データ公開なら Zenodo/Figshare/Dryad の "
            "DOI を、非公開ならその理由 (機密 / IRB 制約) を明示する。"
            "Nature 系では「Data Availability Statement」セクション必須。"
        )

    if "competing_interests" in not_covered:
        comments.append(
            "COI / Funding 宣言なし。投稿前の最終チェックで必ず追加。"
            "「The authors declare no competing interests.」+ "
            "「This work was supported by grant number XXX.」。"
        )

    if (
        "random_seed_disclosure" in not_covered
        and re.search(r"machine learning|deep learning|neural network|"
                       r"機械学習|ニューラル|RNN|CNN|Transformer|LLM",
                       text, re.IGNORECASE)
    ):
        comments.append(
            "機械学習関連の記述があるが random seed の disclosure なし。"
            "ML 論文では再現性のため seed 値明示が標準 "
            "(NeurIPS/ICML reproducibility checklist 必須項目)。"
        )

    if missing_required:
        comments.append(
            f"⚠ Journal 必須項目で未対応: {missing_required}。"
            "投稿前に必ず追加。submit 即 desk reject の典型理由。"
        )

    if score >= 8:
        comments.append(
            f"Reproducibility & Open Science score {score}/10 — "
            "主要 journal の必須要件をほぼ充足。"
        )

    comments = comments[:8]

    hint = (
        "本ツールは regex/keyword baseline で 6 軸 detect。"
        "code/data availability は実装系・実験系で必須化が進む一方、"
        "純理論系では preregistration/seed は適用外。"
        "journal_required 引数で対象 journal の必須項目を渡すと "
        "missing_required で critical 扱いされる。"
        "**重要**: URL の存在検証は別途 (本ツールは記述の有無のみ)。"
    )

    return {
        "score": score,
        "score_max": 10,
        "axes_covered": n_covered,
        "axes_total": 6,
        "per_axis": axes,
        "missing_required": missing_required,
        "comments": comments,
        "hint": hint,
        "source": (
            "TOP Guidelines (Transparency and Openness Promotion, "
            "Nosek et al. 2015 Science 348:1422) + "
            "NeurIPS Reproducibility Checklist + "
            "Nature Author Tutorials (Reporting Standards 2024) + "
            "ICMJE recommendations (COI/funding disclosure)。"
            "(paper-writing 2026-04 採用)。"
        ),
    }
