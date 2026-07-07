# grant-writing

Use this skill for grant proposal drafts, recommendation letters, and final
pre-submission checks.  It is intentionally parallel to the other document
servers:

- `paper-writing`: journal / digest manuscripts
- `figure`: publication-quality figures
- `grant-writing`: proposal logic, prose, feasibility, budget, social impact
- `presentation`: slides and oral-delivery artifacts

## Core Rules

1. Start from the reviewer's question: what public or technical problem is
   solved, why now, and why this applicant can finish it.
2. Convert "we will consider" into verifiable deliverables: a dataset, demo,
   prototype, public repository, report, measurement, workshop, or presentation.
3. Keep social implementation concrete: user, field, PoC, evaluation method,
   and handoff path.
4. Budget text must explain why each cost is necessary for the proposal's
   verification loop, not for generic productivity.  Do not make the requested
   budget artificially small: apply close to the program ceiling when the plan
   genuinely needs it, and make that look natural by showing careful itemized
   calculations.
5. Recommendation letters should say only what the recommender can plausibly
   attest: applicant ability, institutional fit, support, feasibility, and
   significance.

## Budget Policy

予算はほぼ上限いっぱいで申請してよい。重要なのは、上限近くであることを
遠慮して隠すのではなく、研究計画、PoC、計算資源、AI agent 運用、評価、
発表・社会実装に必要な経費として精査済みに見せることである。

予算欄では、単価 x 数量 x 月数/回数、年度配分、見積根拠、検証ループとの
対応を具体的に書く。読み手が「この計画なら上限いっぱいになっても不思議
ではない」と判断できる粒度まで積算する。

## Useful Tools

- `grant_writing_usage()`
- `grant_writing_health_report(text_or_path, program="generic")`
- `grant_writing_section_presence(text, program="generic")`
- `grant_writing_kddi_digital_check(text)`
- `grant_writing_kddi_power_electronics_focus_check(text)`
- `grant_writing_budget_alignment_check(text)`
- `grant_writing_analyze_sentences(text)`
- `grant_writing_count_weak_expressions(text)`
- `grant_writing_lint_bedrock(text)`
- `grant_writing_recommendation_letter_template(program="kddi_digital")`

For KDDI Foundation Digital Innovation / social implementation proposals,
use `program="kddi_digital"` so the report checks social issue, digital use,
PoC, schedule, budget, feasibility, and implementation outcomes.

For the power-electronics-board CAE-AI proposal, also use
`grant_writing_kddi_power_electronics_focus_check(text)`.  It checks that the
main subject is the power-electronics-board circuit / electromagnetic /
thermal CAE-AI environment, while companies below 1000 employees are kept as
the first users and implementation field.  It also warns when commercial CAE
is framed as an adversarial replacement instead of a powerful but hard-to-access
tool category for which AI/MCP provides an entry point.
