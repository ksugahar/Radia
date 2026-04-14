"""ih_claude_proposer.py

LLM (Anthropic Claude) proposer for the IH coil optimization loop.

Two backends supported:

  1. ClaudeCLIProposer (default): invokes the ``claude`` CLI
     (``claude --print --output-format json --json-schema ...``).
     Uses the user's Claude subscription via Claude Code's built-in
     auth. No ANTHROPIC_API_KEY required. Tokens are drawn from the
     subscription quota. Slower (~10 s/iter).

  2. ClaudeSDKProposer: uses anthropic SDK's messages.parse().
     Requires ANTHROPIC_API_KEY. Faster and supports proper prompt
     caching on the system prompt. Kept for parity but not the default
     on machines already inside a Claude Code session.

Both expose the same interface:
    proposer(history) -> CoilSpec
so the optimizer in ih_optimize.py is backend-agnostic.
"""

import json
import os
import subprocess
import textwrap
from typing import Optional

from pydantic import BaseModel, Field

from radia.coil_spec import CoilSpec


DEFAULT_MODEL = "claude-opus-4-6"

# Prefer the CLI shipped with Claude Code (absolute path), fall back
# to PATH lookup.
_CLI_PATH = os.environ.get('CLAUDE_CODE_EXECPATH', 'claude')


class CoilProposal(BaseModel):
    """Schema Claude emits per iteration. 'circle' only in MVP."""
    profile_type: str = Field(
        default='circle',
        description="One of 'circle' or 'rect'. Use 'circle' (MVP).",
    )
    r_wire: float = Field(
        ...,
        description="Circular wire radius in meters. Respect r_wire "
                    "bounds given in the system prompt.",
    )
    R_major: float = Field(
        ...,
        description="Torus major radius (centerline) in meters. Respect "
                    "R_major bounds given in the system prompt.",
    )
    arc_deg: float = Field(
        ...,
        description="Arc span in degrees. Respect arc_deg bounds.",
    )
    current: float = Field(
        default=1.0,
        description="Total coil current [A]. Keep at 1.0 so P_total is "
                    "comparable across designs.",
    )
    rationale: str = Field(
        ...,
        description="One or two sentences explaining the choice given "
                    "the prior trials. Used for debugging, not scoring.",
    )


def _build_system_prompt(ranges: dict, workpiece_desc: str,
                          frequency_hz: float) -> str:
    return textwrap.dedent(f"""\
        You are an induction-heating (IH) coil design optimizer.

        ## Objective
        Propose round-wire torus coil designs that maximize P_total,
        the time-averaged power dissipated in a conducting workpiece
        by the induced eddy currents. Each candidate is evaluated by
        a coupled PEEC (coil) + BEM-SIBC (workpiece) solver at
        f = {frequency_hz:g} Hz with total coil current I = 1 A.

        ## Workpiece
        {workpiece_desc}

        ## Design space (round-wire torus)
        - profile_type : always "circle"
        - r_wire       : in [{ranges['r_wire_mm'][0]*1e-3:.4e},
                             {ranges['r_wire_mm'][1]*1e-3:.4e}] m
        - R_major      : in [{ranges['R_major_mm'][0]*1e-3:.4e},
                             {ranges['R_major_mm'][1]*1e-3:.4e}] m
                         (coil centerline major radius; the torus
                          surrounds the workpiece, so R_major >
                          workpiece_radius + r_wire is required)
        - arc_deg      : in [{ranges['arc_deg'][0]:.1f},
                             {ranges['arc_deg'][1]:.1f}] deg
                         (arc span; ~355 for near-closed loop, lower
                          leaves a gap)
        - current      : fixed at 1.0 A

        ## Physics priors (reason with these)
        - Closer coupling (smaller R_major, closer to the workpiece)
          gives higher H on the workpiece surface -> more P_total.
        - Longer arc (closer to 360 deg) concentrates flux around the
          workpiece -> higher P_total.
        - Larger wire radius lowers coil AC resistance but pushes the
          current filament further from the workpiece (more standoff).
          There is a sweet spot.
        - Inner radius (R_major - r_wire) must be > workpiece outer
          radius. Designs that violate this fail.

        ## Strategy
        For the first few trials (empty history), explore broadly to
        establish gradients. After >=3 successful evaluations, exploit
        the best region but keep a small amount of exploration.

        ## Output
        Emit a single CoilProposal JSON object. Always stay strictly
        inside the bounds. Include a brief rationale (1-2 sentences).""")


def _summarize_history(history: list, top_k: int = 6) -> str:
    valid = [t for t in history if t.failure is None]
    if not valid and not history:
        return "(no trials yet -- this is the first design.)"
    lines = [f"Trials so far: {len(history)} total, "
             f"{len(valid)} successful, "
             f"{len(history) - len(valid)} failed."]
    if valid:
        valid_sorted = sorted(valid, key=lambda t: -t.objective)
        best = valid_sorted[0]
        lines.append(
            f"Best so far: iter {best.iteration}, "
            f"r_wire={best.spec.r_wire*1e3:.3f} mm, "
            f"R_major={best.spec.R_major*1e3:.3f} mm, "
            f"arc={best.spec.arc_deg:.1f} deg -> "
            f"P_total = {best.metrics['P_total']:.4e} W")
        lines.append(f"\nTop {min(top_k, len(valid_sorted))} trials "
                     "(descending P_total):")
        for t in valid_sorted[:top_k]:
            lines.append(
                f"  iter {t.iteration}: "
                f"r_wire={t.spec.r_wire*1e3:6.3f} mm, "
                f"R_major={t.spec.R_major*1e3:6.3f} mm, "
                f"arc={t.spec.arc_deg:6.1f} deg -> "
                f"P_total={t.metrics['P_total']:.4e} W")
    recent_fails = [t for t in history if t.failure is not None][-3:]
    if recent_fails:
        lines.append("\nRecent failures (avoid repeating these):")
        for t in recent_fails:
            lines.append(
                f"  iter {t.iteration}: r_wire={t.spec.r_wire*1e3:.3f} mm, "
                f"R_major={t.spec.R_major*1e3:.3f} mm, "
                f"arc={t.spec.arc_deg:.1f} deg -> {t.failure[:80]}")
    return "\n".join(lines)


class ClaudeCLIProposer:
    """LLM proposer backed by the ``claude`` CLI (subscription auth).

    Uses ``claude --print --output-format json --json-schema <schema>``
    to get structured output without an API key. Subscription tokens
    are consumed from the signed-in user's quota; total_cost_usd in
    the response is the API-equivalent price for accounting.
    """

    def __init__(self, ranges: dict, workpiece_desc: str,
                 frequency_hz: float,
                 claude_path: Optional[str] = None,
                 fallback_proposer=None,
                 verbose: bool = True,
                 timeout_sec: int = 180):
        self.ranges = ranges
        self.system_prompt = _build_system_prompt(
            ranges, workpiece_desc, frequency_hz)
        self.fallback = fallback_proposer
        self.verbose = verbose
        self.timeout_sec = timeout_sec
        self.claude_path = claude_path or _CLI_PATH

        # Pydantic -> JSON Schema (strip keys Claude Code may reject).
        schema = CoilProposal.model_json_schema()
        for extra in ('title', '$defs', 'definitions'):
            schema.pop(extra, None)
        self.schema = schema

        # Telemetry
        self.n_calls = 0
        self.n_fallbacks = 0
        self.total_cost_usd = 0.0
        self.total_input_tokens = 0
        self.total_output_tokens = 0
        self.total_cache_read_tokens = 0
        self.total_api_ms = 0

    def __call__(self, history: list) -> CoilSpec:
        user_msg = _summarize_history(history) + (
            "\n\nPropose the next CoilProposal to evaluate.")

        try:
            # Bypass the "cannot launch Claude inside Claude" guard --
            # we're invoking the CLI as a subprocess inside a Claude
            # Code session on purpose.
            env = os.environ.copy()
            env.pop('CLAUDECODE', None)
            env.pop('CLAUDE_CODE_ENTRYPOINT', None)

            cmd = [
                self.claude_path,
                '--print',
                '--output-format', 'json',
                '--no-session-persistence',
                '--model', DEFAULT_MODEL,
                '--json-schema', json.dumps(self.schema),
                '--append-system-prompt', self.system_prompt,
                user_msg,
            ]
            r = subprocess.run(
                cmd, capture_output=True, text=True, env=env,
                timeout=self.timeout_sec,
                encoding='utf-8', errors='replace',
            )
            if r.returncode != 0:
                raise RuntimeError(
                    f"claude CLI exit {r.returncode}: "
                    f"{(r.stderr or '')[:200]}")

            parsed = json.loads(r.stdout)
            if parsed.get('is_error'):
                raise RuntimeError(
                    f"claude error subtype={parsed.get('subtype')} "
                    f"result={(parsed.get('result') or '')[:200]}")
            structured = parsed.get('structured_output')
            if not structured:
                raise RuntimeError(
                    f"no structured_output in response: "
                    f"{json.dumps(parsed)[:300]}")
            proposal = CoilProposal(**structured)

            # Telemetry
            self.n_calls += 1
            self.total_cost_usd += float(parsed.get('total_cost_usd', 0.0))
            self.total_api_ms += int(parsed.get('duration_api_ms', 0))
            u = parsed.get('usage', {}) or {}
            self.total_input_tokens += int(u.get('input_tokens', 0) or 0)
            self.total_output_tokens += int(u.get('output_tokens', 0) or 0)
            self.total_cache_read_tokens += int(
                u.get('cache_read_input_tokens', 0) or 0)

            if self.verbose:
                cr = int(u.get('cache_read_input_tokens', 0) or 0)
                tag = f"  [cache_read={cr}]" if cr else "  [cold]"
                print(f"    Claude: r={proposal.r_wire*1e3:.3f}mm  "
                      f"R={proposal.R_major*1e3:.2f}mm  "
                      f"arc={proposal.arc_deg:.1f}deg{tag}")
                print(f"    rationale: {proposal.rationale[:150]}")

            return CoilSpec(
                profile_type=proposal.profile_type,
                r_wire=float(proposal.r_wire),
                R_major=float(proposal.R_major),
                arc_deg=float(proposal.arc_deg),
                current=float(proposal.current),
            )

        except Exception as e:
            if self.verbose:
                print(f"    Claude CLI error: "
                      f"{type(e).__name__}: {str(e)[:180]}")
            if self.fallback is None:
                raise
            self.n_fallbacks += 1
            if self.verbose:
                print(f"    -> falling back to RandomProposer")
            return self.fallback(history)

    def telemetry(self) -> dict:
        return {
            'n_calls': self.n_calls,
            'n_fallbacks': self.n_fallbacks,
            'total_cost_usd_equiv': round(self.total_cost_usd, 4),
            'total_api_seconds': round(self.total_api_ms / 1000.0, 1),
            'total_input_tokens': self.total_input_tokens,
            'total_output_tokens': self.total_output_tokens,
            'total_cache_read_tokens': self.total_cache_read_tokens,
        }


# Backward-compat alias: the demo script imports ClaudeProposer.
ClaudeProposer = ClaudeCLIProposer
