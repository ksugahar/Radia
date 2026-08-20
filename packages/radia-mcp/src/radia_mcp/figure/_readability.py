"""Refuse to write a figure nobody can read.

`lab_savefig` already checked one thing -- the on-page font size. Everything
else that makes a figure unreadable was left to the eye, and the eye is not
reliably present. Three of them showed up in one afternoon on one deck:

- panels flattened into strips, because the figure was authored at the aspect
  of the slot it was pasted into rather than the aspect the data needs;
- an annotation sitting on top of the curves it was annotating;
- that same annotation running off the bottom of the axes.

None of those is a rendering fault, so no paste-scale or dpi check finds them.
They are found here, before the file exists, and they raise -- a figure that is
merely warned about gets shipped.
"""

from __future__ import annotations

# A panel narrower than this fraction of its width has had its vertical
# structure squeezed out.  0.28 is below the flattest panel anyone has
# defended; the one that prompted this was 0.31 and was called unreadable.
MIN_AXES_ASPECT = 0.28


def _renderer(fig):
    fig.canvas.draw()
    return fig.canvas.get_renderer()


def _squashed(fig) -> list[str]:
    out = []
    for i, ax in enumerate(fig.axes, 1):
        box = ax.get_position()
        w = box.width * fig.get_figwidth()
        h = box.height * fig.get_figheight()
        if w <= 0:
            continue
        aspect = h / w
        if aspect < MIN_AXES_ASPECT:
            out.append(
                "axes %d is %.2f tall for its width (floor %.2f): the panel is "
                "flattened into a strip, so differences between curves cannot "
                "be read" % (i, aspect, MIN_AXES_ASPECT))
    return out


def _escaping_text(fig, renderer) -> list[str]:
    out = []
    fbox = fig.bbox
    for ax in fig.axes:
        for t in list(ax.texts):
            if not t.get_text().strip():
                continue
            try:
                bb = t.get_window_extent(renderer)
            except Exception:
                continue
            if (bb.x0 < fbox.x0 - 1 or bb.y0 < fbox.y0 - 1
                    or bb.x1 > fbox.x1 + 1 or bb.y1 > fbox.y1 + 1):
                out.append("text %r runs outside the figure" % t.get_text()[:28])
    return out


def _text_on_data(fig, renderer) -> list[str]:
    """An annotation sitting on the curve it annotates."""
    out = []
    for ax in fig.axes:
        pts = []
        for line in ax.lines:
            try:
                xy = line.get_xydata()
            except Exception:
                continue
            if xy is None or len(xy) == 0 or len(xy) > 4000:
                continue
            pts.extend(ax.transData.transform(xy))
        if not pts:
            continue
        for t in list(ax.texts):
            if not t.get_text().strip():
                continue
            try:
                bb = t.get_window_extent(renderer)
            except Exception:
                continue
            hits = sum(1 for x, y in pts
                       if bb.x0 <= x <= bb.x1 and bb.y0 <= y <= bb.y1)
            if hits:
                out.append(
                    "text %r covers %d plotted point(s)"
                    % (t.get_text()[:28], hits))
    return out


def figure_readability_problems(fig) -> list[str]:
    """Everything about this figure that would make it unreadable."""
    renderer = _renderer(fig)
    return _squashed(fig) + _escaping_text(fig, renderer) + _text_on_data(fig, renderer)


def enforce_readable(fig, *, allow_unreadable: bool = False) -> list[str]:
    """Raise unless the figure can be read.

    `allow_unreadable=True` is the caller explicitly choosing to write one
    anyway -- it must be passed at the call site, so the decision is visible
    where the figure is made rather than buried in a default.
    """
    problems = figure_readability_problems(fig)
    if problems and not allow_unreadable:
        raise ValueError(
            "refusing to write an unreadable figure:\n  - "
            + "\n  - ".join(problems)
            + "\n\nFix the figure, or pass allow_unreadable=True at this call "
              "site to record that you decided otherwise."
        )
    return problems
