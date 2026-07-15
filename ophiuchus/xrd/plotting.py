from __future__ import annotations

from pathlib import Path

from .models import CandidateScore, Peak, XrdPattern


def select_dashboard_plot_scores(top_scores: list[CandidateScore]) -> list[tuple[str, CandidateScore]]:
    labels = ["目标模拟峰", "可能杂质 1", "可能杂质 2", "可能杂质 3"]
    return list(zip(labels, top_scores[:4]))


def write_xrd_plot(
    path: str | Path,
    pattern: XrdPattern,
    experimental_peaks: list[Peak],
    top_scores: list[CandidateScore],
    max_candidates: int = 4,
) -> str:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    return _write_matplotlib_plot(out, pattern, experimental_peaks, top_scores, max_candidates)


def _write_matplotlib_plot(
    out: Path,
    pattern: XrdPattern,
    experimental_peaks: list[Peak],
    top_scores: list[CandidateScore],
    max_candidates: int,
) -> str:
    import matplotlib.pyplot as plt

    figure = build_xrd_figure(pattern, experimental_peaks, top_scores, max_candidates)
    figure.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(figure)
    return str(out)


def build_xrd_figure(
    pattern: XrdPattern,
    experimental_peaks: list[Peak],
    top_scores: list[CandidateScore],
    max_candidates: int = 4,
):
    import warnings

    warnings.filterwarnings("ignore", message=".*deprecated.*")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs = [point.two_theta for point in pattern.points]
    ys_raw = [point.intensity for point in pattern.points]
    if not xs or not ys_raw:
        fig, ax = plt.subplots(figsize=(11, 5), dpi=300)
        return fig

    max_y = max(ys_raw) or 1.0
    ys = [value / max_y * 100.0 for value in ys_raw]
    x_min, x_max = min(xs), max(xs)
    selected = select_dashboard_plot_scores(top_scores)[:max_candidates]

    from ophiuchus.theme import configure_matplotlib_fonts

    configure_matplotlib_fonts()
    plt.rcParams["axes.linewidth"] = 1.0
    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(11, 6.2),
        dpi=300,
        sharex=True,
        gridspec_kw={
            "height_ratios": [1.45, 1.0],
            "hspace": 0.035,
        },
    )

    exp_ax = axes[0]
    exp_ax.plot(xs, ys, color="black", linewidth=0.9)
    exp_ax.text(
        0.012,
        0.82,
        "Experimental",
        transform=exp_ax.transAxes,
        ha="left",
        va="center",
        fontsize=10,
        color="black",
    )
    exp_ax.set_ylim(-4, 108)
    exp_ax.set_yticks([])
    exp_ax.tick_params(axis="x", labelbottom=False)

    colors = ["red", "blue", "green", "#ff9900"]
    target_reference_peaks: list[Peak] = []
    guide_peaks: list[tuple[Peak, str]] = []
    ref_ax = axes[1]
    legend_handles = []
    if selected:
        target = selected[0][1]
        target_peaks = [peak for peak in target.candidate.theory_peaks if x_min <= peak.two_theta <= x_max]
        simulated = target.candidate.simulated_pattern
        profile_x = tuple(getattr(simulated, "profile_two_theta_deg", ()) or ())
        profile_y = tuple(getattr(simulated, "profile_normalized_intensity", ()) or ())
        profile_rows = [(x, y) for x, y in zip(profile_x, profile_y) if x_min <= x <= x_max]
        if profile_rows:
            target_xs = [row[0] for row in profile_rows]
            target_ys = [row[1] for row in profile_rows]
        else:
            target_xs, target_ys = _broaden_peaks(target_peaks, x_min, x_max, step=0.02, sigma=0.055)
        target_reference_peaks = _strongest_peaks(target_peaks, 8)
        target_label = f"Target simulated: {target.candidate.formula_pretty}"
        ref_ax.plot(target_xs, target_ys, color="red", linewidth=1.0, label=target_label)
        ref_ax.vlines(
            [peak.two_theta for peak in target_peaks],
            [0.0 for _peak in target_peaks],
            [peak.intensity for peak in target_peaks],
            color="red",
            linewidth=0.45,
            alpha=0.22,
            label="_ophi_target_pattern",
        )
        for peak in target_reference_peaks:
            ref_ax.plot(peak.two_theta, peak.intensity, "o", color="red", markersize=3)
        _label_peaks(ref_ax, target_reference_peaks, color="red", y_offset=4.0, y_limit=112)
        guide_peaks.extend((peak, "red") for peak in target_reference_peaks)
        legend_handles.append(plt.Line2D([0], [0], color="red", linewidth=1.2, label=target_label))

    impurity_offsets = [-18.0, -36.0, -54.0]
    for idx, (_role, score) in enumerate(selected[1:4]):
        color = colors[idx + 1]
        offset = impurity_offsets[idx]
        peaks = [peak for peak in score.candidate.theory_peaks if x_min <= peak.two_theta <= x_max]
        ref_ax.hlines(offset, x_min, x_max, color=color, linewidth=1.0, alpha=0.95)
        ref_ax.text(x_min + 0.4, offset + 3.0, score.candidate.formula_pretty, ha="left", va="bottom", fontsize=9, color=color)
        for peak in peaks[:160]:
            ref_ax.vlines(peak.two_theta, offset, offset + 12.0 * peak.intensity / 100.0, color=color, linewidth=0.8, alpha=0.95)
        strong = _strongest_peaks(peaks, 3)
        for peak in strong:
            y = offset + 12.0 * peak.intensity / 100.0
            ref_ax.plot(peak.two_theta, y, "o", color=color, markersize=2.8)
        guide_peaks.extend((peak, color) for peak in strong)
        legend_handles.append(plt.Line2D([0], [0], color=color, linewidth=1.0, label=f"Impurity {idx + 1}: {score.candidate.formula_pretty}"))

    ref_ax.text(
        0.012,
        0.88,
        "Reference simulated",
        transform=ref_ax.transAxes,
        ha="left",
        va="center",
        fontsize=10,
        color="red",
    )
    ref_ax.set_ylim(-64, 112)
    ref_ax.set_yticks([])
    if legend_handles:
        ref_ax.legend(handles=legend_handles, loc="upper right", frameon=True, fontsize=8)

    for ax in axes:
        ax.set_xlim(x_min, x_max)
        ax.spines["right"].set_visible(True)
        ax.spines["top"].set_visible(True)

    ref_ax.set_xlabel(r"2$\theta$ (degree)", fontsize=11)
    fig.text(0.025, 0.5, "Intensity / a.u.", rotation=90, va="center", ha="center", fontsize=11)
    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.095, top=0.985, hspace=0.04)

    fig.canvas.draw()
    y0 = axes[-1].get_position().y0
    y1 = axes[0].get_position().y1
    for peak, color in guide_peaks:
        x_display = axes[-1].transData.transform((peak.two_theta, 0))[0]
        x_fig = fig.transFigure.inverted().transform((x_display, 0))[0]
        fig.add_artist(
            plt.Line2D(
                [x_fig, x_fig],
                [y0, y1],
                transform=fig.transFigure,
                color=color,
                linestyle="--",
                linewidth=0.75,
                alpha=0.46 if color == "red" else 0.32,
                zorder=0.5,
            )
        )

    return fig


def _strongest_peaks(peaks: list[Peak], top_n: int) -> list[Peak]:
    return sorted(sorted(peaks, key=lambda peak: peak.intensity, reverse=True)[:top_n], key=lambda peak: peak.two_theta)


def _broaden_peaks(peaks: list[Peak], x_min: float, x_max: float, step: float = 0.02, sigma: float = 0.055) -> tuple[list[float], list[float]]:
    import math

    count = max(2, int((x_max - x_min) / step) + 1)
    xs = [x_min + i * step for i in range(count)]
    ys = [0.0 for _ in xs]
    window = 4.0 * sigma
    for peak in peaks:
        if peak.two_theta < x_min or peak.two_theta > x_max:
            continue
        start = max(0, int((peak.two_theta - window - x_min) / step))
        end = min(count - 1, int((peak.two_theta + window - x_min) / step) + 1)
        for idx in range(start, end + 1):
            dx = xs[idx] - peak.two_theta
            ys[idx] += peak.intensity * math.exp(-(dx * dx) / (2.0 * sigma * sigma))
    max_y = max(ys) if ys else 0.0
    if max_y > 0:
        ys = [y / max_y * 100.0 for y in ys]
    return xs, ys


def _label_peaks(
    ax,
    peaks: list[Peak],
    color: str,
    y_offset: float,
    y_limit: float,
    base_offset: float = 0.0,
    scale: float = 100.0,
) -> None:
    placed: list[float] = []
    for peak in peaks:
        if any(abs(peak.two_theta - old) < 0.45 for old in placed):
            continue
        placed.append(peak.two_theta)
        y = base_offset + scale * peak.intensity / 100.0
        ax.text(
            peak.two_theta,
            min(y + y_offset, base_offset + y_limit),
            f"{peak.two_theta:.2f}",
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=7,
            color=color,
        )


def _write_pillow_plot(
    out: Path,
    pattern: XrdPattern,
    experimental_peaks: list[Peak],
    top_scores: list[CandidateScore],
    max_candidates: int,
) -> str:
    from PIL import Image, ImageDraw, ImageFont

    width, height = 1800, 980
    margin_l, margin_r, margin_t, margin_b = 95, 55, 45, 70
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = _load_font(24)
    small_font = _load_font(18)
    label_font = _load_font(30)
    xs = [point.two_theta for point in pattern.points]
    raw_ys = [point.intensity for point in pattern.points]
    if not xs or not raw_ys:
        img.save(out)
        return str(out)
    x_min, x_max = min(xs), max(xs)
    max_y = max(raw_ys) or 1.0
    ys = [value / max_y * 100.0 for value in raw_ys]
    plot_h = height - margin_t - margin_b
    plot_w = width - margin_l - margin_r
    selected = select_dashboard_plot_scores(top_scores)[:max_candidates]
    y_top = 135.0
    y_bottom = min(-8.0, 6 - max(0, len(selected) - 2) * 10 - 5)

    def px(x: float) -> int:
        return int(margin_l + (x - x_min) / max(x_max - x_min, 1e-9) * plot_w)

    def py(y: float) -> int:
        return int(margin_t + (y_top - y) / max(y_top - y_bottom, 1e-9) * plot_h)

    draw.rectangle((margin_l, margin_t, width - margin_r, height - margin_b), outline="#111111")
    points = [(px(x), py(y * 0.72 + 55)) for x, y in zip(xs, ys)]
    if len(points) >= 2:
        draw.line(points, fill="black", width=2)

    colors = ["red", "blue", "green", "#ff9900"]
    legend: list[tuple[str, str]] = [("Experimental", "black")]
    for idx, (_role, score) in enumerate(selected):
        color = colors[idx % len(colors)]
        peaks = [peak for peak in score.candidate.theory_peaks if x_min <= peak.two_theta <= x_max]
        if not peaks:
            continue
        if idx == 0:
            offset = 18.0
            label = f"Target simulated: {score.candidate.formula_pretty}"
            draw.line((px(x_min), py(offset), px(x_max), py(offset)), fill=color, width=1)
            for peak in peaks[:120]:
                y1 = offset + peak.intensity / 100.0 * 45.0
                draw.line((px(peak.two_theta), py(offset), px(peak.two_theta), py(y1)), fill=color, width=2)
            for peak in _strongest_peaks(peaks, 8):
                y1 = offset + peak.intensity / 100.0 * 45.0
                x = px(peak.two_theta)
                _draw_dashed_vertical(draw, x, margin_t, height - margin_b, "#ff9a9a", width=1)
                draw.ellipse((x - 3, py(y1) - 3, x + 3, py(y1) + 3), fill=color)
                _draw_rotated_text(img, f"{peak.two_theta:.2f}", (x - 10, py(y1) - 76), color, small_font)
            legend.append((label, color))
        else:
            offset = 6.0 - (idx - 1) * 10.0
            label = score.candidate.formula_pretty
            draw.line((px(x_min), py(offset), px(x_max), py(offset)), fill=color, width=2)
            draw.text((px(x_min) + 8, py(offset) - 20), label, fill=color, font=font)
            for peak in peaks[:120]:
                y1 = offset + peak.intensity / 100.0 * 7.0
                draw.line((px(peak.two_theta), py(offset), px(peak.two_theta), py(y1)), fill=color, width=2)
            for peak in _strongest_peaks(peaks, 4):
                y1 = offset + peak.intensity / 100.0 * 7.0
                x = px(peak.two_theta)
                _draw_dashed_vertical(draw, x, margin_t, height - margin_b, _soft_line_color(color), width=1)
                draw.ellipse((x - 3, py(y1) - 3, x + 3, py(y1) + 3), fill=color)
                _draw_rotated_text(img, f"{peak.two_theta:.2f}", (x - 10, py(y1) - 66), color, small_font)
            legend.append((label, color))

    legend_x = width - margin_r - 300
    legend_y = margin_t + 15
    draw.rectangle((legend_x - 12, legend_y - 10, width - margin_r - 10, legend_y + 24 * len(legend) + 8), outline="#d0d0d0", fill="white")
    for i, (label, color) in enumerate(legend):
        y = legend_y + i * 24
        draw.line((legend_x, y + 8, legend_x + 48, y + 8), fill=color, width=3)
        draw.text((legend_x + 62, y), label, fill="black", font=font)

    draw.text((width // 2 - 95, height - 46), "2theta (degree)", fill="#111111", font=label_font)
    y_label = Image.new("RGBA", (260, 48), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_label)
    y_draw.text((0, 0), "Intensity / a.u.", fill="#111111", font=label_font)
    rotated = y_label.rotate(90, expand=True)
    img.paste(rotated, (18, margin_t + plot_h // 2 - rotated.height // 2), rotated)
    img.save(out)
    return str(out)


def _load_font(size: int):
    from PIL import ImageFont

    candidates = [
        r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\segoeui.ttf",
        r"C:\Windows\Fonts\msyh.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except Exception:
            continue
    return ImageFont.load_default()


def _draw_dashed_vertical(draw, x: int, y0: int, y1: int, color: str, width: int = 1, dash: int = 12, gap: int = 8) -> None:
    y = y0
    while y < y1:
        draw.line((x, y, x, min(y + dash, y1)), fill=color, width=width)
        y += dash + gap


def _draw_rotated_text(img, text: str, xy: tuple[int, int], fill: str, font) -> None:
    from PIL import Image, ImageDraw

    bbox = font.getbbox(text)
    width = max(1, bbox[2] - bbox[0] + 8)
    height = max(1, bbox[3] - bbox[1] + 8)
    patch = Image.new("RGBA", (width, height), (255, 255, 255, 0))
    draw = ImageDraw.Draw(patch)
    draw.text((4, 4), text, fill=fill, font=font)
    rotated = patch.rotate(90, expand=True)
    img.paste(rotated, xy, rotated)


def _soft_line_color(color: str) -> str:
    mapping = {
        "blue": "#a7b8ff",
        "green": "#a7d8a7",
        "#ff9900": "#ffd799",
        "red": "#ff9a9a",
    }
    return mapping.get(color, color)
