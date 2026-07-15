from __future__ import annotations

import json
from pathlib import Path


def render_xrd_from_analysis_json(analysis_json: str | Path, out_png: str | Path, out_svg: str | Path | None = None) -> dict[str, str]:
    payload = json.loads(Path(analysis_json).read_text(encoding="utf-8"))
    outputs = render_xrd_payload(payload, out_png)
    if out_svg:
        outputs["svg"] = _write_simple_svg(payload, out_svg)
    return outputs


def render_xrd_payload(payload: dict[str, object], out_png: str | Path) -> dict[str, str]:
    from PIL import Image, ImageDraw, ImageFont

    out = Path(out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    width, height = 2100, 1250
    margin_l, margin_r, margin_t, margin_b = 120, 70, 70, 95
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = _font(28)
    small = _font(22)
    title_font = _font(34)

    peaks = payload.get("experimental_peaks") or []
    top = payload.get("top_candidates") or []
    input_summary = payload.get("input") or {}
    if not isinstance(peaks, list) or not peaks:
        img.save(out, dpi=(300, 300))
        return {"png": str(out)}
    xs = [float(p["two_theta"]) for p in peaks if isinstance(p, dict)]
    ys = [float(p["intensity"]) for p in peaks if isinstance(p, dict)]
    x_range = input_summary.get("two_theta_range") if isinstance(input_summary, dict) else None
    if isinstance(x_range, list) and len(x_range) == 2:
        x_min, x_max = float(x_range[0]), float(x_range[1])
    else:
        x_min, x_max = min(xs), max(xs)
    y_max = max(ys) or 1.0
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b
    y_top = 140.0
    y_bottom = -10.0 - max(0, min(len(top), 8) - 4) * 9.0

    def px(x: float) -> int:
        return int(margin_l + (x - x_min) / max(x_max - x_min, 1e-9) * plot_w)

    def py(y: float) -> int:
        return int(margin_t + (y_top - y) / max(y_top - y_bottom, 1e-9) * plot_h)

    draw.text((margin_l, 22), "Ophiuchus XRD candidate evidence", fill="#1f2937", font=title_font)
    draw.rectangle((margin_l, margin_t, width - margin_r, height - margin_b), outline="#d1d5db", width=2)
    _axis_ticks(draw, x_min, x_max, px, py, height - margin_b, small)

    exp_points = sorted([(float(p["two_theta"]), float(p["intensity"]) / y_max * 72.0 + 58.0) for p in peaks if isinstance(p, dict)])
    if len(exp_points) > 1:
        draw.line([(px(x), py(y)) for x, y in exp_points], fill="#111827", width=3)
    for peak in sorted(peaks, key=lambda p: float(p.get("intensity", 0)), reverse=True)[:12]:
        x = float(peak["two_theta"])
        y = float(peak["intensity"]) / y_max * 72.0 + 58.0
        draw.ellipse((px(x) - 4, py(y) - 4, px(x) + 4, py(y) + 4), fill="#111827")

    colors = ["#dc2626", "#2563eb", "#059669", "#d97706", "#7c3aed", "#0891b2", "#be123c", "#4b5563"]
    legend = [("Experimental peaks", "#111827")]
    for idx, item in enumerate(top[:8] if isinstance(top, list) else []):
        if not isinstance(item, dict):
            continue
        color = colors[idx % len(colors)]
        offset = 38.0 - idx * 9.5
        label = f"{idx + 1}. {item.get('formula', '')} ({float(item.get('score', 0)):.2f})"
        draw.line((px(x_min), py(offset), px(x_max), py(offset)), fill=color, width=2)
        simulated = item.get("simulated_pattern")
        if isinstance(simulated, dict):
            theory_x = simulated.get("two_theta_deg") or []
            theory_y = simulated.get("normalized_intensity") or []
            for x_value, intensity_value in zip(theory_x, theory_y):
                x = float(x_value)
                if x_min <= x <= x_max:
                    intensity = float(intensity_value)
                    draw.line((px(x), py(offset), px(x), py(offset + 7.0 * intensity / 100.0)), fill=color, width=2)
        else:
            for match in item.get("matched_peaks", []) or []:
                if not isinstance(match, dict):
                    continue
                x = float(match["theory_two_theta"])
                inten = float(match.get("theory_intensity", 50.0))
                draw.line((px(x), py(offset), px(x), py(offset + 7.0 * inten / 100.0)), fill=color, width=3)
        legend.append((label, color))

    legend_x, legend_y = width - margin_r - 520, margin_t + 24
    draw.rectangle((legend_x - 16, legend_y - 14, width - margin_r - 12, legend_y + 34 * min(len(legend), 9)), fill="#ffffff", outline="#e5e7eb")
    for i, (label, color) in enumerate(legend[:9]):
        y = legend_y + i * 34
        draw.line((legend_x, y + 13, legend_x + 54, y + 13), fill=color, width=4)
        draw.text((legend_x + 70, y), label[:44], fill="#111827", font=small)

    draw.text((width // 2 - 100, height - 62), "2θ (degree)", fill="#111827", font=font)
    draw.text((24, margin_t + 20), "Intensity / offset a.u.", fill="#111827", font=small)
    img.save(out, dpi=(300, 300))
    return {"png": str(out)}


def _write_simple_svg(payload: dict[str, object], out_svg: str | Path) -> str:
    out = Path(out_svg)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        "<svg xmlns='http://www.w3.org/2000/svg' width='1200' height='700'>"
        "<rect width='100%' height='100%' fill='white'/>"
        "<text x='40' y='50' font-family='Arial' font-size='24'>Ophiuchus XRD plot SVG placeholder</text>"
        "<text x='40' y='85' font-family='Arial' font-size='14'>Use PNG for full rendered evidence in this build.</text>"
        "</svg>\n",
        encoding="utf-8",
    )
    return str(out)


def _axis_ticks(draw, x_min, x_max, px, py, y_axis, font) -> None:
    import math

    span = max(x_max - x_min, 1.0)
    step = 10.0 if span > 45 else 5.0 if span > 20 else 2.0
    start = math.ceil(x_min / step) * step
    value = start
    while value <= x_max + 1e-9:
        x = px(value)
        draw.line((x, y_axis, x, y_axis + 10), fill="#6b7280", width=2)
        draw.text((x - 20, y_axis + 16), f"{value:g}", fill="#374151", font=font)
        value += step


def _dashed_line(draw, x: int, y0: int, y1: int, color: str) -> None:
    y = min(y0, y1)
    end = max(y0, y1)
    while y < end:
        draw.line((x, y, x, min(y + 8, end)), fill=color, width=2)
        y += 13


def _font(size: int):
    from PIL import ImageFont

    for candidate in [r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\segoeui.ttf", r"C:\Windows\Fonts\msyh.ttc"]:
        try:
            return ImageFont.truetype(candidate, size=size)
        except Exception:
            pass
    return ImageFont.load_default()
