from __future__ import annotations


FONT_STACK = ("Segoe UI", "Microsoft YaHei UI", "Microsoft YaHei", "SF Pro Text", "PingFang SC", "Noto Sans CJK SC", "Noto Sans")

COLORS = {
    "background": "#eef3f8",
    "panel": "#ffffff",
    "panel_alt": "#f8fbff",
    "glass": "#ffffff",
    "border": "#d8e2ee",
    "border_soft": "#e6edf5",
    "text": "#172033",
    "muted": "#66758a",
    "subtle": "#8a98aa",
    "accent": "#3974d8",
    "accent_hover": "#2f65bf",
    "accent_soft": "#e8f0ff",
    "danger": "#b42318",
    "input": "#ffffff",
    "shadow": "#d6e0ec",
}

FONTS = {
    "body": (FONT_STACK[0], 10),
    "body_zh": (FONT_STACK[1], 10),
    "small": (FONT_STACK[0], 9),
    "title": (FONT_STACK[0], 24, "bold"),
    "section": (FONT_STACK[0], 12, "bold"),
    "button": (FONT_STACK[0], 10, "bold"),
}

SPACING = {
    "window_pad": 24,
    "panel_pad": 18,
    "card_pad": 12,
    "gap": 12,
    "small_gap": 6,
}


def text_font(size: int = 10, weight: str = "normal") -> tuple[str, int, str]:
    return (FONT_STACK[0], size, weight)


def configure_matplotlib_fonts() -> None:
    """Apply the same Chinese-capable fallback stack to Matplotlib figures."""
    from matplotlib import rcParams

    rcParams["font.family"] = "sans-serif"
    # Matplotlib does not reliably fall back per glyph on Windows, so start
    # with the installed YaHei face that contains both Latin and CJK glyphs.
    rcParams["font.sans-serif"] = [FONT_STACK[2], FONT_STACK[1], FONT_STACK[0], *FONT_STACK[3:]]
    rcParams["axes.unicode_minus"] = False
