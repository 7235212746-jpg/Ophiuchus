from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import tkinter as tk
from tkinter import ttk
from typing import Iterable

from .theme import COLORS, FONTS, SPACING


@dataclass(frozen=True)
class Element:
    atomic_number: int
    symbol: str
    period: int
    group: int
    display_row: int
    display_column: int


_PERIOD_ROWS: tuple[tuple[int, tuple[tuple[str, int], ...]], ...] = (
    (1, (("H", 1), ("He", 18))),
    (2, (("Li", 1), ("Be", 2), ("B", 13), ("C", 14), ("N", 15), ("O", 16), ("F", 17), ("Ne", 18))),
    (3, (("Na", 1), ("Mg", 2), ("Al", 13), ("Si", 14), ("P", 15), ("S", 16), ("Cl", 17), ("Ar", 18))),
    (4, tuple(zip(("K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"), range(1, 19)))),
    (5, tuple(zip(("Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"), range(1, 19)))),
    (6, (("Cs", 1), ("Ba", 2), ("Hf", 4), ("Ta", 5), ("W", 6), ("Re", 7), ("Os", 8), ("Ir", 9), ("Pt", 10), ("Au", 11), ("Hg", 12), ("Tl", 13), ("Pb", 14), ("Bi", 15), ("Po", 16), ("At", 17), ("Rn", 18))),
    (7, (("Fr", 1), ("Ra", 2), ("Rf", 4), ("Db", 5), ("Sg", 6), ("Bh", 7), ("Hs", 8), ("Mt", 9), ("Ds", 10), ("Rg", 11), ("Cn", 12), ("Nh", 13), ("Fl", 14), ("Mc", 15), ("Lv", 16), ("Ts", 17), ("Og", 18))),
)

_LANTHANIDES = ("La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu")
_ACTINIDES = ("Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr")


def _build_elements() -> tuple[Element, ...]:
    by_symbol: dict[str, tuple[int, int, int]] = {}
    for period, row in _PERIOD_ROWS:
        for symbol, group in row:
            by_symbol[symbol] = (period, group, period)
    for column, symbol in enumerate(_LANTHANIDES, start=3):
        by_symbol[symbol] = (6, 3, 8)
    for column, symbol in enumerate(_ACTINIDES, start=3):
        by_symbol[symbol] = (7, 3, 9)

    atomic_symbols = (
        "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
        "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg "
        "Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
    ).split()
    elements: list[Element] = []
    lanthanide_columns = {symbol: column for column, symbol in enumerate(_LANTHANIDES, start=3)}
    actinide_columns = {symbol: column for column, symbol in enumerate(_ACTINIDES, start=3)}
    for atomic_number, symbol in enumerate(atomic_symbols, start=1):
        period, group, display_row = by_symbol[symbol]
        display_column = lanthanide_columns.get(symbol, actinide_columns.get(symbol, group))
        elements.append(Element(atomic_number, symbol, period, group, display_row, display_column))
    return tuple(elements)


ELEMENTS = _build_elements()
ELEMENT_BY_SYMBOL = {element.symbol: element for element in ELEMENTS}
ATOMIC_ORDER = {element.symbol: element.atomic_number for element in ELEMENTS}


def parse_element_symbols(value: str | Iterable[str]) -> tuple[str, ...]:
    raw = re.split(r"[\s,，;；/]+", value.strip()) if isinstance(value, str) else [str(item).strip() for item in value]
    parsed: list[str] = []
    invalid: list[str] = []
    for token in raw:
        if not token:
            continue
        normalized = token[:1].upper() + token[1:].lower()
        if normalized not in ELEMENT_BY_SYMBOL:
            invalid.append(token)
        elif normalized not in parsed:
            parsed.append(normalized)
    if invalid:
        raise ValueError(f"invalid element symbol(s): {', '.join(invalid)}")
    return tuple(parsed)


def _elements_from_formula_token(token: str) -> tuple[str, ...]:
    if not re.fullmatch(r"(?:[A-Z][a-z]?\d*(?:\.\d+)?)+", token):
        return ()
    symbols: list[str] = []
    position = 0
    matcher = re.compile(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)")
    while position < len(token):
        match = matcher.match(token, position)
        if match is None or match.group(1) not in ELEMENT_BY_SYMBOL:
            return ()
        symbol = match.group(1)
        if symbol not in symbols:
            symbols.append(symbol)
        position = match.end()
    return tuple(symbols) if len(symbols) >= 2 else ()


def infer_elements_from_xrd_path(path: str | Path) -> tuple[str, ...]:
    source = Path(path)
    for text in (source.stem, source.parent.name):
        for token in re.split(r"[\s_+\-()\[\]]+", text):
            inferred = _elements_from_formula_token(token)
            if inferred:
                return inferred
    return ()


def element_scope_mismatch(
    selected: str | Iterable[str], inferred: str | Iterable[str]
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    selected_symbols = parse_element_symbols(selected)
    inferred_symbols = parse_element_symbols(inferred)
    selected_set = set(selected_symbols)
    inferred_set = set(inferred_symbols)
    missing = tuple(symbol for symbol in inferred_symbols if symbol not in selected_set)
    extra = tuple(symbol for symbol in selected_symbols if symbol not in inferred_set)
    return missing, extra


class PeriodicTableDialog(tk.Toplevel):
    def __init__(self, parent: tk.Misc, selected: str | Iterable[str] = ()) -> None:
        super().__init__(parent)
        self.title("选择主元素")
        self.configure(bg=COLORS["background"])
        self.resizable(False, False)
        self.transient(parent)
        self.result: tuple[str, ...] | None = None
        try:
            initial = parse_element_symbols(selected)
        except ValueError:
            initial = ()
        self._selected = set(initial)
        self._buttons: dict[str, tk.Button] = {}
        self._summary_var = tk.StringVar()
        self._build()
        self._refresh()
        self.protocol("WM_DELETE_WINDOW", self._cancel)
        self.bind("<Escape>", lambda _event: self._cancel())
        self.bind("<Return>", lambda _event: self._confirm())
        self.grab_set()
        self.update_idletasks()
        x = parent.winfo_rootx() + max(20, (parent.winfo_width() - self.winfo_reqwidth()) // 2)
        y = parent.winfo_rooty() + max(20, (parent.winfo_height() - self.winfo_reqheight()) // 2)
        self.geometry(f"+{x}+{y}")

    def _build(self) -> None:
        shell = ttk.Frame(self, style="Shell.TFrame", padding=SPACING["panel_pad"])
        shell.pack(fill="both", expand=True)
        ttk.Label(shell, text="元素周期表", style="Title.TLabel").grid(row=0, column=0, columnspan=18, sticky="w")
        ttk.Label(shell, text="点击元素进行多选，确认后会按原子序写入主元素范围。", style="Muted.TLabel").grid(
            row=1, column=0, columnspan=18, sticky="w", pady=(2, 12)
        )
        for group in range(1, 19):
            ttk.Label(shell, text=str(group), style="Muted.TLabel", anchor="center").grid(row=2, column=group - 1, padx=2)
        for element in ELEMENTS:
            button = tk.Button(
                shell,
                text=f"{element.atomic_number}\n{element.symbol}",
                width=4,
                height=2,
                font=FONTS["small"],
                bd=1,
                relief="flat",
                command=lambda symbol=element.symbol: self._toggle(symbol),
            )
            button.grid(row=element.display_row + 2, column=element.display_column - 1, padx=2, pady=2, sticky="nsew")
            self._buttons[element.symbol] = button
        ttk.Label(shell, textvariable=self._summary_var, style="Muted.TLabel").grid(
            row=12, column=0, columnspan=18, sticky="w", pady=(12, 8)
        )
        actions = ttk.Frame(shell, style="Shell.TFrame")
        actions.grid(row=13, column=0, columnspan=18, sticky="e")
        ttk.Button(actions, text="清空", command=self._clear).pack(side="left")
        ttk.Button(actions, text="取消", command=self._cancel).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="确认", style="Primary.TButton", command=self._confirm).pack(side="left", padx=(8, 0))

    def _toggle(self, symbol: str) -> None:
        if symbol in self._selected:
            self._selected.remove(symbol)
        else:
            self._selected.add(symbol)
        self._refresh()

    def _clear(self) -> None:
        self._selected.clear()
        self._refresh()

    def _refresh(self) -> None:
        ordered = sorted(self._selected, key=ATOMIC_ORDER.__getitem__)
        self._summary_var.set("已选：" + (" ".join(ordered) if ordered else "无"))
        for symbol, button in self._buttons.items():
            selected = symbol in self._selected
            button.configure(
                bg=COLORS["accent"] if selected else COLORS["panel"],
                fg="#ffffff" if selected else COLORS["text"],
                activebackground=COLORS["accent_hover"] if selected else COLORS["accent_soft"],
                highlightbackground=COLORS["border"],
            )

    def _confirm(self) -> None:
        self.result = tuple(sorted(self._selected, key=ATOMIC_ORDER.__getitem__))
        self.destroy()

    def _cancel(self) -> None:
        self.result = None
        self.destroy()


def show_periodic_table(parent: tk.Misc, selected: str | Iterable[str] = ()) -> tuple[str, ...] | None:
    dialog = PeriodicTableDialog(parent, selected)
    parent.wait_window(dialog)
    return dialog.result
