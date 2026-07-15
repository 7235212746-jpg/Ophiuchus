from __future__ import annotations

import tkinter as tk
from tkinter import ttk
from typing import Callable

from ophiuchus.refinement.conclusion import ImpuritySignalEstimate, SampleConclusion
from ophiuchus.theme import COLORS, SPACING


class CandidateEvidenceDialog(tk.Toplevel):
    def __init__(
        self,
        parent: tk.Misc,
        conclusion: SampleConclusion,
        *,
        on_exclude_formula: Callable[[str], None] | None = None,
    ) -> None:
        super().__init__(parent)
        self.title("Ophi 候选证据")
        self.configure(bg=COLORS["background"])
        self.transient(parent)
        self.minsize(760, 500)
        self.geometry(_centered_geometry(self, 920, 620))
        self._on_exclude_formula = on_exclude_formula
        self._estimates = {item.candidate_id: item for item in conclusion.impurity_estimates}

        shell = ttk.Frame(self, style="Shell.TFrame", padding=SPACING["window_pad"])
        shell.pack(fill="both", expand=True)
        shell.columnconfigure(0, weight=1)
        shell.rowconfigure(2, weight=1)
        ttk.Label(shell, text="候选证据", style="Title.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(
            shell,
            text="只展示目标相残差上的增量证据；缺失理论强峰和边界参数会降低可信度。",
            style="Muted.TLabel",
        ).grid(row=1, column=0, sticky="w", pady=(3, 10))

        body = ttk.Panedwindow(shell, orient="horizontal")
        body.grid(row=2, column=0, sticky="nsew")
        table_panel = ttk.Frame(body, style="Panel.TFrame", padding=8)
        detail_panel = ttk.Frame(body, style="Panel.TFrame", padding=8)
        body.add(table_panel, weight=5)
        body.add(detail_panel, weight=5)

        self.tree = ttk.Treeview(
            table_panel,
            columns=("formula", "status", "gain"),
            show="headings",
            height=12,
        )
        self.tree.heading("formula", text="候选")
        self.tree.heading("status", text="状态")
        self.tree.heading("gain", text="残差改善")
        self.tree.column("formula", width=60, minwidth=55, anchor="w", stretch=True)
        self.tree.column("status", width=70, minwidth=65, anchor="w", stretch=True)
        self.tree.column("gain", width=45, minwidth=42, anchor="e", stretch=False)
        self.tree.pack(fill="both", expand=True)
        for estimate in conclusion.impurity_estimates:
            self.tree.insert(
                "",
                "end",
                iid=estimate.candidate_id,
                values=(estimate.formula, _short_status(estimate), f"{estimate.rwp_improvement:.2f}"),
            )
        self.tree.bind("<<TreeviewSelect>>", self._show_selected)

        self.detail_text = tk.Text(
            detail_panel,
            wrap="word",
            relief="flat",
            borderwidth=0,
            padx=12,
            pady=10,
            bg=COLORS["panel"],
            fg=COLORS["text"],
            font=("Microsoft YaHei UI", 10),
            spacing3=5,
        )
        self.detail_text.pack(fill="both", expand=True)

        actions = ttk.Frame(shell, style="Shell.TFrame")
        actions.grid(row=3, column=0, sticky="ew", pady=(10, 0))
        self.exclude_button = ttk.Button(actions, text="排除此候选并重算", command=self._exclude_selected)
        self.exclude_button.pack(side="left")
        if on_exclude_formula is None or not conclusion.impurity_estimates:
            self.exclude_button.state(["disabled"])
        ttk.Button(actions, text="关闭", command=self.destroy).pack(side="right")

        children = self.tree.get_children()
        if children:
            self.tree.selection_set(children[0])
            self.tree.focus(children[0])
            self._show_selected()
        self.lift()

    def _selected_estimate(self) -> ImpuritySignalEstimate | None:
        selected = self.tree.selection()
        return self._estimates.get(selected[0]) if selected else None

    def _show_selected(self, _event=None) -> None:
        estimate = self._selected_estimate()
        text = estimate.to_chinese_evidence_text() if estimate is not None else "请选择一个候选相。"
        self.detail_text.configure(state="normal")
        self.detail_text.delete("1.0", "end")
        self.detail_text.insert("1.0", text)
        self.detail_text.configure(state="disabled")

    def _exclude_selected(self) -> None:
        estimate = self._selected_estimate()
        if estimate is None or self._on_exclude_formula is None:
            return
        self.exclude_button.state(["disabled"])
        self._on_exclude_formula(estimate.formula)


def _centered_geometry(window: tk.Misc, width: int, height: int) -> str:
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()
    actual_width = min(width, screen_width - 80)
    actual_height = min(height, screen_height - 80)
    left = max(0, (screen_width - actual_width) // 2)
    top = max(0, (screen_height - actual_height) // 2)
    return f"{actual_width}x{actual_height}+{left}+{top}"


def _short_status(estimate: ImpuritySignalEstimate) -> str:
    if estimate.peak_width_at_boundary:
        return "不稳定"
    if estimate.rwp_improvement >= 2.0:
        return "较强候选"
    if estimate.rwp_improvement >= 0.75:
        return "可能存在"
    return "弱候选"
