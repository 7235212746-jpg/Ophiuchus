from __future__ import annotations

import math
import os
from pathlib import Path
import subprocess
import time
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Callable, Iterable

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from ophiuchus.theme import COLORS, FONTS, SPACING, configure_matplotlib_fonts
from ophiuchus.xrd.models import Candidate, Peak

from .background import estimate_xrd_background
from .composition import inspect_peak_composition
from .export import export_phase_stripping_session
from .models import AnalysisContext
from .ranking import CandidateEvidence, extract_residual_peaks
from .session import PhasePreview, PhaseStrippingSession


class PhaseStrippingWindow(tk.Toplevel):
    """Interactive, non-destructive manual phase stripping workbench."""

    def __init__(
        self,
        parent: tk.Misc,
        context: AnalysisContext,
        candidates: Iterable[Candidate],
        *,
        element_scope: Iterable[str] = (),
        on_closed: Callable[[], None] | None = None,
    ) -> None:
        super().__init__(parent)
        self.title("手动物相剥离 / 残差谱分析")
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        window_width = min(1520, max(1180, screen_width - 80))
        window_height = min(860, max(700, screen_height - 100))
        self.geometry(f"{window_width}x{window_height}")
        self.minsize(
            max(900, min(1120, screen_width - 120)),
            max(620, min(680, screen_height - 140)),
        )
        self.configure(bg=COLORS["background"])
        background = estimate_xrd_background(context.x, context.intensity)
        self.session = PhaseStrippingSession(
            context,
            background_y=background.values,
            background_method=background.method,
            background_parameters=background.parameters,
        )
        self.candidates = list(candidates)
        self.element_scope = tuple(element_scope)
        self._on_closed = on_closed
        self._closing = False
        self._rankings: list[CandidateEvidence] = []
        self._evidence_by_candidate_id: dict[str, CandidateEvidence] = {}
        self._row_candidates: dict[str, Candidate] = {}
        self._selected_candidate: Candidate | None = None
        self._preview: PhasePreview | None = None
        self._last_export_path: Path | None = None
        self._default_sashes_positioned = False

        self.search_var = tk.StringVar()
        self.sort_var = tk.StringVar(value="综合得分")
        self.show_sticks_var = tk.BooleanVar(value=True)
        self.scale_var = tk.StringVar(value="0")
        self.shift_var = tk.StringVar(value="0")
        self.sigma_var = tk.StringVar(value="0.060")
        self.range_min_var = tk.StringVar(value=f"{context.two_theta_range[0]:.3f}")
        self.range_max_var = tk.StringVar(value=f"{context.two_theta_range[1]:.3f}")
        self.readout_var = tk.StringVar(value="在图中点击可读取 2theta 与强度")
        self.warning_var = tk.StringVar(value="拟合比例仅用于候选证据，不是 Rietveld 定量相含量。")
        self.status_var = tk.StringVar(
            value=f"继承：{context.radiation} / {context.wavelength_angstrom:.6f} A / "
            f"{context.two_theta_range[0]:.2f}-{context.two_theta_range[1]:.2f} deg / 背景：AsLS（固定）"
        )

        self._build_layout()
        self._bind_local_shortcuts()
        self.protocol("WM_DELETE_WINDOW", self._request_close)
        self._refresh_all()

    def _build_layout(self) -> None:
        shell = ttk.Frame(self, style="Shell.TFrame", padding=SPACING["window_pad"])
        shell.pack(fill="both", expand=True)
        header = ttk.Frame(shell, style="Shell.TFrame")
        header.pack(fill="x", pady=(0, 12))
        ttk.Label(header, text="手动物相剥离 / 残差谱分析", style="Title.TLabel").pack(side="left")
        ttk.Label(header, textvariable=self.status_var, style="Muted.TLabel").pack(side="right", anchor="e")

        self.panes = ttk.Panedwindow(shell, orient="horizontal")
        self.panes.pack(fill="both", expand=True)
        self.left_pane = ttk.Frame(self.panes, style="Panel.TFrame", padding=12)
        self.center_pane = ttk.Frame(self.panes, style="Panel.TFrame", padding=12)
        self.right_pane = ttk.Frame(self.panes, style="Panel.TFrame", padding=12)
        self.panes.add(self.left_pane, weight=2)
        self.panes.add(self.center_pane, weight=5)
        self.panes.add(self.right_pane, weight=3)
        self._build_candidate_pane()
        self._build_chart_pane()
        self._build_control_pane()
        self.panes.bind("<Configure>", self._position_default_sashes, add="+")

    def _position_default_sashes(self, _event: tk.Event | None = None) -> None:
        if self._default_sashes_positioned:
            return
        total_width = self.panes.winfo_width()
        if total_width < 1000:
            return
        left_width = min(340, max(300, total_width // 5))
        right_width = min(440, max(410, total_width // 4))
        if total_width - left_width - right_width < 520:
            return
        self._default_sashes_positioned = True
        self.panes.sashpos(0, left_width)
        self.panes.sashpos(1, total_width - right_width)

    def _build_candidate_pane(self) -> None:
        pane = self.left_pane
        pane.rowconfigure(3, weight=1)
        pane.columnconfigure(0, weight=1)
        ttk.Label(pane, text="候选相", style="Section.TLabel").grid(row=0, column=0, sticky="w")
        search = ttk.Entry(pane, textvariable=self.search_var)
        search.grid(row=1, column=0, sticky="ew", pady=(8, 6))
        search.bind("<KeyRelease>", lambda _event: self._populate_candidate_tree())
        sort = ttk.Combobox(
            pane,
            textvariable=self.sort_var,
            values=("综合得分", "独立峰证据", "化学式", "来源"),
            state="readonly",
        )
        sort.grid(row=2, column=0, sticky="ew", pady=(0, 6))
        sort.bind("<<ComboboxSelected>>", lambda _event: self._populate_candidate_tree())

        columns = ("formula", "sg", "source", "state", "matched", "missing", "independent", "score")
        self.candidate_tree = ttk.Treeview(pane, columns=columns, show="headings", height=18)
        headings = {
            "formula": "化学式",
            "sg": "空间群",
            "source": "来源",
            "state": "状态",
            "matched": "命中",
            "missing": "缺失",
            "independent": "独立",
            "score": "得分",
        }
        widths = {"formula": 90, "sg": 75, "source": 90, "state": 62, "matched": 48, "missing": 48, "independent": 48, "score": 58}
        for column in columns:
            self.candidate_tree.heading(column, text=headings[column])
            self.candidate_tree.column(column, width=widths[column], minwidth=42, stretch=column in {"formula", "source"})
        ybar = ttk.Scrollbar(pane, orient="vertical", command=self.candidate_tree.yview)
        xbar = ttk.Scrollbar(pane, orient="horizontal", command=self.candidate_tree.xview)
        self.candidate_tree.configure(yscrollcommand=ybar.set, xscrollcommand=xbar.set)
        self.candidate_tree.grid(row=3, column=0, sticky="nsew")
        ybar.grid(row=3, column=1, sticky="ns")
        xbar.grid(row=4, column=0, sticky="ew")
        self.candidate_tree.bind("<<TreeviewSelect>>", lambda _event: self._on_candidate_selected())

        actions = ttk.Frame(pane, style="Panel.TFrame")
        actions.grid(row=5, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        for index in range(2):
            actions.columnconfigure(index, weight=1)
        self.preview_button = ttk.Button(actions, text="预览", command=self._auto_fit_preview)
        self.preview_button.grid(row=0, column=0, sticky="ew", padx=(0, 4))
        self.accept_button = ttk.Button(actions, text="接受", style="Primary.TButton", command=self._accept_preview)
        self.accept_button.grid(row=0, column=1, sticky="ew", padx=(4, 0))
        self.exclude_button = ttk.Button(actions, text="排除", command=self._exclude_selected)
        self.exclude_button.grid(row=1, column=0, sticky="ew", padx=(0, 4), pady=(6, 0))
        self.view_cif_button = ttk.Button(actions, text="查看 CIF", command=self._view_cif)
        self.view_cif_button.grid(row=1, column=1, sticky="ew", padx=(4, 0), pady=(6, 0))

    def _build_chart_pane(self) -> None:
        pane = self.center_pane
        pane.rowconfigure(1, weight=1)
        pane.columnconfigure(0, weight=1)
        bar = ttk.Frame(pane, style="Panel.TFrame")
        bar.grid(row=0, column=0, sticky="ew")
        ttk.Label(bar, text="峰重构：原谱、物相贡献与残差", style="Section.TLabel").pack(side="left")
        ttk.Checkbutton(bar, text="峰棒", variable=self.show_sticks_var, command=self._refresh_chart).pack(side="right")
        configure_matplotlib_fonts()
        self.figure = Figure(figsize=(8.0, 6.0), dpi=100, facecolor=COLORS["panel"])
        axes = self.figure.subplots(3, 1, sharex=True, gridspec_kw={"height_ratios": (3.2, 2.0, 1.5)})
        self.pattern_axis, self.contribution_axis, self.residual_axis = axes
        self.axis = self.pattern_axis
        self.canvas = FigureCanvasTkAgg(self.figure, master=pane)
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        toolbar_frame = ttk.Frame(pane, style="Panel.TFrame")
        toolbar_frame.grid(row=2, column=0, sticky="ew")
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side="left", fill="x")
        self.readout_label = ttk.Label(pane, textvariable=self.readout_var, style="PanelMuted.TLabel", anchor="w")
        self.readout_label.grid(row=3, column=0, sticky="ew", pady=(4, 0))
        self.canvas.mpl_connect("button_press_event", self._on_chart_click)

    def _build_control_pane(self) -> None:
        pane = self.right_pane
        pane.rowconfigure(2, weight=1)
        pane.columnconfigure(0, weight=1)
        ttk.Label(pane, text="峰组成与拟合控制", style="Section.TLabel").grid(row=0, column=0, sticky="w")
        summary = ttk.Label(pane, textvariable=self.warning_var, style="PanelMuted.TLabel", wraplength=390, justify="left")
        summary.grid(row=1, column=0, sticky="ew", pady=(6, 8))

        self.detail_tabs = ttk.Notebook(pane)
        self.detail_tabs.grid(row=2, column=0, sticky="nsew")
        composition_tab = ttk.Frame(self.detail_tabs, style="Panel.TFrame", padding=4)
        peaks_tab = ttk.Frame(self.detail_tabs, style="Panel.TFrame", padding=4)
        self.detail_tabs.add(composition_tab, text="峰组成")
        self.detail_tabs.add(peaks_tab, text="候选峰")
        composition_tab.rowconfigure(0, weight=1)
        composition_tab.columnconfigure(0, weight=1)
        peaks_tab.rowconfigure(0, weight=1)
        peaks_tab.columnconfigure(0, weight=1)

        composition_columns = ("label", "intensity", "share", "hkl", "reflection", "delta")
        self.composition_tree = ttk.Treeview(composition_tab, columns=composition_columns, show="headings", height=11)
        composition_headings = {
            "label": "组成",
            "intensity": "该点强度",
            "share": "已解释占比",
            "hkl": "hkl",
            "reflection": "最近峰位",
            "delta": "差值",
        }
        composition_widths = {"label": 80, "intensity": 60, "share": 60, "hkl": 42, "reflection": 60, "delta": 48}
        for column in composition_columns:
            self.composition_tree.heading(column, text=composition_headings[column])
            self.composition_tree.column(column, width=composition_widths[column], minwidth=44, anchor="center")
        composition_scroll = ttk.Scrollbar(composition_tab, orient="vertical", command=self.composition_tree.yview)
        self.composition_tree.configure(yscrollcommand=composition_scroll.set)
        self.composition_tree.grid(row=0, column=0, sticky="nsew")
        composition_scroll.grid(row=0, column=1, sticky="ns")

        columns = ("two_theta", "d", "hkl", "intensity", "hit", "overlap")
        self.peak_tree = ttk.Treeview(peaks_tab, columns=columns, show="headings", height=11)
        headings = {"two_theta": "2theta", "d": "d(A)", "hkl": "hkl", "intensity": "强度", "hit": "残差命中", "overlap": "重叠"}
        for column in columns:
            self.peak_tree.heading(column, text=headings[column])
            self.peak_tree.column(column, width=62, minwidth=45, anchor="center")
        peak_scroll = ttk.Scrollbar(peaks_tab, orient="vertical", command=self.peak_tree.yview)
        self.peak_tree.configure(yscrollcommand=peak_scroll.set)
        self.peak_tree.grid(row=0, column=0, sticky="nsew")
        peak_scroll.grid(row=0, column=1, sticky="ns")

        controls = ttk.Frame(pane, style="Panel.TFrame")
        controls.grid(row=3, column=0, sticky="ew", pady=(10, 0))
        controls.columnconfigure(1, weight=1)
        fields = (
            ("比例 >= 0", self.scale_var),
            ("整体偏移 (deg)", self.shift_var),
            ("峰宽 sigma (deg)", self.sigma_var),
            ("拟合范围起点", self.range_min_var),
            ("拟合范围终点", self.range_max_var),
        )
        for row, (label, variable) in enumerate(fields):
            ttk.Label(controls, text=label, style="Panel.TLabel").grid(row=row, column=0, sticky="w", pady=2)
            ttk.Entry(controls, textvariable=variable, width=12).grid(row=row, column=1, sticky="ew", padx=(8, 0), pady=2)
        fit_actions = ttk.Frame(controls, style="Panel.TFrame")
        fit_actions.grid(row=len(fields), column=0, columnspan=2, sticky="ew", pady=(8, 0))
        fit_actions.columnconfigure(0, weight=1)
        fit_actions.columnconfigure(1, weight=1)
        ttk.Button(fit_actions, text="自动拟合", command=self._auto_fit_preview).grid(row=0, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(fit_actions, text="按参数预览", command=self._manual_preview).grid(row=0, column=1, sticky="ew", padx=(4, 0))

        history = ttk.Frame(pane, style="Panel.TFrame")
        history.grid(row=4, column=0, sticky="ew", pady=(10, 0))
        for column in range(3):
            history.columnconfigure(column, weight=1)
        self.cancel_button = ttk.Button(history, text="取消预览", command=self._cancel_preview)
        self.cancel_button.grid(row=0, column=0, sticky="ew", padx=(0, 3))
        self.undo_button = ttk.Button(history, text="撤销", command=self._undo)
        self.undo_button.grid(row=0, column=1, sticky="ew", padx=3)
        self.redo_button = ttk.Button(history, text="重做", command=self._redo)
        self.redo_button.grid(row=0, column=2, sticky="ew", padx=(3, 0))
        ttk.Button(history, text="重置", command=self._reset).grid(row=1, column=0, columnspan=2, sticky="ew", padx=(0, 3), pady=(6, 0))
        ttk.Button(history, text="导出", command=self._export).grid(row=1, column=2, sticky="ew", padx=(3, 0), pady=(6, 0))

    def _bind_local_shortcuts(self) -> None:
        self.bind("<Control-z>", lambda _event: self._undo())
        self.bind("<Control-y>", lambda _event: self._redo())
        self.bind("<Escape>", lambda _event: self._cancel_preview())

    def _refresh_all(self) -> None:
        self._rankings = self.session.rank_candidates(self.candidates, element_scope=self.element_scope)
        self._evidence_by_candidate_id = {item.candidate.candidate_id: item for item in self._rankings}
        self._populate_candidate_tree()
        self._populate_peak_tree()
        self._refresh_chart()
        self._refresh_controls()

    def _populate_candidate_tree(self) -> None:
        selected_id = self._selected_candidate.candidate_id if self._selected_candidate else None
        query = self.search_var.get().strip().lower()
        rankings = [
            evidence
            for evidence in self._rankings
            if not query
            or query in evidence.candidate.formula_pretty.lower()
            or query in evidence.candidate.source.lower()
            or query in evidence.candidate.space_group_symbol.lower()
        ]
        sort_key = self.sort_var.get()
        if sort_key == "独立峰证据":
            rankings.sort(key=lambda item: (-item.independent_peak_evidence, -item.final_score))
        elif sort_key == "化学式":
            rankings.sort(key=lambda item: (item.candidate.formula_pretty.lower(), -item.final_score))
        elif sort_key == "来源":
            rankings.sort(key=lambda item: (item.candidate.source.lower(), -item.final_score))
        else:
            rankings.sort(key=lambda item: -item.final_score)
        self.candidate_tree.delete(*self.candidate_tree.get_children())
        self._row_candidates.clear()
        accepted_ids = {operation.phase_fit.candidate_id for operation in self.session.accepted_operations}
        selection_row = None
        for index, evidence in enumerate(rankings):
            candidate = evidence.candidate
            matched, missing, independent = self._candidate_peak_counts(candidate)
            state = "已接受" if candidate.candidate_id in accepted_ids else candidate.simulation_state or "可预览"
            row_id = f"candidate-{index}"
            self._row_candidates[row_id] = candidate
            self.candidate_tree.insert(
                "",
                "end",
                iid=row_id,
                values=(
                    candidate.formula_pretty,
                    candidate.space_group_symbol or "-",
                    candidate.source,
                    state,
                    matched,
                    missing,
                    independent,
                    f"{evidence.final_score:.3f}",
                ),
            )
            if candidate.candidate_id == selected_id:
                selection_row = row_id
        if selection_row:
            self.candidate_tree.selection_set(selection_row)
            self.candidate_tree.see(selection_row)

    def _candidate_peak_counts(self, candidate: Candidate) -> tuple[int, int, int]:
        residual_peaks = extract_residual_peaks(self.session.context.x, self.session.residual_y)
        tolerance = self.session.context.tolerance_deg
        strong = [peak for peak in candidate.theory_peaks if peak.intensity >= 10.0] or list(candidate.theory_peaks)
        matched = [
            peak for peak in strong if any(abs(peak.two_theta - residual.two_theta) <= tolerance for residual in residual_peaks)
        ]
        independent = [
            peak
            for peak in matched
            if not any(abs(peak.two_theta - accepted) <= tolerance for accepted in self.session.accepted_peak_positions)
        ]
        return len(matched), len(strong) - len(matched), len(independent)

    def _on_candidate_selected(self) -> None:
        selected = self.candidate_tree.selection()
        self._selected_candidate = self._row_candidates.get(selected[0]) if selected else None
        self.session.cancel_preview()
        self._preview = None
        if self._selected_candidate is not None:
            evidence = self._evidence_by_candidate_id.get(self._selected_candidate.candidate_id)
            if evidence is not None:
                detail = [*evidence.explanations[:2], *evidence.warnings]
                self.warning_var.set("\n".join(detail) if detail else "候选证据已加载。")
        self._populate_peak_tree()
        self._refresh_chart()
        self._refresh_controls()

    def _fit_mask(self) -> np.ndarray:
        start = float(self.range_min_var.get())
        end = float(self.range_max_var.get())
        lower, upper = self.session.context.two_theta_range
        if not (lower <= start < end <= upper):
            raise ValueError(f"拟合范围必须位于 {lower:.3f}-{upper:.3f} deg，且起点小于终点。")
        return (self.session.context.x >= start) & (self.session.context.x <= end)

    def _auto_fit_preview(self) -> None:
        candidate = self._selected_candidate
        if candidate is None:
            self._show_warning("请先在左侧选择候选相。")
            return
        if self._candidate_already_accepted(candidate):
            self._show_warning("该候选已经接受。请先撤销对应操作，避免重复扣除同一相。")
            return
        try:
            preview = self.session.preview(candidate, fit_mask=self._fit_mask())
        except Exception as exc:
            self._show_warning(f"自动拟合失败：{exc}")
            return
        self._preview = preview
        self.scale_var.set(f"{preview.phase_fit.scale:.8g}")
        self.shift_var.set(f"{preview.phase_fit.shift_deg:.6f}")
        self.sigma_var.set(f"{preview.phase_fit.sigma_deg:.6f}")
        self._after_preview_changed()

    def _manual_preview(self) -> None:
        candidate = self._selected_candidate
        if candidate is None:
            self._show_warning("请先在左侧选择候选相。")
            return
        if self._candidate_already_accepted(candidate):
            self._show_warning("该候选已经接受。请先撤销对应操作，避免重复扣除同一相。")
            return
        try:
            self._fit_mask()
            preview = self.session.preview_with_parameters(
                candidate,
                scale=float(self.scale_var.get()),
                shift_deg=float(self.shift_var.get()),
                sigma_deg=float(self.sigma_var.get()),
            )
        except Exception as exc:
            self._show_warning(f"参数预览失败：{exc}")
            return
        self._preview = preview
        self._after_preview_changed()

    def _after_preview_changed(self) -> None:
        self._populate_peak_tree()
        self._refresh_chart()
        self._refresh_controls()
        warnings = list(self._preview.warnings if self._preview else ())
        warnings.append("比例仅表示当前残差上的拟合尺度，不是定量相含量。")
        self.warning_var.set("\n".join(warnings))

    def _accept_preview(self) -> None:
        if self._preview is None:
            return
        try:
            self.session.accept_preview(self._preview)
        except Exception as exc:
            self._show_warning(f"接受候选失败：{exc}")
            return
        self._preview = None
        self.warning_var.set("已接受该贡献并从不可变原谱重算残差；候选已按新残差重新排序。")
        self._refresh_all()

    def _cancel_preview(self) -> None:
        self.session.cancel_preview()
        self._preview = None
        self.warning_var.set("预览已取消；原谱与已接受历史没有改变。")
        self._populate_peak_tree()
        self._refresh_chart()
        self._refresh_controls()

    def _exclude_selected(self) -> None:
        if self._selected_candidate is None:
            return
        self.session.exclude(self._selected_candidate)
        self._selected_candidate = None
        self._preview = None
        self._refresh_all()

    def _undo(self) -> None:
        if self.session.undo():
            self._preview = None
            self.warning_var.set("已撤销上一次接受操作。")
            self._refresh_all()

    def _redo(self) -> None:
        if self.session.redo():
            self._preview = None
            self.warning_var.set("已重做上一次接受操作。")
            self._refresh_all()

    def _reset(self) -> None:
        if self.session.accepted_operations and not messagebox.askyesno("重置剥离", "清空全部物相接受历史并恢复扣背景后的晶相信号？\n原始谱和固定背景模型不会改变。", parent=self):
            return
        self.session.reset()
        self._preview = None
        self.warning_var.set("已清空物相操作；原始谱与固定背景模型保持不变。")
        self._refresh_all()

    def _populate_peak_tree(self) -> None:
        self.peak_tree.delete(*self.peak_tree.get_children())
        candidate = self._selected_candidate
        if candidate is None:
            return
        residual_peaks = extract_residual_peaks(self.session.context.x, self.session.residual_y)
        tolerance = self.session.context.tolerance_deg
        shift = self._preview.phase_fit.shift_deg if self._preview else 0.0
        for index, peak in enumerate(sorted(candidate.theory_peaks, key=lambda item: item.two_theta)):
            position = peak.two_theta + shift
            hit = any(abs(position - residual.two_theta) <= tolerance for residual in residual_peaks)
            overlap = any(abs(position - accepted) <= tolerance for accepted in self.session.accepted_peak_positions)
            d_spacing = peak.d_spacing if peak.d_spacing is not None else self._d_spacing(position)
            self.peak_tree.insert(
                "",
                "end",
                iid=f"peak-{index}",
                values=(f"{position:.4f}", f"{d_spacing:.4f}" if d_spacing else "-", peak.hkl or "-", f"{peak.intensity:.3f}", "是" if hit else "否", "是" if overlap else "否"),
            )

    def _d_spacing(self, two_theta: float) -> float | None:
        sine = math.sin(math.radians(two_theta / 2.0))
        return None if sine <= 0.0 else self.session.context.wavelength_angstrom / (2.0 * sine)

    def _refresh_chart(self) -> None:
        old_x_limits = self.pattern_axis.get_xlim() if self.pattern_axis.lines else None
        for axis in (self.pattern_axis, self.contribution_axis, self.residual_axis):
            axis.clear()
        x = self.session.context.x
        experimental = self.session.context.intensity
        background = self.session.background_y
        explained = self.session.fitted_total
        reconstructed = self.session.reconstructed_y
        residual = self.session.residual_y
        self.pattern_axis.plot(x, experimental, color="#111827", linewidth=1.05, label="实验谱")
        self.pattern_axis.plot(x, background, color="#d48a17", linewidth=1.05, label="估计背景")
        self.pattern_axis.plot(x, reconstructed, color="#3974d8", linewidth=1.15, label="背景 + 已接受相")
        colors = ("#2a9d8f", "#8f5bd7", "#d48a17", "#0086b3", "#b84a62")
        candidates_by_id = {candidate.candidate_id: candidate for candidate in self.candidates}
        contribution_ticks: list[float] = []
        contribution_labels: list[str] = []
        for index, (operation, contribution) in enumerate(zip(self.session.accepted_operations, self.session.accepted_contributions)):
            candidate = candidates_by_id.get(operation.phase_fit.candidate_id)
            label = candidate.formula_pretty if candidate is not None else operation.phase_fit.candidate_id
            baseline = float(index)
            maximum = max(float(np.max(contribution)), np.finfo(float).eps)
            self.contribution_axis.plot(
                x,
                baseline + 0.72 * contribution / maximum,
                color=colors[index % len(colors)],
                linewidth=1.0,
            )
            self.contribution_axis.fill_between(
                x,
                baseline,
                baseline + 0.72 * contribution / maximum,
                color=colors[index % len(colors)],
                alpha=0.12,
            )
            contribution_ticks.append(baseline + 0.32)
            contribution_labels.append(label)
        if self._preview is not None:
            preview_total = background + explained + self._preview.contribution
            self.pattern_axis.plot(x, preview_total, color="#15a34a", linewidth=1.0, linestyle="--", label="含预览的重构")
            baseline = float(len(contribution_ticks))
            maximum = max(float(np.max(self._preview.contribution)), np.finfo(float).eps)
            self.contribution_axis.plot(
                x,
                baseline + 0.72 * self._preview.contribution / maximum,
                color="#15a34a",
                linewidth=1.0,
                linestyle="--",
            )
            candidate = candidates_by_id.get(self._preview.phase_fit.candidate_id)
            label = candidate.formula_pretty if candidate is not None else self._preview.phase_fit.candidate_id
            contribution_ticks.append(baseline + 0.32)
            contribution_labels.append(f"{label} (预览)")
        if self.show_sticks_var.get() and self._selected_candidate is not None:
            shift = self._preview.phase_fit.shift_deg if self._preview else 0.0
            positions = [peak.two_theta + shift for peak in self._selected_candidate.theory_peaks]
            ymax = max(float(np.max(experimental)), 1.0)
            heights = [0.08 * ymax * peak.intensity / 100.0 for peak in self._selected_candidate.theory_peaks]
            stick_bases = np.interp(positions, x, background)
            self.pattern_axis.vlines(positions, stick_bases, stick_bases + heights, color="#15a34a", linewidth=0.8, alpha=0.55, label="所选候选峰")

        self.residual_axis.plot(x, residual, color="#b42318", linewidth=0.9, label="有符号残差")
        self.residual_axis.fill_between(x, 0.0, residual, where=residual >= 0.0, color="#e05a47", alpha=0.23)
        self.residual_axis.fill_between(x, 0.0, residual, where=residual < 0.0, color="#3974d8", alpha=0.20)
        self.residual_axis.axhline(0.0, color="#7b8798", linewidth=0.75)

        self.pattern_axis.set_ylabel("强度 (a.u.)")
        self.contribution_axis.set_ylabel("各相峰形")
        self.residual_axis.set_ylabel("残差")
        self.residual_axis.set_xlabel("2theta (degree)")
        self.contribution_axis.set_yticks(contribution_ticks, contribution_labels, fontsize=8)
        if contribution_ticks:
            self.contribution_axis.set_ylim(-0.08, len(contribution_ticks) - 0.05)
        else:
            self.contribution_axis.set_ylim(-0.1, 0.9)
            self.contribution_axis.text(0.5, 0.5, "接受或预览候选后，这里会逐相显示峰形", transform=self.contribution_axis.transAxes, ha="center", va="center", color="#7b8798", fontsize=9)
        for axis in (self.pattern_axis, self.contribution_axis, self.residual_axis):
            axis.grid(color=COLORS["border_soft"], linewidth=0.55, alpha=0.62)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
        self.pattern_axis.legend(loc="upper right", fontsize=8, frameon=False, ncols=2)
        self.pattern_axis.set_xlim(self.session.context.two_theta_range)
        if old_x_limits and old_x_limits[1] > old_x_limits[0]:
            current_range = self.session.context.two_theta_range
            if old_x_limits[0] >= current_range[0] and old_x_limits[1] <= current_range[1]:
                self.pattern_axis.set_xlim(old_x_limits)
        self.figure.tight_layout()
        self.canvas.draw()
        self._populate_composition_tree()

    def _on_chart_click(self, event) -> None:
        if event.inaxes in (self.pattern_axis, self.contribution_axis, self.residual_axis) and event.xdata is not None:
            result = self._populate_composition_tree(float(event.xdata))
            self.readout_var.set(
                f"2theta {result.two_theta:.4f} | 实验 {self._format_point_intensity(result.experimental)} | "
                f"背景 {self._format_point_intensity(result.background)} | "
                f"扣背景 {self._format_point_intensity(result.corrected)} | "
                f"解释 {self._format_point_intensity(result.explained)} | "
                f"残差 {self._format_point_intensity(result.residual)}"
            )

    def _populate_composition_tree(self, two_theta: float | None = None):
        if two_theta is None:
            two_theta = float(self.session.context.x[int(np.argmax(self.session.context.intensity))])
        candidates_by_id = {candidate.candidate_id: candidate for candidate in self.candidates}
        result = inspect_peak_composition(
            self.session,
            candidates_by_id,
            two_theta,
            preview=self._preview,
        )
        self.composition_tree.delete(*self.composition_tree.get_children())
        for index, row in enumerate(result.rows):
            share = "-" if row.explained_share_percent is None else f"{row.explained_share_percent:.1f}%"
            reflection = "-" if row.reflection_two_theta is None else f"{row.reflection_two_theta:.4f}"
            delta = "-" if row.reflection_delta is None else f"{row.reflection_delta:+.4f}"
            self.composition_tree.insert(
                "",
                "end",
                iid=f"composition-{index}",
                values=(
                    row.label,
                    self._format_point_intensity(row.intensity),
                    share,
                    row.hkl or "-",
                    reflection,
                    delta,
                ),
            )
        return result

    @staticmethod
    def _format_point_intensity(value: float) -> str:
        if abs(float(value)) < 1e-8:
            return "0"
        return f"{float(value):.6g}"

    def _refresh_controls(self) -> None:
        selected = self._selected_candidate is not None
        accepted_ids = {operation.phase_fit.candidate_id for operation in self.session.accepted_operations}
        selected_is_accepted = bool(
            self._selected_candidate and self._selected_candidate.candidate_id in accepted_ids
        )
        preview = self._preview is not None
        self._set_enabled(self.preview_button, selected and not selected_is_accepted)
        self._set_enabled(self.exclude_button, selected)
        self._set_enabled(self.view_cif_button, selected and bool(self._selected_candidate and self._selected_candidate.source_path))
        self._set_enabled(self.accept_button, preview)
        self._set_enabled(self.cancel_button, preview)
        self._set_enabled(self.undo_button, bool(self.session.accepted_operations))
        self._set_enabled(self.redo_button, bool(self.session.to_dict()["redo"]))

    @staticmethod
    def _set_enabled(widget: ttk.Widget, enabled: bool) -> None:
        widget.state(["!disabled"] if enabled else ["disabled"])

    def _candidate_already_accepted(self, candidate: Candidate) -> bool:
        return any(
            operation.phase_fit.candidate_id == candidate.candidate_id
            for operation in self.session.accepted_operations
        )

    def _view_cif(self) -> None:
        if self._selected_candidate is None:
            return
        path = Path(self._selected_candidate.source_path)
        if not path.is_file():
            self._show_warning(f"CIF 文件不存在：{path}")
            return
        try:
            if os.name == "nt":
                os.startfile(path)  # type: ignore[attr-defined]
            else:
                subprocess.Popen(["xdg-open", str(path)])
        except Exception as exc:
            self._show_warning(f"无法打开 CIF：{exc}")

    def _export(self) -> bool:
        parent = filedialog.askdirectory(title="选择残差谱导出的上级文件夹", parent=self)
        if not parent:
            return False
        base = Path(parent) / time.strftime("phase_stripping_%Y%m%d_%H%M%S")
        target = base
        counter = 2
        while target.exists():
            target = base.with_name(f"{base.name}_{counter}")
            counter += 1
        try:
            outputs = export_phase_stripping_session(self.session, target)
        except Exception as exc:
            self._show_warning(f"导出失败：{exc}")
            return False
        self._last_export_path = target
        self.status_var.set(f"已显式导出：{outputs['json']}")
        messagebox.showinfo("导出完成", f"残差 CSV、图像和会话 JSON 已保存到：\n{target}", parent=self)
        return True

    def _show_warning(self, text: str) -> None:
        self.warning_var.set(text)

    def _request_close(self) -> None:
        if self._closing:
            self.destroy()
            return
        if self.session.accepted_operations:
            answer = messagebox.askyesnocancel(
                "关闭残差谱分析",
                "当前有已接受的物相操作。关闭前是否显式导出残差谱会话？\n选择“否”将只丢弃本窗口会话，不写入永久结果。",
                parent=self,
            )
            if answer is None:
                return
            if answer and not self._export():
                return
        self._closing = True
        self.destroy()

    def destroy(self) -> None:
        existed = bool(self.winfo_exists())
        super().destroy()
        if existed and self._on_closed is not None:
            callback, self._on_closed = self._on_closed, None
            callback()


def open_or_raise_phase_stripping_window(
    parent: tk.Misc,
    existing: PhaseStrippingWindow | None,
    *,
    context: AnalysisContext,
    candidates: Iterable[Candidate],
    element_scope: Iterable[str] = (),
    on_closed: Callable[[], None] | None = None,
) -> PhaseStrippingWindow:
    if existing is not None:
        try:
            if existing.winfo_exists():
                existing.deiconify()
                existing.lift()
                existing.focus_force()
                return existing
        except tk.TclError:
            pass
    return PhaseStrippingWindow(
        parent,
        context,
        candidates,
        element_scope=element_scope,
        on_closed=on_closed,
    )
