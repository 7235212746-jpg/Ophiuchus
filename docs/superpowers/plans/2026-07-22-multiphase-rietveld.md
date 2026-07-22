# Multiphase Rietveld Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add guarded two-to-four-phase RIETAN-FP refinement with experimental wt%, repeat-run stability checks, and explicit suppression of untrustworthy quantitative results.

**Architecture:** Keep the existing single-phase backend intact. Add a multiphase domain layer, reproduce the official `cif2ins + combins` input contract in a focused builder, run RIETAN-FP as the sole final-profile engine, and place a quantification gate between parsed engine results and the UI. Ophi's existing fast simulator remains a screening aid only.

**Tech Stack:** Python 3.11+, dataclasses, NumPy, SciPy, tkinter/ttk, matplotlib, VESTA, RIETAN-FP 3.12, unittest.

## Global Constraints

- Support Cu Kalpha1+2 only in the first release.
- Accept two to four CIF-backed phases; phase 1 is the selected target.
- Keep atom positions, occupancies, displacement parameters, preferred orientation, and anisotropic broadening locked.
- Do not derive wt% from peak area, candidate score, or the existing impurity signal-share estimator.
- Hide exact wt% whenever the quantification gate returns `FAIL`.
- Label non-failing quantitative output `实验性定量` until a measured known-composition standard is registered.
- Use operating-system temporary directories and persist RIETAN artifacts only on explicit export.
- Preserve the existing single-phase workflow and public interfaces.

---

## File Map

- Create `ophiuchus/xrd/multiphase_models.py`: immutable multiphase inputs and results.
- Create `ophiuchus/xrd/multiphase_input.py`: official-template validation, `cif2ins` invocation, and deterministic phase combination.
- Create `ophiuchus/xrd/multiphase_refinement.py`: RIETAN execution, output parsing, staged repeated runs, and provenance.
- Create `ophiuchus/xrd/quantification.py`: Hill-Howard weight calculation and `PASS/WARNING/FAIL` gate.
- Create `tests/test_multiphase_models.py`: model and weight-fraction tests.
- Create `tests/test_multiphase_input.py`: official marker and combined-input tests.
- Create `tests/test_multiphase_refinement.py`: parser and backend contract tests.
- Create `tests/test_quantification_gate.py`: hard-failure and stability tests.
- Modify `ophiuchus/refinement/window.py`: phase selection, multiphase execution, contribution plot, and guarded table.
- Modify `tests/test_refinement_window.py`: UI state and result rendering tests.
- Modify `ophiuchus/xrd/rietan_backend.py`: discover `cif2ins.exe` and the official multiphase template.
- Modify `README.md` and `docs/Ophiuchus_操作手册.md`: installation, workflow, interpretation, and citation requirements.

### Task 1: Domain Models And Weight Fractions

**Files:**
- Create: `ophiuchus/xrd/multiphase_models.py`
- Create: `ophiuchus/xrd/quantification.py`
- Test: `tests/test_multiphase_models.py`

**Interfaces:**
- Produces: `PhaseRefinementInput`, `MultiphaseRefinementSettings`, `PhaseRefinementResult`, `MultiphaseRefinementResult`, `GateLevel`, `GateFinding`, `QuantificationAssessment`, `weight_fractions_from_zmv()`.

- [ ] **Step 1: Write failing tests for phase validation and Hill-Howard normalization**

```python
def test_weight_fractions_use_scale_times_zmv():
    values = weight_fractions_from_zmv([
        ("main", 2.0, 4.0, 100.0, 50.0),
        ("imp", 1.0, 2.0, 50.0, 25.0),
    ])
    self.assertAlmostEqual(values["main"], 94.117647, places=6)
    self.assertAlmostEqual(values["imp"], 5.882353, places=6)

def test_phase_input_requires_a_real_cif(tmp_path):
    with self.assertRaises(FileNotFoundError):
        PhaseRefinementInput("x", "X", self.root / "missing.cif", "target")
```

- [ ] **Step 2: Run `python -m unittest tests.test_multiphase_models -v` and verify missing imports fail**

- [ ] **Step 3: Implement immutable models and a pure, finite-positive ZMV calculation**

```python
def weight_fractions_from_zmv(rows):
    weighted = {phase_id: scale * z * molar_mass * volume for phase_id, scale, z, molar_mass, volume in rows}
    if not weighted or any(not np.isfinite(value) or value <= 0 for value in weighted.values()):
        raise ValueError("Every quantitative phase requires positive finite scale, Z, M, and V values.")
    total = sum(weighted.values())
    return {phase_id: value / total * 100.0 for phase_id, value in weighted.items()}
```

- [ ] **Step 4: Run the focused tests, then `python -m unittest tests.test_rietan_refinement -v`**

- [ ] **Step 5: Commit `Add multiphase refinement domain models`**

### Task 2: Official Multiphase Input Contract

**Files:**
- Create: `ophiuchus/xrd/multiphase_input.py`
- Modify: `ophiuchus/xrd/rietan_backend.py`
- Create: `tests/test_multiphase_input.py`
- Create: `tests/fixtures/rietan_multiphase/template_excerpt.ins`
- Create: `tests/fixtures/rietan_multiphase/phase1.ins`
- Create: `tests/fixtures/rietan_multiphase/phase2.ins`

**Interfaces:**
- Consumes: `PhaseRefinementInput`.
- Produces: `discover_cif2ins_executable()`, `discover_multiphase_template()`, `export_phase_input()`, `combine_phase_inputs(template_text, phase_texts)`.

- [ ] **Step 1: Write failing tests for required official markers and deterministic phase order**

```python
def test_combines_numbered_official_phase_sections():
    combined = combine_phase_inputs(TEMPLATE, [PHASE1, PHASE2])
    assert "NPHASE@ = 2:" in combined
    assert combined.index("! Phase @1") < combined.index("! Phase @2")
    assert combined.index("! Parameters @1") < combined.index("! Parameters @2")
    assert combined.count("Zr1@1/Zr") == 1
    assert combined.count("Sn1@2/Sn") == 1

def test_rejects_a_vesta_template_without_combins_markers():
    with self.assertRaisesRegex(ValueError, "multiphase template"):
        combine_phase_inputs("NMODE = 0", [PHASE1, PHASE2])
```

- [ ] **Step 2: Run the focused test and verify it fails because the builder is absent**

- [ ] **Step 3: Implement official marker extraction without ad-hoc phase renumbering**

```python
PHASE_RE = re.compile(r"(?ms)^\s*!\s*Phase @\d+.*?^\s*# End Phase @\d+\s*$")
PARAM_RE = re.compile(r"(?ms)^\s*!\s*Parameters @\d+.*?^\s*# End Parameters @\d+\s*$")
ELEMENT_RE = re.compile(r"(?ms)^\s*!\s*Elements @\d+\s*$\n(.*?)^\s*# End Elements @\d+\s*$")
```

The implementation follows the official `combins.command` order: template prefix, all phase blocks, template bridge, all parameter blocks, constraints for phases 2..N, template suffix, merged unique element symbols, and final `NPHASE@` replacement.

- [ ] **Step 4: Implement `cif2ins.exe` invocation with ASCII temporary names `phase@1.cif` through `phase@4.cif`**

```python
command = [
    str(cif2ins_exe), "0", cif.name, template.name, ins.name,
    "report.tex", "report.pdf", "structure.pdf", "result.lst", "mscs.pdf", "density.pdf",
]
```

- [ ] **Step 5: Run focused tests and compare the generated structure with the official Cu3Fe4P6 combins example contract**

- [ ] **Step 6: Commit `Build official RIETAN multiphase inputs`**

### Task 3: Multiphase RIETAN Backend And Parser

**Files:**
- Create: `ophiuchus/xrd/multiphase_refinement.py`
- Create: `tests/test_multiphase_refinement.py`

**Interfaces:**
- Consumes: `PhaseRefinementInput`, `MultiphaseRefinementSettings`, input-builder functions, `build_rietan_command()`, `write_xy_intensity()`.
- Produces: `parse_multiphase_output(gpd_text, lst_text, phase_inputs)`, `RietanMultiphaseBackend.refine(phases, x, intensity, settings)`.

- [ ] **Step 1: Write failing parser tests with two repeated scale factors and phase summaries**

```python
def test_parser_assigns_scale_and_zmv_by_phase_order():
    result = parse_multiphase_output(GPD, LST, PHASES)
    self.assertAlmostEqual(result.phases[0].scale, 1.25e-5)
    self.assertAlmostEqual(result.phases[1].scale, 3.50e-6)
    self.assertAlmostEqual(result.phases[0].z, 3.0)
    self.assertAlmostEqual(result.phases[1].volume_angstrom3, 410.647)
```

- [ ] **Step 2: Run the parser test and verify the missing parser failure**

- [ ] **Step 3: Implement parser state machines for final statistics, per-phase scale/ZMV, reflection groups, and profile columns**

Do not associate phases by formula text alone. Use RIETAN phase number first and cross-check formula/CIF provenance after parsing.

- [ ] **Step 4: Write a failing backend contract test that asserts VESTA is not used for final profile calculation**

```python
with patch("ophiuchus.xrd.multiphase_refinement.subprocess.run", side_effect=fake_run) as run:
    result = backend.refine(PHASES, x, y, settings)
assert Path(run.call_args_list[-1].args[0][0]).name.lower() == "rietan.exe"
assert result.provenance["profile_engine"] == "RIETAN-FP multiphase refinement"
```

- [ ] **Step 5: Implement a temporary-directory backend with timeout, cancellation-ready process boundaries, and no permanent writes**

- [ ] **Step 6: Run focused parser/backend tests and existing RIETAN tests**

- [ ] **Step 7: Commit `Run guarded RIETAN multiphase refinement`**

### Task 4: Stability And Quantification Gate

**Files:**
- Modify: `ophiuchus/xrd/quantification.py`
- Create: `tests/test_quantification_gate.py`

**Interfaces:**
- Consumes: one or more `MultiphaseRefinementResult` values.
- Produces: `assess_quantification(primary, repeated_results, validation_states)`.

- [ ] **Step 1: Write failing tests for hard failures**

```python
def test_negative_scale_hides_weight_fraction():
    assessment = assess_quantification(result_with(scale=-1.0), [], all_validated())
    assert assessment.level is GateLevel.FAIL
    assert assessment.allow_weight_percent is False

def test_unexplained_residual_group_is_a_hard_failure():
    assessment = assess_quantification(result_with_residual_group(45.0, 46.0), [], all_validated())
    assert assessment.level is GateLevel.FAIL
```

- [ ] **Step 2: Run focused tests and verify expected failures**

- [ ] **Step 3: Implement explicit findings for engine failure, invalid ZMV, unexplained peak groups, non-identifiable phases, and validation failure**

- [ ] **Step 4: Write failing repeat-run stability tests using multiple initial scale vectors**

```python
def test_unstable_repeat_runs_fail_quantification():
    assessment = assess_quantification(primary_70_30(), [result_55_45(), result_82_18()], all_validated())
    assert assessment.level is GateLevel.FAIL
    assert any("initial" in item.code for item in assessment.findings)
```

- [ ] **Step 5: Implement median, range, and maximum absolute wt% spread reporting**

The backend runs at least three deterministic initial scale vectors: equal, target-dominant, and impurity-dominant. The first release does not invent a universal scientific threshold; the configured threshold is marked as a software guard and included in provenance.

- [ ] **Step 6: Run all pure scientific tests**

- [ ] **Step 7: Commit `Gate experimental phase quantification`**

### Task 5: Refinement Window Integration

**Files:**
- Modify: `ophiuchus/refinement/window.py`
- Modify: `tests/test_refinement_window.py`

**Interfaces:**
- Consumes: `RietanMultiphaseBackend`, `QuantificationAssessment`, existing candidate lists and `AnalysisContext`.
- Produces: user-selectable target plus up to three CIF-backed impurities, staged run controls, contribution plot, guarded wt% table, explicit export.

- [ ] **Step 1: Write failing UI tests for default phase selection and disabled invalid candidates**

```python
def test_multiphase_mode_selects_target_and_top_three_cif_impurities():
    window = RefinementWindow(root, context, candidates, multiphase_backend=FakeMultiphaseBackend())
    assert window.target_var.get().startswith("ZrFe6Ge4")
    assert len(window.selected_impurity_ids()) <= 3
    assert all(Path(item.source_path).suffix.lower() == ".cif" for item in window.selected_phases())
```

- [ ] **Step 2: Run the UI test and verify the missing controls failure**

- [ ] **Step 3: Add a `单相确认 / 多相实验定量` segmented mode and a scrollable phase checklist**

Keep target and impurity roles explicit. Reject duplicate structure hashes and more than four total phases before starting a worker thread.

- [ ] **Step 4: Write failing result-rendering tests for hidden and visible wt% states**

```python
def test_failed_gate_never_renders_numeric_weight_percent():
    window.show_multiphase_result(result, failed_assessment())
    assert "不可定量" in window.metrics_var.get()
    assert "32.4%" not in window.metrics_var.get()
```

- [ ] **Step 5: Render observed, calculated, background, difference, per-phase contribution lines, and phase reflection ticks**

- [ ] **Step 6: Add explicit export of PNG, CSV, JSON, CIF hashes, engine versions, settings, gate findings, and stability runs**

- [ ] **Step 7: Run refinement-window and app-experience tests**

- [ ] **Step 8: Commit `Expose guarded multiphase Rietveld workflow`**

### Task 6: Real-Engine Validation, Documentation, And Release Gate

**Files:**
- Create: `tests/test_rietan_multiphase_integration.py`
- Modify: `README.md`
- Modify: `docs/Ophiuchus_操作手册.md`
- Modify: `docs/version_history/CHANGELOG.md`

**Interfaces:**
- Consumes: installed official `cif2ins.exe`, multiphase template, RIETAN.exe, Zr3V3GeSn4 CIF, and ZrFe6Ge4 CIF.
- Produces: skipped-when-unavailable integration tests and a documented user workflow.

- [ ] **Step 1: Add an opt-in real-engine test that builds and runs a two-phase official-example job**

```python
@unittest.skipUnless(os.environ.get("OPHI_RUN_RIETAN_INTEGRATION") == "1", "real RIETAN integration disabled")
def test_official_two_phase_job_runs_without_format_errors(self):
    result = backend.refine(phases, x, y, settings)
    self.assertEqual(len(result.phases), 2)
    self.assertTrue(np.isfinite(result.rwp_percent))
```

- [ ] **Step 2: Run the test without the environment flag and confirm it skips cleanly**

- [ ] **Step 3: Run the real-engine test with `OPHI_RUN_RIETAN_INTEGRATION=1` on the local installation**

- [ ] **Step 4: Add 60-90 degree reflection regression checks for Zr3V3GeSn4 and ZrFe6Ge4**

- [ ] **Step 5: Document official RIETAN asset discovery, experimental wt% meaning, non-crystalline limitations, and required citations**

- [ ] **Step 6: Run `run_tests.bat`, then start Ophiuchus and manually exercise target selection, multiphase execution, failed-gate hiding, and explicit export**

- [ ] **Step 7: Run `git diff --check` and inspect `git status --short`**

- [ ] **Step 8: Commit `Validate multiphase Rietveld workflow`**
