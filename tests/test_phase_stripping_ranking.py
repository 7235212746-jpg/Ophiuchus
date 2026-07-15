import unittest

import numpy as np

from ophiuchus.phase_stripping.models import AnalysisContext, FitBounds
from ophiuchus.phase_stripping.profile import project_candidate_profile
from ophiuchus.phase_stripping.ranking import deduplicate_candidates, rank_candidates
from ophiuchus.phase_stripping.session import PhaseStrippingSession
from ophiuchus.xrd.models import Candidate, Peak


def make_candidate(
    candidate_id: str,
    peaks: list[Peak],
    *,
    formula: str = "FeGe",
    structure_hash: str | None = None,
    space_group: str = "P 1",
    phase_entry_ids: list[str] | None = None,
    elements: list[str] | None = None,
) -> Candidate:
    return Candidate(
        candidate_id=candidate_id,
        formula_pretty=formula,
        source="test",
        source_path=f"{candidate_id}.cif",
        elements=elements or ["Fe", "Ge"],
        structure_hash=structure_hash,
        theory_peaks=peaks,
        phase_entry_ids=phase_entry_ids or [],
        space_group_symbol=space_group,
    )


def make_session(
    x: np.ndarray,
    intensity: np.ndarray,
    *,
    tolerance_deg: float = 0.15,
    bounds: FitBounds | None = None,
) -> PhaseStrippingSession:
    context = AnalysisContext(
        x=x,
        intensity=intensity,
        radiation="CuKalpha12",
        wavelength_angstrom=1.54056,
        two_theta_range=(float(x[0]), float(x[-1])),
        tolerance_deg=tolerance_deg,
        source_path="sample.xy",
        source_fingerprint="source-sha256",
        data_fingerprint="data-sha256",
    )
    return PhaseStrippingSession(
        context,
        bounds=bounds or FitBounds(shift_deg=(-0.01, 0.01), sigma_deg=(0.049, 0.051), scale=(0.0, 2.0)),
    )


class PhaseStrippingRankingTests(unittest.TestCase):
    def test_accidental_one_peak_candidate_ranks_below_complete_multi_peak_candidate(self):
        x = np.arange(20.0, 29.001, 0.002)
        complete = make_candidate(
            "complete",
            [Peak("a", 21.0, 100.0), Peak("b", 24.0, 75.0), Peak("c", 27.0, 60.0)],
        )
        accidental = make_candidate("accidental", [Peak("one", 21.0, 100.0)])
        intensity = project_candidate_profile(x, complete, shift_deg=0.0, sigma_deg=0.05)

        ranked = rank_candidates(make_session(x, intensity), [accidental, complete])

        self.assertEqual(ranked[0].candidate.candidate_id, "complete")
        isolated = next(item for item in ranked if item.candidate.candidate_id == "accidental")
        self.assertLess(isolated.final_score, 0.4)
        self.assertTrue(any("isolated" in warning.lower() for warning in isolated.warnings))

    def test_missing_strong_theoretical_peaks_warn_and_reduce_score(self):
        x = np.arange(20.0, 27.001, 0.002)
        incomplete = make_candidate(
            "incomplete",
            [Peak("a", 21.0, 100.0), Peak("b", 24.0, 90.0), Peak("c", 26.0, 80.0)],
        )
        observed_only = make_candidate("observed", [Peak("a", 21.0, 100.0)])
        intensity = project_candidate_profile(x, observed_only, shift_deg=0.0, sigma_deg=0.05)

        evidence = rank_candidates(make_session(x, intensity), [incomplete])[0]

        self.assertGreater(evidence.missing_strong_penalty, 0.0)
        self.assertLess(evidence.strong_peak_coverage, 0.5)
        self.assertTrue(any("Missing strong" in warning for warning in evidence.warnings))

    def test_candidate_evidence_explanations_are_readable_chinese_ui_copy(self):
        x = np.arange(20.0, 25.001, 0.002)
        candidate = make_candidate(
            "candidate",
            [Peak("a", 21.0, 100.0), Peak("b", 24.0, 65.0)],
        )
        intensity = project_candidate_profile(x, candidate, 0.0, 0.05)

        evidence = rank_candidates(make_session(x, intensity), [candidate])[0]

        self.assertIn("残差离散峰得分", evidence.explanations[0])
        self.assertIn("强峰覆盖率", evidence.explanations[1])
        self.assertIn("全谱改善", evidence.explanations[2])

    def test_deduplication_collapses_equivalent_profiles_but_preserves_polymorphs_and_provenance(self):
        peaks = [Peak("a", 21.0, 100.0), Peak("b", 24.0, 60.0)]
        duplicate_one = make_candidate("c1", peaks, structure_hash="same", phase_entry_ids=["mp-1"])
        duplicate_two = make_candidate("c2", peaks, structure_hash="same", phase_entry_ids=["mp-2"])
        duplicate_profile = make_candidate("c3", peaks, structure_hash="other", phase_entry_ids=["mp-3"])
        polymorph = make_candidate(
            "c4",
            [Peak("a", 21.4, 100.0), Peak("b", 24.5, 60.0)],
            structure_hash="polymorph",
            space_group="P 21/c",
            phase_entry_ids=["mp-4"],
        )

        collapsed = deduplicate_candidates([duplicate_one, duplicate_two, duplicate_profile, polymorph])
        evidence = rank_candidates(
            make_session(np.arange(20.0, 26.001, 0.002), np.zeros(3001)),
            [duplicate_one, duplicate_two, duplicate_profile, polymorph],
        )

        self.assertEqual([candidate.candidate_id for candidate in collapsed], ["c1", "c4"])
        self.assertEqual(collapsed[0].phase_entry_ids, ["mp-1", "mp-2", "mp-3"])
        collapsed_evidence = next(item for item in evidence if item.candidate.candidate_id == "c1")
        self.assertEqual(
            collapsed_evidence.provenance_ids,
            ("c1", "c2", "c3", "mp-1", "mp-2", "mp-3"),
        )

    def test_ranking_reacts_to_the_current_residual_after_accepting_phase_a(self):
        x = np.arange(20.0, 29.001, 0.002)
        phase_a = make_candidate("phase-a", [Peak("a1", 21.0, 100.0), Peak("a2", 23.0, 65.0)])
        phase_b = make_candidate("phase-b", [Peak("b1", 26.0, 100.0), Peak("b2", 28.0, 70.0)])
        intensity = project_candidate_profile(x, phase_a, 0.0, 0.05)
        intensity += 0.8 * project_candidate_profile(x, phase_b, 0.0, 0.05)
        session = make_session(x, intensity)

        before = session.rank_candidates([phase_a, phase_b])
        session.accept_preview(session.preview(phase_a))
        after = session.rank_candidates([phase_a, phase_b])

        self.assertEqual(before[0].candidate.candidate_id, "phase-a")
        self.assertEqual(after[0].candidate.candidate_id, "phase-b")
        self.assertGreater(
            next(item for item in after if item.candidate.candidate_id == "phase-b").final_score,
            next(item for item in after if item.candidate.candidate_id == "phase-a").final_score,
        )

    def test_independent_evidence_excludes_peaks_explained_by_an_accepted_phase(self):
        x = np.arange(20.0, 29.001, 0.002)
        phase_a = make_candidate("phase-a", [Peak("shared-a", 21.0, 100.0), Peak("a-only", 22.5, 50.0)])
        phase_b = make_candidate(
            "phase-b",
            [Peak("shared-b", 21.0, 100.0), Peak("b-one", 25.0, 80.0), Peak("b-two", 27.0, 60.0)],
        )
        shared_only = make_candidate("shared-only", [Peak("shared", 21.0, 100.0)])
        intensity = project_candidate_profile(x, phase_a, 0.0, 0.05)
        intensity += 0.8 * project_candidate_profile(x, phase_b, 0.0, 0.05)
        session = make_session(x, intensity)

        session.accept_preview(session.preview(phase_a))
        ranked = session.rank_candidates([phase_a, phase_b, shared_only])
        phase_b_evidence = next(item for item in ranked if item.candidate.candidate_id == "phase-b")
        shared_only_evidence = next(item for item in ranked if item.candidate.candidate_id == "shared-only")

        self.assertAlmostEqual(phase_b_evidence.independent_peak_evidence, 140.0 / 240.0, places=3)
        self.assertEqual(shared_only_evidence.independent_peak_evidence, 0.0)
        self.assertGreater(phase_b_evidence.independent_peak_evidence, shared_only_evidence.independent_peak_evidence)

    def test_deduplication_compares_against_all_group_members_regardless_of_input_order(self):
        first_profile = [Peak("first-a", 21.0, 100.0), Peak("first-b", 24.0, 60.0)]
        second_profile = [Peak("second-a", 30.0, 100.0), Peak("second-b", 33.0, 60.0)]
        first = make_candidate("first", first_profile, structure_hash="shared-hash", phase_entry_ids=["entry-1"])
        bridge = make_candidate("bridge", second_profile, structure_hash="shared-hash", phase_entry_ids=["entry-2"])
        profile_duplicate = make_candidate("profile-duplicate", second_profile, structure_hash="other-hash", phase_entry_ids=["entry-3"])

        forward = deduplicate_candidates([first, bridge, profile_duplicate])
        reverse = deduplicate_candidates([bridge, profile_duplicate, first])

        self.assertEqual(len(forward), 1)
        self.assertEqual(len(reverse), 1)
        self.assertEqual(set(forward[0].phase_entry_ids), {"entry-1", "entry-2", "entry-3"})
        self.assertEqual(set(reverse[0].phase_entry_ids), {"entry-1", "entry-2", "entry-3"})

    def test_element_scope_rejects_candidates_with_elements_outside_the_allowed_range(self):
        x = np.arange(20.0, 25.001, 0.002)
        outside_scope = make_candidate(
            "outside-scope",
            [Peak("a", 21.0, 100.0), Peak("b", 23.0, 65.0)],
            elements=["Fe", "O"],
        )
        intensity = project_candidate_profile(x, outside_scope, 0.0, 0.05)

        evidence = rank_candidates(make_session(x, intensity), [outside_scope], element_scope=["Fe", "Ge"])[0]

        self.assertEqual(evidence.final_score, 0.0)
        scope_warnings = [warning for warning in evidence.warnings if "outside current element scope" in warning]
        self.assertEqual(scope_warnings, ["Candidate contains elements outside current element scope: O."])

    def test_low_intensity_patterns_use_all_positive_peaks_as_the_effective_strong_set(self):
        x = np.arange(20.0, 25.001, 0.002)
        low_intensity = make_candidate("low-intensity", [Peak("a", 21.0, 5.0), Peak("b", 24.0, 4.0)])
        observed_only = make_candidate("observed", [Peak("a", 21.0, 5.0)])
        intensity = project_candidate_profile(x, observed_only, 0.0, 0.05)

        evidence = rank_candidates(make_session(x, intensity), [low_intensity])[0]

        self.assertAlmostEqual(evidence.strong_peak_coverage, 5.0 / 9.0, places=3)
        self.assertAlmostEqual(evidence.missing_strong_penalty, 0.7 * 4.0 / 9.0, places=3)
        self.assertTrue(any("Missing strong" in warning for warning in evidence.warnings))

    def test_independent_evidence_uses_shifted_accepted_contributions_when_accepted_candidate_is_omitted(self):
        x = np.arange(20.0, 30.001, 0.002)
        phase_a = make_candidate("phase-a", [Peak("a-shared", 21.0, 100.0), Peak("a-only", 23.0, 50.0)])
        phase_b = make_candidate(
            "phase-b",
            [Peak("b-shared", 21.12, 100.0), Peak("b-one", 26.0, 80.0), Peak("b-two", 28.0, 60.0)],
        )
        shared_only = make_candidate("shared-only", [Peak("shared", 21.12, 100.0)])
        intensity = project_candidate_profile(x, phase_a, 0.12, 0.05)
        intensity += 0.8 * project_candidate_profile(x, phase_b, 0.0, 0.05)
        session = make_session(
            x,
            intensity,
            tolerance_deg=0.04,
            bounds=FitBounds(shift_deg=(0.11, 0.13), sigma_deg=(0.049, 0.051), scale=(0.0, 2.0)),
        )

        session.accept_preview(session.preview(phase_a))
        ranked = rank_candidates(session, [phase_b, shared_only])
        phase_b_evidence = next(item for item in ranked if item.candidate.candidate_id == "phase-b")
        shared_only_evidence = next(item for item in ranked if item.candidate.candidate_id == "shared-only")

        self.assertTrue(any(abs(position - 21.12) < 0.001 for position in session.accepted_peak_positions))
        self.assertAlmostEqual(phase_b_evidence.independent_peak_evidence, 140.0 / 240.0, places=3)
        self.assertEqual(shared_only_evidence.independent_peak_evidence, 0.0)


if __name__ == "__main__":
    unittest.main()
