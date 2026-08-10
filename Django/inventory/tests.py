from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

from django.core.files.uploadedfile import SimpleUploadedFile
from django.contrib.auth.models import User
from django.test import RequestFactory
from django.test import SimpleTestCase
from django.test import TestCase
from django.urls import reverse
from Bio.Seq import Seq

from inventory.custom.primer_access import visible_primers_for_user
from inventory.custom.primer_dimers import PrimerDimerConditions
from inventory.custom.primer_dimers import PrimerDimerInputError
from inventory.custom.primer_dimers import PrimerDimerThresholds
from inventory.custom.primer_dimers import PrimerInput
from inventory.custom.primer_dimers import analyze_pair
from inventory.custom.primer_dimers import analyze_primers
from inventory.custom.primer_dimers import primer_combinations
from inventory.custom.primer_dimers import read_primers
from inventory.custom.primer_dimers import validate_primers
from inventory.custom.primer_import import import_primers_from_fasta
from inventory.custom.primer_import import primer_entries_from_fasta
from inventory.custom.pcr import classify_ytk_overhang
from inventory.custom.pcr import find_primer_binding_hits
from inventory.custom.pcr import infer_type_iis_overhang
from inventory.custom.pcr import inferred_primer_parts
from inventory.custom.pcr import matching_amplicon_annotations
from inventory.custom.pcr import primer_pair_amplicons
from inventory.custom.pcr import primer_pair_complementarity
from inventory.custom.pcr import select_non_overlapping_amplicons
from inventory.custom.pcr import suggest_pcr_primers
from inventory.custom.restriction_digest import DigestConstraints
from inventory.custom.restriction_digest import LabEnzyme
from inventory.custom.restriction_digest import compatible_buffers
from inventory.custom.restriction_digest import compatible_temperature
from inventory.custom.restriction_digest import digest_candidates
from inventory.custom.restriction_digest import enzymes_with_effective_cuts
from inventory.custom.restriction_digest import effective_cut_sites
from inventory.custom.restriction_digest import evaluate_digest
from inventory.custom.restriction_digest import fragment_sizes
from inventory.custom.restriction_digest import min_band_difference
from inventory.custom.restriction_digest import normalize_regions
from inventory.custom.restriction_digest import region_contains_position
from inventory.custom.restriction_digest import serialize_digest_response
from inventory.custom.sanger import detect_format
from inventory.custom.sanger import detect_confidence_regions
from inventory.custom.sanger import normalized_group_name
from inventory.custom.sanger import parse_phd1
from inventory.custom.sanger import parse_seq
from inventory.custom.sanger import process_sanger_files
from inventory.models import Primer
from inventory.models import Plasmid
from inventory.models import SangerVerificationRun
from inventory.views import fasta_alignment_result
from inventory.views import fasta_record_from_text
from inventory.views import fasta_records_from_text
from inventory.views import amplicon_contains_region
from inventory.views import amplicon_matches_any_primer_id
from inventory.views import amplicon_matches_primer_id
from inventory.views import optional_int_query_param
from inventory.views import sanger_feature_color
from organization.models import Membership
from organization.models import Project


def primer(name, sequence_3, direction):
    return SimpleNamespace(
        id=uuid4(),
        name=name,
        sequence_3=sequence_3,
        sequence_5="",
        fwd_or_rev=direction,
    )


def primer_with_overhang(name, sequence_3, sequence_5, direction):
    return SimpleNamespace(
        id=uuid4(),
        name=name,
        sequence_3=sequence_3,
        sequence_5=sequence_5,
        fwd_or_rev=direction,
    )


class RestrictionEnzymeStub(SimpleNamespace):
    def __str__(self):
        return self.name


def restriction_enzyme(name, **activities):
    defaults = {
        "activity_buffer_1_1": 100,
        "activity_buffer_2_1": 100,
        "activity_buffer_3_1": 100,
        "activity_buffer_CS": 100,
        "activity_buffer_aari": 100,
    }
    defaults.update(activities)
    return RestrictionEnzymeStub(name=name, hf_version=False, **defaults)


def lab_enzyme(name, temperature=37, activities=None):
    return LabEnzyme(
        name=name,
        display_name=name,
        recognition_site="TEST",
        fcut=1,
        rcut=1,
        temperature=temperature,
        activities=activities or {
            "buffer_1_1": 100,
            "buffer_2_1": 100,
            "buffer_3_1": 100,
            "buffer_CS": 100,
            "buffer_aari": 100,
        },
    )


class FastaParsingTests(SimpleTestCase):
    def test_fasta_records_accept_comment_lines_before_header(self):
        records = fasta_records_from_text("; exported from aligner\n# comment\n>read1\nACGT\n")

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "read1")
        self.assertEqual(str(records[0].seq), "ACGT")

    def test_fasta_record_from_text_wraps_plain_sequence(self):
        record = fasta_record_from_text("ACGTACGT", name="Amplicon")

        self.assertEqual(record.id, "Amplicon")
        self.assertEqual(str(record.seq), "ACGTACGT")

    def test_fasta_alignment_result_uses_sanger_result_shape(self):
        record = fasta_record_from_text(">query\nGTACGT\n")
        result = fasta_alignment_result("AACGTACGTT", record)

        self.assertEqual(result["combined"]["useful_reads"], 1)
        self.assertEqual(result["reads"][0]["formats"], ["FASTA"])
        self.assertEqual(result["reads"][0]["alignment"]["start_display"], 4)
        self.assertIn("reference_projection", result["reads"][0]["alignment"])

    def test_fasta_alignment_result_accepts_multiple_records_and_detects_reverse(self):
        records = fasta_records_from_text(">forward\nACGTTGC\n>reverse\nGCAACGT\n")
        result = fasta_alignment_result("AAAACGTTGCAAAA", records)

        self.assertEqual(result["combined"]["useful_reads"], 2)
        self.assertEqual(result["reads"][0]["alignment"]["best_orientation"], "forward")
        self.assertEqual(result["reads"][1]["alignment"]["best_orientation"], "reverse")


class SangerVerificationServiceTests(SimpleTestCase):
    def test_normalizes_complementary_sanger_files_to_one_group(self):
        base = "sample_C02"

        self.assertEqual(normalized_group_name(base + ".ab1"), base)
        self.assertEqual(normalized_group_name(base + ".phd.1"), base)
        self.assertEqual(normalized_group_name(base + ".seq"), base)
        self.assertEqual(detect_format(base + ".phd.1"), "phd1")

    def test_parse_seq_accepts_plain_or_fasta_sequence(self):
        parsed = parse_seq(b">read one\nacgt n\nRYSW\n")

        self.assertFalse(parsed.errors)
        self.assertEqual(parsed.sequence, "ACGTNRYSW")

    def test_parse_seq_reports_invalid_characters(self):
        parsed = parse_seq(b"ACGTZ\n")

        self.assertTrue(parsed.errors)
        self.assertIn("Z", parsed.errors[0])

    def test_parse_phd1_reads_bases_qualities_and_peaks(self):
        text = b"""BEGIN_SEQUENCE read1\r\nBEGIN_COMMENT\r\nCHROMAT_FILE: read1.ab1\r\nEND_COMMENT\r\nBEGIN_DNA\r\nA 31 10\r\nC 20 22\r\nEND_DNA\r\nEND_SEQUENCE\r\n"""

        parsed = parse_phd1(text)

        self.assertFalse(parsed.errors)
        self.assertEqual(parsed.sequence, "AC")
        self.assertEqual(parsed.qualities, [31, 20])
        self.assertEqual(parsed.peak_positions, [10, 22])
        self.assertEqual(parsed.metadata["CHROMAT_FILE"], "read1.ab1")

    def test_duplicate_same_format_in_group_is_rejected(self):
        files = [
            SimpleUploadedFile("sample.seq", b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            SimpleUploadedFile("sample.seq", b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ]

        result = process_sanger_files(files, "ACGT" * 30)

        self.assertEqual(len(result["reads"]), 1)
        self.assertFalse(result["reads"][0]["is_usable"])
        self.assertIn("Multiple SEQ files", result["reads"][0]["errors"][0])

    def test_detects_reverse_complement_orientation(self):
        reference = "ATGCGTACCTAGGATCCGATTAACCGGTTAGCTAGGCTAATCGGATCCGGAATTCGCGGCCGC"
        reverse_read = str(Seq(reference[8:58]).reverse_complement())
        files = [SimpleUploadedFile("reverse.seq", reverse_read.encode())]

        result = process_sanger_files(files, reference)

        self.assertTrue(result["reads"][0]["is_usable"])
        self.assertEqual(result["reads"][0]["alignment"]["best_orientation"], "reverse")

    def test_partial_sanger_coverage_does_not_prevent_pass_by_itself(self):
        reference = "ACGT" * 500

        from inventory.custom.sanger import SangerProcessingParameters
        from inventory.custom.sanger import align_read
        from inventory.custom.sanger import classify_run
        from inventory.custom.sanger import combined_metrics

        alignment = align_read(reference, "NNNNNN" + reference[100:600], [40] * 506, 0, SangerProcessingParameters())
        reads = [{
            "name": "partial",
            "is_usable": True,
            "alignment": alignment,
            "warnings": [],
            "errors": [],
        }]
        combined = combined_metrics(reference, reads, SangerProcessingParameters())
        classification = classify_run(combined, reads, SangerProcessingParameters())

        self.assertLess(combined["combined_coverage"], 95)
        self.assertEqual(classification["state"], "PASS")
        self.assertEqual(alignment["query_start"], 6)
        self.assertEqual(alignment["reference_projection_base_indices"][alignment["start"]], 6)

    def test_seq_without_quality_reduces_confidence_without_high_quality_variants(self):
        reference = "ACGT" * 80
        read_sequence = reference[:40] + "T" + reference[41:80]
        files = [SimpleUploadedFile("no-quality.seq", read_sequence.encode())]

        result = process_sanger_files(files, reference)

        self.assertEqual(result["combined"]["high_quality_variant_count"], 0)
        self.assertEqual(result["classification"]["state"], "REVIEW")
        self.assertIn("Quality scores unavailable", result["classification"]["reasons"][0])

    def test_confidence_regions_detect_progressive_3prime_drop(self):
        sequence = "A" * 120
        qualities = [32] * 80 + [18] * 10 + [9] * 30

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertTrue(metrics["low_confidence_regions"])
        self.assertGreaterEqual(metrics["low_confidence_regions"][0]["start"], 75)
        self.assertEqual(metrics["low_confidence_regions"][-1]["end"], 119)

    def test_terminal_quality_collapse_backtracks_to_first_nearby_warning(self):
        sequence = "A" * 1200
        qualities = [32] * 1200
        qualities[977] = 15
        qualities[1056:] = [8] * (1200 - 1056)

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(metrics["low_confidence_regions"][-1]["start"], 977)
        self.assertEqual(metrics["low_confidence_regions"][-1]["end"], 1199)
        self.assertEqual(metrics["alignment_blocks"][-1][1], 977)
        self.assertIn("Terminal quality collapse", " ".join(metrics["low_confidence_regions"][-1]["reasons"]))

    def test_confidence_regions_detect_bad_5prime_start(self):
        sequence = "A" * 120
        qualities = [8] * 25 + [31] * 95

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(metrics["low_confidence_regions"][0]["start"], 0)
        self.assertGreater(metrics["accepted_blocks"][0][0], 0)

    def test_confidence_regions_ignore_single_bad_base_inside_good_read(self):
        sequence = "A" * 100
        qualities = [32] * 100
        qualities[50] = 3

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(metrics["low_confidence_regions"], [])

    def test_confidence_regions_merge_bad_regions_separated_by_short_gap(self):
        sequence = "A" * 130
        qualities = [32] * 130
        qualities[30:45] = [8] * 15
        qualities[52:67] = [8] * 15

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(len(metrics["low_confidence_regions"]), 1)
        self.assertLessEqual(metrics["low_confidence_regions"][0]["start"], 30)
        self.assertGreaterEqual(metrics["low_confidence_regions"][0]["end"], 66)

    def test_confidence_regions_short_recovery_does_not_split_region(self):
        sequence = "A" * 150
        qualities = [32] * 150
        qualities[30:60] = [8] * 30
        qualities[60:70] = [32] * 10
        qualities[70:100] = [8] * 30

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(len(metrics["low_confidence_regions"]), 1)

    def test_confidence_regions_sustained_recovery_closes_region(self):
        sequence = "A" * 180
        qualities = [32] * 180
        qualities[25:45] = [8] * 20
        qualities[80:100] = [8] * 20

        metrics = detect_confidence_regions(sequence, qualities)

        self.assertEqual(len(metrics["low_confidence_regions"]), 2)

    def test_high_phred_low_signal_is_warning_not_low_region(self):
        sequence = "A" * 100
        qualities = [32] * 100
        base_pos = [i * 10 for i in range(100)]
        trace_len = base_pos[-1] + 20
        a_trace = [100] * trace_len
        for index in range(35, 50):
            a_trace[base_pos[index]] = 10
        chromatogram = {
            "basePos": base_pos,
            "aTrace": a_trace,
            "cTrace": [1] * trace_len,
            "gTrace": [1] * trace_len,
            "tTrace": [1] * trace_len,
        }

        metrics = detect_confidence_regions(sequence, qualities, chromatogram)

        self.assertEqual(metrics["low_confidence_regions"], [])
        self.assertTrue(metrics["signal"]["warnings"])

    def test_secondary_peaks_with_low_phred_support_low_region(self):
        sequence = "A" * 100
        qualities = [32] * 100
        qualities[35:55] = [10] * 20
        base_pos = [i * 10 for i in range(100)]
        trace_len = base_pos[-1] + 20
        a_trace = [100] * trace_len
        c_trace = [1] * trace_len
        for index in range(35, 55):
            c_trace[base_pos[index]] = 70
        chromatogram = {
            "basePos": base_pos,
            "aTrace": a_trace,
            "cTrace": c_trace,
            "gTrace": [1] * trace_len,
            "tTrace": [1] * trace_len,
        }

        metrics = detect_confidence_regions(sequence, qualities, chromatogram)

        self.assertTrue(metrics["low_confidence_regions"])
        self.assertIn("Chromatogram signal morphology", metrics["low_confidence_regions"][0]["reasons"])

    def test_clear_range_intersects_weaver_alignment_blocks(self):
        sequence = "A" * 120
        qualities = [32] * 120
        qualities[:20] = [8] * 20
        qualities[100:] = [8] * 20

        metrics = detect_confidence_regions(sequence, qualities, clear_range=(10, 90))

        self.assertEqual(metrics["file_clear_range"], (10, 90))
        self.assertGreaterEqual(metrics["alignment_blocks"][0][0], 10)
        self.assertLessEqual(metrics["alignment_blocks"][-1][1], 90)

    def test_seq_without_quality_reports_quality_unavailable_and_aligns_with_warning(self):
        reference = "ACGT" * 40
        files = [SimpleUploadedFile("plain.seq", reference[:100].encode())]

        result = process_sanger_files(files, reference)

        self.assertFalse(result["reads"][0]["quality_metrics"]["quality_available"])
        self.assertTrue(result["reads"][0]["is_usable"])
        self.assertIn("Quality scores unavailable", result["reads"][0]["warnings"][0])

    def test_internal_low_quality_region_splits_read_into_aligned_blocks(self):
        reference = "ACGT" * 90
        read = reference[:150] + reference[180:330]
        qualities = [35] * 150 + [35] * 150
        sequence = read[:80] + ("N" * 40) + read[120:]
        qualities[80:120] = [3] * 40
        phd_lines = ["BEGIN_SEQUENCE split", "BEGIN_COMMENT", "END_COMMENT", "BEGIN_DNA"]
        phd_lines.extend("{} {} {}".format(base, quality, index * 10) for index, (base, quality) in enumerate(zip(sequence, qualities)))
        phd_lines.extend(["END_DNA", "END_SEQUENCE"])

        result = process_sanger_files([SimpleUploadedFile("split.phd.1", "\n".join(phd_lines).encode())], reference)

        alignment = result["reads"][0]["alignment"]
        self.assertTrue(result["reads"][0]["quality_metrics"]["low_confidence_regions"])
        self.assertGreaterEqual(len(alignment["segments"]), 2)
        self.assertEqual(result["combined"]["high_quality_variant_count"], 0)

    def test_multiple_forward_and_reverse_reads_combine_high_confidence_coverage(self):
        reference = "ATGCGTAC" * 80
        reverse_read = str(Seq(reference[300:520]).reverse_complement())
        files = [
            SimpleUploadedFile("left.seq", reference[:220].encode()),
            SimpleUploadedFile("right.seq", reverse_read.encode()),
        ]

        result = process_sanger_files(files, reference)

        self.assertEqual(result["combined"]["useful_reads"], 2)
        self.assertIn("reverse", {read["alignment"]["best_orientation"] for read in result["reads"]})
        self.assertGreater(result["combined"]["combined_coverage"], 30)

    def test_real_complementary_files_are_grouped_when_available(self):
        base = Path("/home/tproschle-sticta/Downloads/PUC-6.jul.26/Diego Lagos 060726/6096-Diego-495_P6_Stuffer_c1-457-pL0-chech-F_2026-07-07-13-15-37_C02")
        paths = [Path(str(base) + suffix) for suffix in (".ab1", ".phd.1", ".seq")]
        if not all(path.exists() for path in paths):
            self.skipTest("Real Sanger fixture files are not available on this machine")
        files = [SimpleUploadedFile(path.name, path.read_bytes()) for path in paths]
        reference = parse_seq(paths[2].read_bytes()).sequence

        result = process_sanger_files(files, reference)

        self.assertEqual(len(result["reads"]), 1)
        self.assertIn("ab1", result["reads"][0]["formats"])
        self.assertIn("phd1", result["reads"][0]["formats"])
        self.assertIn("seq", result["reads"][0]["formats"])


class SangerFeatureColorTests(SimpleTestCase):
    def test_known_feature_roles_use_stable_palette(self):
        self.assertEqual(sanger_feature_color("CDS", "CamR", {"color": ["#0000ff"]}), "#2fb344")
        self.assertEqual(sanger_feature_color("misc_feature", "pTDH3 promoter", {}), "#f0b429")
        self.assertEqual(sanger_feature_color("terminator", "CYC1 terminator", {}), "#d94841")
        self.assertEqual(sanger_feature_color("restriction_site", "BsaI", {}), "#8f63d9")
        self.assertEqual(sanger_feature_color("protein_bind", "BsaI", {}), "#8f63d9")
        self.assertEqual(sanger_feature_color("misc_feature", "attB1 recombination site", {}), "#e66fb2")
        self.assertEqual(sanger_feature_color("misc_recomb", "Hr1Chr4_Up", {}), "#e66fb2")
        self.assertEqual(sanger_feature_color("misc_recomb", "HR1-Chr4", {}), "#e66fb2")

    def test_unknown_feature_keeps_explicit_genbank_color(self):
        self.assertEqual(sanger_feature_color("misc_feature", "custom tag", {"ApEinfo_fwdcolor": ["#123456"]}), "#123456")


class PcrSuggestionTests(SimpleTestCase):
    def test_infers_type_iis_overhang_from_full_primer_sequence(self):
        inferred = infer_type_iis_overhang("aaCGTCTCtctccTATGcgtaaaggcgaagag")

        self.assertEqual(inferred["sequence_5"], "AACGTCTCTCTCCTATG")
        self.assertEqual(inferred["sequence_3"], "CGTAAAGGCGAAGAG")
        self.assertEqual(inferred["cloning_overhang"], "TATG")
        self.assertEqual(inferred["ytk"]["key"], "2-3")
        self.assertEqual(inferred["ytk"]["orientation"], "forward")

    def test_infers_reverse_primer_type_iis_overhang_from_full_sequence(self):
        inferred = infer_type_iis_overhang("aaCGTCTCtctcgGGATtttgtacagttcatccataccatg")

        self.assertEqual(inferred["sequence_5"], "AACGTCTCTCTCGGGAT")
        self.assertEqual(inferred["sequence_3"], "TTTGTACAGTTCATCCATACCATG")
        self.assertEqual(inferred["cloning_overhang"], "GGAT")
        self.assertEqual(inferred["ytk"]["key"], "3-4")
        self.assertEqual(inferred["ytk"]["canonical_overhang"], "ATCC")
        self.assertEqual(inferred["ytk"]["orientation"], "reverse_complement")

    def test_classifies_non_ytk_overhang_as_empty(self):
        self.assertEqual(classify_ytk_overhang("AAAA")["key"], "")

    def test_amplicon_finder_uses_inferred_hybridizing_sequence(self):
        sequence = (
            "CGTAAAGGCGAAGAG" +
            "AAAAC" +
            str(Seq("TTTGTACAGTTCATCCATACCATG").reverse_complement())
        )
        primers = [
            primer("693-L0-P3-sfGFP-F", "aaCGTCTCtctccTATGcgtaaaggcgaagag", "f"),
            primer("694-L0-P3-sfGFP-R", "aaCGTCTCtctcgGGATtttgtacagttcatccataccatg", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_product_size=999,
            max_tm_difference=99,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["notes"]["fwd_primer"], ["L0-P3-sfGFP-F"])
        self.assertEqual(amplicons[0]["notes"]["rev_primer"], ["L0-P3-sfGFP-R"])
        self.assertEqual(
            amplicons[0]["notes"]["amplicon_sequence"],
            ["AACGTCTCTCTCCTATG" + sequence + str(Seq("AACGTCTCTCTCGGGAT").reverse_complement())],
        )

    def test_existing_declared_overhang_is_not_reinferred(self):
        parts = inferred_primer_parts(
            primer_with_overhang("declared-F", "aaCGTCTCtctccTATGcgtaaaggcgaagag", "TT", "f")
        )

        self.assertFalse(parts["inferred"])
        self.assertEqual(parts["sequence_5"], "TT")
        self.assertEqual(parts["sequence_3"], "AACGTCTCTCTCCTATGCGTAAAGGCGAAGAG")

    def test_primer_binding_hits_can_match_from_3prime_end(self):
        hits = find_primer_binding_hits(
            "AAAACCCCGGGGTTTT",
            "GGGGGAAAACCCCGGGGTTTT",
            min_binding_length=15,
        )

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["start"], 0)
        self.assertEqual(hits[0]["binding_sequence"], "AAAACCCCGGGGTTTT")
        self.assertEqual(hits[0]["unmatched_5"], "GGGGG")

    def test_amplicon_finder_allows_3prime_binding_without_full_primer_containment(self):
        fwd_binding = "AAAACCCCGGGGTTTT"
        rev_binding = "TTTTGGGGCCCCAAAA"
        sequence = fwd_binding + "GG" + str(Seq(rev_binding).reverse_complement())
        primers = [
            primer_with_overhang("partial-F", "GGGGG" + fwd_binding, "AA", "f"),
            primer_with_overhang("partial-R", "CCCCC" + rev_binding, "TT", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_product_size=999,
            max_tm_difference=99,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["notes"]["partial_primer_binding"], ["true"])
        self.assertEqual(amplicons[0]["notes"]["fwd_partial_binding"], ["true"])
        self.assertEqual(amplicons[0]["notes"]["rev_partial_binding"], ["true"])
        self.assertEqual(amplicons[0]["notes"]["fwd_binding_length"], ["16"])
        self.assertEqual(amplicons[0]["notes"]["rev_binding_length"], ["16"])
        self.assertEqual(amplicons[0]["notes"]["fwd_unmatched_5"], ["GGGGG"])
        self.assertEqual(amplicons[0]["notes"]["rev_unmatched_5"], ["CCCCC"])
        self.assertIn("FWD partial 3' binding: 16 bp aligned, 5 bp 5' unaligned", amplicons[0]["notes"]["warnings"])
        self.assertIn("REV partial 3' binding: 16 bp aligned, 5 bp 5' unaligned", amplicons[0]["notes"]["warnings"])
        self.assertEqual(
            amplicons[0]["notes"]["amplicon_sequence"],
            ["AA" + "GGGGG" + sequence + str(Seq("CCCCC").reverse_complement()) + "AA"],
        )

    def test_primer_pair_amplicons_flags_partial_binding(self):
        fwd_binding = "AAAACCCCGGGGTTTT"
        rev_binding = "TTTTGGGGCCCCAAAA"
        sequence = fwd_binding + "GG" + str(Seq(rev_binding).reverse_complement())
        primer_f = primer("partial-F", "GGGGG" + fwd_binding, "f")
        primer_r = primer("partial-R", "CCCCC" + rev_binding, "r")

        amplicons = primer_pair_amplicons(
            sequence,
            primer_f,
            primer_r,
            max_product_size=999,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertTrue(amplicons[0]["partial_primer_binding"])
        self.assertTrue(amplicons[0]["fwd_partial_binding"])
        self.assertTrue(amplicons[0]["rev_partial_binding"])

    def test_amplicon_region_filter_requires_region_inside_amplicon(self):
        amplicon = {"start": 10, "end": 30}
        inside = SimpleNamespace(start=12, end=20)
        outside = SimpleNamespace(start=5, end=12)

        self.assertTrue(amplicon_contains_region(amplicon, inside, 100))
        self.assertFalse(amplicon_contains_region(amplicon, outside, 100))

    def test_amplicon_region_filter_allows_configured_flank(self):
        amplicon = {"start": 6374, "end": 6844}
        region = SimpleNamespace(start=6369, end=6743)

        self.assertFalse(amplicon_contains_region(amplicon, region, 9171))
        self.assertTrue(amplicon_contains_region(amplicon, region, 9171, flank_bp=30))

    def test_amplicon_region_filter_supports_circular_amplicons(self):
        amplicon = {"start": 90, "end": 10, "overlapsSelf": True}
        inside = SimpleNamespace(start=95, end=5)
        outside = SimpleNamespace(start=20, end=30)

        self.assertTrue(amplicon_contains_region(amplicon, inside, 100))
        self.assertFalse(amplicon_contains_region(amplicon, outside, 100))

    def test_amplicon_primer_id_filter_matches_either_primer(self):
        amplicon = {"notes": {"fwd_primer_id": ["693"], "rev_primer_id": ["694"]}}

        self.assertTrue(amplicon_matches_primer_id(amplicon, "693"))
        self.assertTrue(amplicon_matches_primer_id(amplicon, "694"))
        self.assertFalse(amplicon_matches_primer_id(amplicon, "695"))
        self.assertTrue(amplicon_matches_any_primer_id(amplicon, ("695", "694")))
        self.assertFalse(amplicon_matches_any_primer_id(amplicon, ("695", "696")))

    def test_amplicon_size_filter_limits_are_optional(self):
        factory = RequestFactory()

        empty_request = factory.get("/amplicons", {"min_size": "", "max_size": ""})
        min_request = factory.get("/amplicons", {"min_size": "200", "max_size": ""})
        max_request = factory.get("/amplicons", {"min_size": "", "max_size": "1200"})

        self.assertEqual(optional_int_query_param(empty_request, "min_size", 100), 100)
        self.assertIsNone(optional_int_query_param(empty_request, "max_size"))
        self.assertEqual(optional_int_query_param(min_request, "min_size", 100), 200)
        self.assertEqual(optional_int_query_param(max_request, "max_size"), 1200)

    def test_suggests_pair_covering_selection(self):
        sequence = "AAAACCCCGGGGTTTT"
        primers = [
            primer("1-left-F", "AAAA", "f"),
            primer("2-right-R", "CCCC", "r"),
        ]

        suggestions = suggest_pcr_primers(
            sequence,
            primers,
            4,
            7,
            margin=8,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(len(suggestions), 1)
        self.assertEqual(suggestions[0]["category"], "Full selection")
        self.assertEqual(suggestions[0]["selection_coverage"], 4)
        self.assertEqual(suggestions[0]["product_range"], "1..12")

    def test_suggests_circular_pair_crossing_origin(self):
        sequence = "AAACCCGGGTTT"
        primers = [
            primer("1-origin-F", "TTT", "f"),
            primer("2-origin-R", "GGG", "r"),
        ]

        suggestions = suggest_pcr_primers(
            sequence,
            primers,
            10,
            2,
            margin=5,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(len(suggestions), 1)
        self.assertIn("Circular", suggestions[0]["category"])
        self.assertEqual(suggestions[0]["selection_coverage"], 5)
        self.assertEqual(suggestions[0]["product_range"], "10..6 (circular)")

    def test_default_minimum_product_size_excludes_short_amplicons(self):
        sequence = "AAAACCCCGGGGTTTT"
        primers = [
            primer("1-left-F", "AAAA", "f"),
            primer("2-right-R", "CCCC", "r"),
        ]

        suggestions = suggest_pcr_primers(sequence, primers, 4, 7, margin=8)

        self.assertEqual(suggestions, [])

    def test_default_tm_filter_excludes_dissimilar_pairs(self):
        sequence = "GGGGAAAACCCCTTTTGGGGCCCC"
        primers = [
            primer("1-low-tm-F", "AAAA", "f"),
            primer("2-high-tm-R", "GGGG", "r"),
        ]

        suggestions = suggest_pcr_primers(sequence, primers, 4, 7, margin=20, min_product_size=1)

        self.assertEqual(suggestions, [])

    def test_matching_amplicons_are_returned_as_parts(self):
        sequence = "AAAACCCCGGGGTTTT"
        primers = [
            primer("1-left-F", "AAAA", "f"),
            primer("2-right-R", "CCCC", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["annotationTypePlural"], "parts")
        self.assertTrue(amplicons[0]["id"].startswith("weaver-amplicon-"))
        self.assertIn("tm_difference", amplicons[0]["notes"])
        self.assertEqual(amplicons[0]["notes"]["fwd_tm"], ["8.0"])
        self.assertEqual(amplicons[0]["notes"]["rev_tm"], ["16.0"])
        self.assertEqual(amplicons[0]["notes"]["product_tm"], ["40.0"])
        self.assertEqual(amplicons[0]["notes"]["recommended_annealing_tm"], ["15.5"])

    def test_matching_amplicons_include_copyable_amplicon_sequence(self):
        sequence = "AAAACCCCGGGGTTTT"
        primers = [
            primer_with_overhang("1-left-F", "AAAA", "TT", "f"),
            primer_with_overhang("2-right-R", "CCCC", "GG", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(amplicons[0]["notes"]["amplicon_sequence"], ["TTAAAACCCCGGGGCC"])
        self.assertEqual(amplicons[0]["notes"]["product_size"], ["16"])
        self.assertIn("primer_complementarity_alignment_fwd", amplicons[0]["notes"])
        self.assertIn("primer3_dimer_risk", amplicons[0]["notes"])
        self.assertIn("primer3_dimer_dg", amplicons[0]["notes"])

    def test_circular_amplicon_visual_interval_matches_real_product(self):
        sequence = "AAACCCGGGTTT"
        primers = [
            primer("1-origin-F", "TTT", "f"),
            primer("2-origin-R", "GGG", "r"),
        ]

        amplicons = matching_amplicon_annotations(
            sequence,
            primers,
            min_product_size=1,
            max_tm_difference=99,
        )

        self.assertEqual(len(amplicons), 1)
        self.assertTrue(amplicons[0]["overlapsSelf"])
        self.assertEqual(amplicons[0]["notes"]["template_size"], ["9"])
        self.assertEqual(amplicons[0]["notes"]["visual_start"], ["9"])
        self.assertEqual(amplicons[0]["notes"]["visual_end"], ["5"])
        self.assertEqual(amplicons[0]["notes"]["visual_size"], ["9"])
        self.assertEqual(amplicons[0]["notes"]["visual_is_alternative"], ["false"])

    def test_selects_non_overlapping_amplicons(self):
        amplicons = [
            {"id": "a", "start": 0, "end": 9},
            {"id": "b", "start": 5, "end": 15},
            {"id": "c", "start": 16, "end": 20},
        ]

        selected = select_non_overlapping_amplicons(amplicons, 30)

        self.assertEqual([amplicon["id"] for amplicon in selected], ["a", "c"])

    def test_selects_non_overlapping_amplicons_across_origin(self):
        amplicons = [
            {"id": "a", "start": 18, "end": 3, "overlapsSelf": True},
            {"id": "b", "start": 2, "end": 5},
            {"id": "c", "start": 8, "end": 10},
        ]

        selected = select_non_overlapping_amplicons(amplicons, 20)

        self.assertEqual([amplicon["id"] for amplicon in selected], ["a", "c"])

    def test_primer_pair_amplicons_finds_product(self):
        sequence = "AAAACCCCGGGGTTTT"
        primer_f = primer("1-left-F", "AAAA", "f")
        primer_r = primer("2-right-R", "CCCC", "r")

        amplicons = primer_pair_amplicons(sequence, primer_f, primer_r)

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["product_range"], "1..12")
        self.assertEqual(amplicons[0]["product_size"], 12)
        self.assertFalse(amplicons[0]["circular"])

    def test_primer_pair_amplicons_finds_circular_product(self):
        sequence = "AAACCCGGGTTT"
        primer_f = primer("1-origin-F", "TTT", "f")
        primer_r = primer("2-origin-R", "GGG", "r")

        amplicons = primer_pair_amplicons(sequence, primer_f, primer_r)

        self.assertEqual(len(amplicons), 1)
        self.assertEqual(amplicons[0]["product_range"], "10..6 (circular)")
        self.assertEqual(amplicons[0]["product_size"], 9)
        self.assertTrue(amplicons[0]["circular"])
        self.assertIn("primer3_dimer", amplicons[0])

    def test_primer_pair_complementarity_flags_3prime_dimer_risk(self):
        primer_f = primer("1-risk-F", "AAAACCCC", "f")
        primer_r = primer("2-risk-R", "GGGGTTTT", "r")

        complementarity = primer_pair_complementarity(primer_f, primer_r)

        self.assertEqual(complementarity["severity"], "high")
        self.assertEqual(complementarity["max_both_3prime_contiguous"], 8)
        self.assertEqual(set(complementarity["alignment"].keys()), {"fwd", "match", "rev"})
        self.assertIn("|", complementarity["alignment"]["match"])
        self.assertIn("FWD 5'", complementarity["alignment"]["fwd"])
        self.assertIn("REV 3'", complementarity["alignment"]["rev"])
        self.assertIn("AAAACCCC", complementarity["alignment"]["fwd"])
        self.assertIn("TTTTGGGG", complementarity["alignment"]["rev"])

    def test_primer_pair_complementarity_allows_unmatched_pair(self):
        primer_f = primer("1-ok-F", "AAAA", "f")
        primer_r = primer("2-ok-R", "AAAA", "r")

        complementarity = primer_pair_complementarity(primer_f, primer_r)

        self.assertEqual(complementarity["severity"], "none")
        self.assertEqual(complementarity["max_contiguous"], 0)
        self.assertEqual(set(complementarity["alignment"].keys()), {"fwd", "match", "rev"})
        self.assertIn("AAAA", complementarity["alignment"]["fwd"])
        self.assertIn("AAAA", complementarity["alignment"]["rev"])

    def test_primer_pair_complementarity_does_not_warn_for_trivial_matches(self):
        primer_f = primer("1-trivial-F", "GCCTTTTGCTGGCCTTTTGC", "f")
        primer_r = primer("2-trivial-R", "CTGTGTTGACATCTGGTTTG", "r")

        complementarity = primer_pair_complementarity(primer_f, primer_r)

        self.assertEqual(complementarity["severity"], "none")
        self.assertEqual(complementarity["max_contiguous"], 2)
        self.assertEqual(complementarity["max_both_3prime_contiguous"], 1)
        self.assertEqual(complementarity["warnings"], [])
        self.assertIn("||", complementarity["alignment"]["match"])


class RestrictionDigestTests(SimpleTestCase):
    def test_circular_plasmid_without_cuts_has_no_fragments(self):
        result = evaluate_digest(
            "AAAACCCCGGGG",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragment_size_bp=0, min_band_difference_bp=0),
            is_circular=True,
        )

        self.assertEqual(result["cut_count"], 0)
        self.assertEqual(result["fragment_count"], 0)
        self.assertFalse(result["exact"])

    def test_circular_plasmid_with_one_cut_returns_full_length_fragment(self):
        result = evaluate_digest(
            "AAAAGAATTCTTT",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0),
            is_circular=True,
        )

        self.assertEqual(result["fragments_map_order"], [13])

    def test_single_fragment_does_not_check_band_separation(self):
        result = evaluate_digest(
            "AAAAGAATTCTTT",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, max_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=500),
            is_circular=True,
        )

        self.assertTrue(result["exact"])
        self.assertFalse(any(item["criterion"] == "minimum band separation bp" for item in result["violations"]))

    def test_circular_plasmid_with_multiple_cuts_includes_origin_fragment(self):
        result = evaluate_digest(
            "AAAAGAATTCTTTGAATTCAAA",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=2, min_fragment_size_bp=1, min_band_difference_bp=1),
            is_circular=True,
        )

        self.assertEqual([site["position"] for site in result["cut_sites"]], [5, 14])
        self.assertEqual(result["fragments_map_order"], [9, 13])

    def test_linear_fragment_generation_zero_one_and_many_cuts(self):
        self.assertEqual(fragment_sizes(10, [], is_circular=False), [10])
        self.assertEqual(fragment_sizes(10, [4], is_circular=False), [4, 6])
        self.assertEqual(fragment_sizes(10, [2, 7], is_circular=False), [2, 5, 3])

    def test_coincident_cuts_keep_enzyme_identity(self):
        result = evaluate_digest(
            "AAAAGAATTCTTT",
            [lab_enzyme("EcoRI"), lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0),
            is_circular=True,
        )

        self.assertEqual(result["cut_count"], 1)
        self.assertEqual(result["cut_sites"][0]["enzymes"], ["EcoRI"])

    def test_enzyme_with_multiple_sites(self):
        sites = effective_cut_sites("AAAAGAATTCTTTGAATTCAAA", lab_enzyme("EcoRI"), is_circular=False)

        self.assertEqual([site["position"] for site in sites], [5, 14])

    def test_iupac_motif_site_is_found(self):
        sites = effective_cut_sites("AAAGAATCAAA", lab_enzyme("HinfI"), is_circular=False)

        self.assertEqual([site["position"] for site in sites], [4])

    def test_type_iis_cut_can_be_outside_recognition_site(self):
        sites = effective_cut_sites("AAAAGGTCTCAAAAATTTT", lab_enzyme("BsaI"), is_circular=False)

        self.assertEqual([site["position"] for site in sites], [11])

    def test_reverse_strand_type_iis_site_on_circular_sequence(self):
        sites = effective_cut_sites("AAAAGAGACCTTTTTAAAA", lab_enzyme("BsaI"), is_circular=True)

        self.assertEqual([site["position"] for site in sites], [18])

    def test_regions_include_normal_and_origin_crossing_intervals(self):
        normal, circular = normalize_regions([
            {"start": 2, "end": 5},
            {"start": 8, "end": 1},
        ], 10)

        self.assertTrue(region_contains_position(normal, 3, 10))
        self.assertFalse(region_contains_position(normal, 7, 10))
        self.assertTrue(region_contains_position(circular, 9, 10))
        self.assertTrue(region_contains_position(circular, 0, 10))
        self.assertFalse(region_contains_position(circular, 4, 10))

    def test_minimum_band_difference_at_threshold_is_accepted(self):
        self.assertEqual(min_band_difference([500, 1000, 1500]), 500)

    def test_equal_fragment_sizes_have_zero_band_separation(self):
        self.assertEqual(min_band_difference([1000, 1000]), 0)

    def test_fragment_count_above_range_is_not_exact(self):
        result = evaluate_digest(
            "AAAAGAATTCTTTGAATTCAAAGAATTCTTT",
            [lab_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, max_fragments=2, min_fragment_size_bp=1, min_band_difference_bp=0),
            is_circular=True,
        )

        self.assertEqual(result["fragment_count"], 3)
        self.assertFalse(result["exact"])
        self.assertIn("above range 1-2", result["violations"][0]["message"])

    def test_best_buffer_maximizes_minimum_then_average_activity(self):
        left = lab_enzyme("Left", activities={
            "buffer_1_1": 100,
            "buffer_2_1": 90,
            "buffer_3_1": 80,
            "buffer_CS": 80,
            "buffer_aari": 10,
        })
        right = lab_enzyme("Right", activities={
            "buffer_1_1": 80,
            "buffer_2_1": 90,
            "buffer_3_1": 100,
            "buffer_CS": 80,
            "buffer_aari": 10,
        })

        buffers = compatible_buffers([left, right], 75)

        self.assertEqual(buffers[0]["name"], "NEB 2.1")
        self.assertEqual(buffers[0]["min_activity"], 90)

    def test_incompatible_temperatures_are_not_exact(self):
        temperature, is_compatible = compatible_temperature([
            lab_enzyme("EcoRI", temperature=37),
            lab_enzyme("BamHI", temperature=42),
        ])

        self.assertEqual(temperature, 37)
        self.assertFalse(is_compatible)

    def test_unknown_activity_or_temperature_is_not_exact(self):
        enzyme = lab_enzyme("EcoRI", temperature=None, activities={
            "buffer_1_1": None,
            "buffer_2_1": None,
            "buffer_3_1": None,
            "buffer_CS": None,
            "buffer_aari": None,
        })

        result = evaluate_digest(
            "AAAAGAATTCTTT",
            [enzyme],
            DigestConstraints(min_fragments=0, min_fragment_size_bp=0, min_band_difference_bp=0),
        )

        self.assertEqual(result["temperature_status"], "unknown")
        self.assertFalse(result["exact"])
        self.assertTrue(any(item["criterion"] == "shared buffer" for item in result["violations"]))

    def test_exact_results_rank_before_closest_matches_deterministically(self):
        enzymes = [
            restriction_enzyme("EcoRI"),
            restriction_enzyme("HinfI"),
        ]

        results = digest_candidates(
            "AAAGAATCAAAAGAATTCAAA",
            enzymes,
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0, limit=10),
            is_circular=True,
        )

        self.assertTrue(results[0]["exact"])
        self.assertEqual(results, sorted(results, key=lambda result: (not result["exact"], len(result["enzymes"]), "+".join(result["enzymes"]))))

    def test_closest_match_reports_concrete_violations(self):
        results = digest_candidates(
            "AAAAGAATTCTTT",
            [restriction_enzyme("EcoRI")],
            DigestConstraints(min_fragments=2, min_fragment_size_bp=500, min_band_difference_bp=500),
            is_circular=True,
        )

        self.assertFalse(results[0]["exact"])
        self.assertTrue(results[0]["violations"])
        self.assertIn("missing", results[0]["violations"][0])

    def test_digest_candidates_skip_enzymes_without_effective_cuts(self):
        results = digest_candidates(
            "AAAAGAATTCTTT",
            [restriction_enzyme("EcoRI"), restriction_enzyme("BamHI")],
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0, limit=10),
            is_circular=True,
        )

        self.assertTrue(results)
        self.assertTrue(all("BamHI" not in result["enzyme_names"] for result in results))

    def test_enzymes_with_effective_cuts_filters_non_cutters(self):
        enzymes = enzymes_with_effective_cuts(
            "AAAAGAATTCTTT",
            [restriction_enzyme("EcoRI"), restriction_enzyme("BamHI")],
            is_circular=True,
        )

        self.assertEqual([enzyme.name for enzyme in enzymes], ["EcoRI"])

    def test_digest_candidates_can_require_an_enzyme(self):
        results = digest_candidates(
            "AAAGAATCAAAAGAATTCAAA",
            [restriction_enzyme("EcoRI"), restriction_enzyme("HinfI")],
            DigestConstraints(
                min_fragments=1,
                min_fragment_size_bp=1,
                min_band_difference_bp=0,
                limit=10,
                required_enzymes=("HinfI",),
            ),
            is_circular=True,
        )

        self.assertTrue(results)
        self.assertTrue(all("HinfI" in result["enzyme_names"] for result in results))

    def test_serialized_contract_contains_frontend_card_fields(self):
        response = serialize_digest_response(
            "AAAAGAATTCTTTGAATTCAAA",
            [restriction_enzyme("EcoRI")],
            DigestConstraints(min_fragments=1, min_fragment_size_bp=1, min_band_difference_bp=0),
            is_circular=True,
        )
        result = response["results"][0]

        for key in [
            "enzymes",
            "best_buffer",
            "cut_sites",
            "fragment_count",
            "fragments_map_order",
            "fragments_by_size",
            "min_band_difference",
            "regions",
            "status",
            "violations",
        ]:
            self.assertIn(key, result)


class PrimerBatchImportTests(TestCase):
    def test_defaults_missing_direction_to_forward(self):
        project = Project.objects.create(name="Sticta", public=False)
        fasta = StringIO(">123-no-direction\nAAAACCCC\n")

        result = import_primers_from_fasta(fasta, project, default_direction="f")

        imported = Primer.objects.get(project=project)
        self.assertEqual(result["created"], 1)
        self.assertEqual(imported.name, "123-no-direction")
        self.assertEqual(imported.sequence_3, "AAAACCCC")
        self.assertEqual(imported.fwd_or_rev, "f")

    def test_parser_infers_reverse_direction_from_name(self):
        entries = primer_entries_from_fasta(StringIO(">826_Primer-R\nAAAACCCC\n"), default_direction="f")

        self.assertEqual(entries[0]["direction"], "r")

    def test_fasta_header_overhang_is_split_from_full_sequence(self):
        entries = primer_entries_from_fasta(
            StringIO(">1001-L0-P3a-ATR-F overhang=aaCGTCTCtctcc\naaCGTCTCtctccTATGACTTCTGCTTTGTATGCATCAG\n")
        )

        self.assertEqual(entries[0]["name"], "1001-L0-P3a-ATR-F")
        self.assertEqual(entries[0]["direction"], "f")
        self.assertEqual(entries[0]["sequence_5"], "aaCGTCTCtctcc")
        self.assertEqual(entries[0]["sequence"], "TATGACTTCTGCTTTGTATGCATCAG")


class PrimerAccessTests(TestCase):
    def test_visible_primers_include_all_readable_projects(self):
        user = User.objects.create_user(username="weaver-user")
        visible_project_a = Project.objects.create(name="Visible A", public=False)
        visible_project_b = Project.objects.create(name="Visible B", public=False)
        hidden_project = Project.objects.create(name="Hidden", public=False)
        Membership.objects.create(member=user, project=visible_project_a, access_policies="r")
        Membership.objects.create(member=user, project=visible_project_b, access_policies="w")
        Primer.objects.create(name="1-visible-a-F", sequence_3="AAAA", fwd_or_rev="f", project=visible_project_a)
        Primer.objects.create(name="2-visible-b-R", sequence_3="CCCC", fwd_or_rev="r", project=visible_project_b)
        Primer.objects.create(name="3-hidden-F", sequence_3="GGGG", fwd_or_rev="f", project=hidden_project)

        primer_names = set(visible_primers_for_user(user).values_list("name", flat=True))

        self.assertEqual(primer_names, {"1-visible-a-F", "2-visible-b-R"})

    def test_fasta_header_direction_can_replace_name_suffix(self):
        entries = primer_entries_from_fasta(
            StringIO(">1001-L0-P3a-ATR overhang=aaCGTCTCtctcc direction=F\naaCGTCTCtctccTATGACTTCTGCTTTGTATGCATCAG\n"),
            require_direction=True,
        )

        self.assertEqual(entries[0]["name"], "1001-L0-P3a-ATR")
        self.assertEqual(entries[0]["direction"], "f")
        self.assertEqual(entries[0]["sequence_5"], "aaCGTCTCtctcc")
        self.assertEqual(entries[0]["sequence"], "TATGACTTCTGCTTTGTATGCATCAG")


class SangerVerificationEntryTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(username="seq-user", password="pw")
        self.project = Project.objects.create(name="Sticta", public=False)
        Membership.objects.create(member=self.user, project=self.project, access_policies="r")
        self.plasmid = Plasmid.objects.create(
            idx=499,
            name="L0-P8b-Hr1Chr4_Up",
            intended_use="test",
            project=self.project,
        )
        self.client.force_login(self.user)

    def test_entry_redirects_to_upload_when_no_sanger_run_exists(self):
        response = self.client.get(reverse("plasmid_seq_verification_entry", kwargs={"weaver_id": self.plasmid.idx}))

        self.assertRedirects(response, reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}))

    def test_entry_with_trailing_slash_redirects_to_upload(self):
        response = self.client.get(reverse("plasmid_seq_verification_entry_slash", kwargs={"weaver_id": self.plasmid.idx}))

        self.assertRedirects(response, reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}))

    def test_uuid_entry_redirects_to_upload_when_no_idx_link_is_available(self):
        response = self.client.get(reverse("plasmid_seq_verification_entry_uuid", kwargs={"plasmid_id": self.plasmid.id}))

        self.assertRedirects(response, reverse("plasmid_align_sanger", kwargs={"plasmid_id": self.plasmid.id}))

    def test_entry_redirects_to_latest_verified_run_first(self):
        SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)
        verified = SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user, manual_decision="VERIFIED")

        response = self.client.get(reverse("plasmid_seq_verification_entry", kwargs={"weaver_id": self.plasmid.idx}))

        self.assertRedirects(response, reverse("sanger_run_detail", kwargs={"plasmid_id": self.plasmid.id, "run_id": verified.id}), fetch_redirect_response=False)

    def test_entry_redirects_to_latest_loaded_run_when_none_verified(self):
        SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)
        latest = SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.user)

        response = self.client.get(reverse("plasmid_seq_verification_entry", kwargs={"weaver_id": self.plasmid.idx}))

        self.assertRedirects(response, reverse("sanger_run_detail", kwargs={"plasmid_id": self.plasmid.id, "run_id": latest.id}), fetch_redirect_response=False)


class SangerDecisionTests(TestCase):
    def setUp(self):
        self.editor = User.objects.create_user(username="seq-editor", password="pw")
        self.reader = User.objects.create_user(username="seq-reader", password="pw")
        self.project = Project.objects.create(name="DecisionProject", public=False)
        Membership.objects.create(member=self.editor, project=self.project, access_policies="w")
        Membership.objects.create(member=self.reader, project=self.project, access_policies="r")
        self.plasmid = Plasmid.objects.create(
            idx=501,
            name="Decision plasmid",
            intended_use="test",
            project=self.project,
            ligation_state=1,
            colonypcr_state=0,
            digestion_state=0,
            sequencing_state=1,
        )
        self.run = SangerVerificationRun.objects.create(plasmid=self.plasmid, created_by=self.editor)

    def test_verified_decision_updates_plasmid_and_evidence(self):
        self.client.force_login(self.editor)

        response = self.client.post(reverse("sanger_run_decision", kwargs={"plasmid_id": self.plasmid.id, "run_id": self.run.id}), {
            "manual_decision": "VERIFIED",
            "manual_decision_effective_date": "2026-07-30",
            "manual_decision_comment": "Clean overlapping reads.",
        })

        self.assertEqual(response.status_code, 302)
        self.run.refresh_from_db()
        self.plasmid.refresh_from_db()
        self.assertEqual(self.run.manual_decision, "VERIFIED")
        self.assertEqual(self.run.manual_decision_by, self.editor)
        self.assertEqual(str(self.run.manual_decision_effective_date), "2026-07-30")
        self.assertEqual(self.run.manual_decision_comment, "Clean overlapping reads.")
        self.assertEqual(self.plasmid.sequencing_state, 2)
        self.assertEqual(str(self.plasmid.sequencing_date), "2026-07-30")

    def test_inconclusive_decision_does_not_mark_plasmid_verified(self):
        self.client.force_login(self.editor)

        self.client.post(reverse("sanger_run_decision", kwargs={"plasmid_id": self.plasmid.id, "run_id": self.run.id}), {
            "manual_decision": "INCONCLUSIVE",
            "manual_decision_effective_date": "2026-07-31",
            "manual_decision_comment": "Coverage gap across junction.",
        })

        self.run.refresh_from_db()
        self.plasmid.refresh_from_db()
        self.assertEqual(self.run.manual_decision, "INCONCLUSIVE")
        self.assertEqual(self.plasmid.sequencing_state, 1)
        self.assertIsNone(self.plasmid.sequencing_date)

    def test_read_only_user_cannot_save_decision(self):
        self.client.force_login(self.reader)

        response = self.client.post(reverse("sanger_run_decision", kwargs={"plasmid_id": self.plasmid.id, "run_id": self.run.id}), {
            "manual_decision": "VERIFIED",
        })

        self.assertNotEqual(response.status_code, 302)
        self.run.refresh_from_db()
        self.plasmid.refresh_from_db()
        self.assertEqual(self.run.manual_decision, "")
        self.assertEqual(self.plasmid.sequencing_state, 1)


class PrimerDimerAnalysisTests(SimpleTestCase):
    def test_normalizes_lowercase_and_spaces(self):
        primers, warnings = validate_primers([
            PrimerInput("p1", "aa aa cccc\n"),
        ])

        self.assertEqual(warnings, [])
        self.assertEqual(primers[0].sequence, "AAAACCCC")

    def test_rejects_invalid_sequence(self):
        with self.assertRaises(PrimerDimerInputError):
            validate_primers([PrimerInput("p1", "AAAN")])

    def test_rejects_duplicate_names(self):
        with self.assertRaises(PrimerDimerInputError):
            validate_primers([PrimerInput("p1", "AAAA"), PrimerInput("p1", "CCCC")])

    def test_reports_identical_sequences_with_different_names(self):
        primers, warnings = validate_primers([
            PrimerInput("p1", "AAAA"),
            PrimerInput("p2", "AAAA"),
        ])

        self.assertEqual(len(primers), 2)
        self.assertEqual(len(warnings), 1)

    def test_multiplex_combination_count(self):
        primers = [PrimerInput("p1", "AAAA"), PrimerInput("p2", "CCCC"), PrimerInput("p3", "TTTT")]

        combinations = list(primer_combinations(primers))

        self.assertEqual(len(combinations), 6)

    def test_reads_fasta(self):
        primers = read_primers(StringIO(">p1\nAAAA\n>p2\nCCCC\n"), "fasta")

        self.assertEqual([primer.name for primer in primers], ["p1", "p2"])

    def test_reads_csv(self):
        primers = read_primers(StringIO("name,sequence\np1,AAAA\np2,CCCC\n"), "csv")

        self.assertEqual([primer.sequence for primer in primers], ["AAAA", "CCCC"])

    def test_high_for_both_3prime_extendable_heterodimer(self):
        conditions = PrimerDimerConditions(annealing_temp_c=55)
        thresholds = PrimerDimerThresholds(high_delta_g_kcal_mol=-999, moderate_delta_g_kcal_mol=-999)

        result = analyze_pair(
            PrimerInput("a", "AAAACCCC"),
            PrimerInput("b", "GGGGTTTT"),
            "heterodimer",
            conditions,
            thresholds,
        )

        self.assertEqual(result["risk"], "HIGH")
        self.assertTrue(result["extendable_3prime"])
        self.assertEqual(result["terminal_3prime_run"], 8)
        self.assertLess(result["delta_g_kcal_mol"], 0)

    def test_four_base_3prime_run_with_weak_tm_is_moderate_for_pcr(self):
        result = analyze_pair(
            PrimerInput("a", "GCCGGAGACCCAGGCGCGGC"),
            PrimerInput("b", "TCGCCGAAGCACCGGAGAGT"),
            "heterodimer",
            PrimerDimerConditions(annealing_temp_c=60),
            PrimerDimerThresholds(high_delta_g_kcal_mol=-999, moderate_delta_g_kcal_mol=-999),
        )

        self.assertEqual(result["risk"], "MODERATE")
        self.assertEqual(result["terminal_3prime_run"], 4)
        self.assertLess(result["dimer_tm_c"], 55)

    def test_category_is_symmetric_for_heterodimer(self):
        conditions = PrimerDimerConditions(annealing_temp_c=55)
        thresholds = PrimerDimerThresholds(high_delta_g_kcal_mol=-999, moderate_delta_g_kcal_mol=-999)
        primer_a = PrimerInput("a", "AAAACCCC")
        primer_b = PrimerInput("b", "GGGGTTTT")

        result_ab = analyze_pair(primer_a, primer_b, "heterodimer", conditions, thresholds)
        result_ba = analyze_pair(primer_b, primer_a, "heterodimer", conditions, thresholds)

        self.assertEqual(result_ab["risk"], result_ba["risk"])

    def test_conditions_change_thermodynamic_result_reproducibly(self):
        primer_a = PrimerInput("a", "AAAACCCC")
        primer_b = PrimerInput("b", "GGGGTTTT")
        low_mg = analyze_pair(
            primer_a,
            primer_b,
            "heterodimer",
            PrimerDimerConditions(dv_conc_mM=0.5, annealing_temp_c=55),
            PrimerDimerThresholds(),
        )
        high_mg = analyze_pair(
            primer_a,
            primer_b,
            "heterodimer",
            PrimerDimerConditions(dv_conc_mM=3.0, annealing_temp_c=55),
            PrimerDimerThresholds(),
        )

        self.assertNotEqual(low_mg["dimer_tm_c"], high_mg["dimer_tm_c"])

    def test_internal_run_contributes_to_risk_even_if_preview_prefers_3prime(self):
        result = analyze_pair(
            PrimerInput("a", "AGGAGATGACTTAAGGCA"),
            PrimerInput("b", "AAAATCTCCTGAATACAA"),
            "heterodimer",
            PrimerDimerConditions(annealing_temp_c=55),
            PrimerDimerThresholds(high_delta_g_kcal_mol=-999, moderate_delta_g_kcal_mol=-999),
        )

        self.assertEqual(result["risk"], "MODERATE")
        self.assertEqual(result["longest_complementary_run"], 7)

    def test_analyze_primers_returns_expected_number_of_rows(self):
        analysis = analyze_primers([
            PrimerInput("p1", "AAAA"),
            PrimerInput("p2", "CCCC"),
            PrimerInput("p3", "TTTT"),
        ], conditions=PrimerDimerConditions(annealing_temp_c=55))

        self.assertEqual(len(analysis["results"]), 6)
