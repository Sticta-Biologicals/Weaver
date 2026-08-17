# Sanger Verification

Sanger verification compares one or more sequencing reads with a readable reference plasmid sequence and saves the parsed files, alignment evidence, metrics, and review decision as a verification run.

## Upload a Run

Open a plasmid and choose `Align > Sanger`, or use `/inventory/plasmid/align/sanger/<plasmid_id>`. The plasmid must have a readable reference sequence. Upload one or more `.ab1`, `.phd.1`, or `.seq` files. The form also accepts FASTA extensions for the separate FASTA alignment mode, but FASTA files cannot be mixed with Sanger trace files in one batch. At least one supported file is required.

Optionally enter a run label and notes, and enable saving a Clustal file. The default processing limit is 96 files, 8 MB per file, and 96 MB for the batch. Invalid extensions, an empty upload, mixed alignment modes, oversized files, or unreadable content are reported as validation or processing errors.

Weaver groups files by base name when a read has more than one representation. For a group it prefers AB1, then PHD.1, then SEQ as the sequence source while retaining the uploaded files. SEQ supplies sequence text only; it cannot supply chromatogram traces or Phred quality values.

## Read Orientation and Alignment

For each read, Weaver records the selected source, raw and trimmed sequence, parsing errors, warnings, quality metrics, and detected orientation. Orientation can be `Forward`, `Reverse complement`, `Ambiguous`, or `Unmapped`. The aligned read is compared with the plasmid reference, and the run aggregates coverage, alignment metrics, variants, gaps, and low-confidence regions. A read can be retained but marked unusable when its sequence or quality data do not support alignment.

Open a saved run at `/inventory/plasmid/align/sanger/<plasmid_id>/run/<run_id>`. The alignment browser provides a reference coordinate view, read coverage and orientation, feature context when the reference is annotated, variants, gaps, and quality regions. Use its read controls to move through the reference and inspect the evidence supporting a possible variant. Coordinates for stored variants are zero-based internally; use the displayed reference coordinate labels when communicating a finding.

## Quality, Coverage, Variants, and Gaps

Coverage shows which reference positions are supported by usable reads. A gap is a reference interval without sufficient aligned read support. A variant is an observed difference from the reference and may carry a quality value and a low-quality flag. Low-confidence and intermediate-confidence regions identify positions where the base call deserves manual review rather than automatic acceptance. Compare a variant with all covering reads and the chromatogram before deciding whether it is biological, a sequencing artifact, or unresolved.

## Chromatograms

For an AB1-backed read, choose its chromatogram action at `/inventory/plasmid/align/sanger/<plasmid_id>/run/<run_id>/read/<read_id>/chromatogram`. The viewer plots A, C, G, and T traces with base calls and quality coloring. Use `Standard` or `Autoscaled`, the left and right controls, the zoom slider, and the base number plus `Go` to inspect a region. Low and intermediate confidence regions are shaded. PHD.1 and SEQ reads do not provide the AB1 trace arrays needed for this viewer.

## Manual Decision and Validation State

The automated run classification can be `PASS`, `REVIEW`, `FAIL`, or `NO_DATA`. It is evidence for review, not the final laboratory decision. On the saved run, choose a manual decision of `Verified`, `Not verified`, `Inconclusive`, or `Pending`; enter an effective date and comment when appropriate, then save.

When a run is marked `Verified`, Weaver sets the associated plasmid's sequencing validation state to verified and stores the manual effective date as the plasmid sequencing date. Other decisions do not mark the plasmid as verified and clear the sequencing date when the saved decision is not verified. The run retains the reviewer, review timestamp, decision, effective date, comment, automated reasons, and any saved Clustal file.

## FASTA Alignment

For sequence text without trace data, choose `Align > Fasta` at `/inventory/plasmid/align/fasta/<plasmid_id>`. Paste FASTA text or upload one or more `.fa`, `.fas`, `.fasta`, or `.txt` files. Choose `Together` or `One at a time` in `View mode`, and optionally save a Clustal file. This mode provides sequence alignment evidence but no chromatogram or Phred interpretation.
