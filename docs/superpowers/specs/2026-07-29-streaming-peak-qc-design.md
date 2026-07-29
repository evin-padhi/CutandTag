# Streaming Peak QC Design

## Problem

`PEAK_QC` is allocated 4 GB as a `process_light` task. Its CLI currently
reads the complete BEDPE file into memory and then creates additional
normalized and deduplicated fragment lists. Large CUT&Tag libraries can
therefore be killed with exit status 137 even though overlap calculation
itself only requires one fragment at a time.

## Design

Add a line-oriented fragment iterator and pass it through normalization,
name-based deduplication, and overlap counting without creating fragment
lists. Preserve the existing in-memory APIs for callers that provide lists,
and preserve the rule that repeated fragment names count once.

The peak index, per-peak counters, peak summaries, and output formats remain
unchanged. The deduplication set remains proportional to the number of unique
fragment names, but full `Fragment` objects and input lines are released as
the stream advances.

## Validation

Add a regression test with a guarded fragment generator that refuses to yield
its second record until overlap processing has occurred for the first. The
existing eager implementation must fail this test; the streaming
implementation must pass it while retaining all existing duplicate-name,
FRiP, per-peak, and rollback tests.
