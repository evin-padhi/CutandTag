# Motif BED Column Normalization Design

## Problem

`PREPARE_MOTIF_SEQUENCES` concatenates final MACS2 peaks with an optional
blacklist before sorting and merging the excluded regions. MACS2 peak files
contain nine columns, while the supplied blacklist contains three. BEDTools
requires a consistent number of fields within one input stream and exits when
it reaches the first row with a different width.

## Design

Normalize final peaks and the optional blacklist independently to BED3 before
concatenating them. The downstream `bedtools sort` and `bedtools merge`
operations use only chromosome, start, and end, so discarding annotation
columns does not change the excluded genomic intervals.

Keep the original staged peak and blacklist files unchanged because other
steps use their additional fields. Add a regression assertion that requires
both exclusion sources to be normalized before they are combined.

## Validation

The motif integration test must fail against the existing implementation,
pass after BED3 normalization is added, and retain all existing motif workflow
assertions. The complete integration suite and Nextflow lint should also pass
when their runtimes are available.
