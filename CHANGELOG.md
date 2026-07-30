# Changelog

All notable changes to this project are documented in this file.

## Unreleased

### Added

- Added a self-contained consolidated QC dashboard with joined per-sample TSV
  and JSON summaries, a TSS scalar and profile export, and a top-ten AME motif
  summary. These reports are descriptive data products and do not make
  biological classifications.
- Added schema-versioned, generator-versioned dashboard JSON with explicit
  target/control counts, per-family availability, structured warnings, and
  insert-size and peak-width distributions.

### Changed

- Corrected native deepTools 3.5.5 profile parsing and made the offline fake
  byte-compatible with the pinned producer format.
- The `QC` subworkflow now accepts AME result and AME status channels as its
  tenth and eleventh inputs. The `TSS_ENRICHMENT.out.profiles` tuple now emits
  the native profile table as its sixth element and the status table as its
  seventh. Direct Nextflow importers must update their calls/destructuring.
- Dashboard TSV booleans now use lowercase machine values, optional-analysis
  states distinguish skipped/empty/missing/failed results, and the five-file
  direct-CLI output is published as a rollback-safe bundle.

## 0.1.0 — 2026-07-23

- Added strict CSV manifest validation and synchronized streaming I2
  demultiplexing.
- Added Nextflow DSL2 modules for read QC, Bowtie2 alignment, SAMtools BAM/QC
  processing, deepTools coverage/TSS profiles, matched-IgG MACS2 broad peaks,
  optional motif-only narrow peaks, blacklist filtering, fragment-based FRiP,
  MEME Suite motif analysis, and MultiQC reporting.
- Added exact pinned Conda definitions and process-specific container images.
- Added the six-library NX701–NX706 example mapping and deterministic
  synthetic test fixtures.
- Added unit, static, direct-fixture, environment, documentation, integration,
  end-to-end, and resume-aware assertions.
- Documented local/shared-HPC Conda caching, all pipeline parameters, output
  interpretation, resume behavior, troubleshooting, and runtime limitations.

### Verification limitation

Nextflow and Conda were unavailable in the implementation workspace. Runtime
DSL2 execution, environment solving, real tool execution, container pulls,
and cache/resume behavior therefore remain to be verified on a supported
system; they are not claimed as passing in this release.
