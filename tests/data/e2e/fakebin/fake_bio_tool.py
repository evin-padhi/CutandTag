#!/usr/bin/env python3
"""Deterministic multi-call fake bioinformatics tools for the offline E2E profile."""

from __future__ import annotations

import gzip
import json
import shutil
import sys
from pathlib import Path


def argument_value(args: list[str], option: str, default: str | None = None) -> str:
    if option not in args:
        if default is None:
            raise SystemExit(f"fake tool: missing required option {option}")
        return default
    return args[args.index(option) + 1]


def read_text(path: str) -> str:
    if path == "-":
        return sys.stdin.read()
    return Path(path).read_text(encoding="utf-8")


def write_text(path: str, text: str) -> None:
    if path == "-":
        sys.stdout.write(text)
    else:
        Path(path).write_text(text, encoding="utf-8")


def read_fastq_names(path: str) -> list[str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = list(handle)
    return [lines[index].strip()[1:].split()[0].removesuffix("/1") for index in range(0, len(lines), 4)]


def fasta_records(path: str) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name is not None:
            records[name].append(line.strip())
    return {key: "".join(value) for key, value in records.items()}


def bed_rows(path: str) -> list[list[str]]:
    text = read_text(path)
    return [
        line.split("\t")
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def overlap(left: list[str], right: list[str]) -> bool:
    return (
        left[0] == right[0]
        and int(left[1]) < int(right[2])
        and int(right[1]) < int(left[2])
    )


def version(tool: str) -> bool:
    if any(value in sys.argv[1:] for value in ("--version", "-version")):
        versions = {
            "fastqc": "FastQC v0.12.1",
            "bowtie2": "bowtie2-align-s version 2.5.4",
            "samtools": "samtools 1.20",
            "bamCoverage": "bamCoverage 3.5.5",
            "computeMatrix": "computeMatrix 3.5.5",
            "plotProfile": "plotProfile 3.5.5",
            "macs2": "macs2 2.2.9.1",
            "bedtools": "bedtools v2.31.1",
            "ame": "5.5.7",
            "streme": "5.5.7",
            "fimo": "5.5.7",
            "multiqc": "multiqc, version 1.25.2",
        }
        print(versions[tool])
        return True
    return False


def fake_fastqc(args: list[str]) -> None:
    outdir = Path(argument_value(args, "--outdir"))
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [value for value in args if value.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))]
    for value in inputs:
        name = Path(value).name
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        (outdir / f"{name}_fastqc.html").write_text("<html>fake FastQC</html>\n")
        (outdir / f"{name}_fastqc.zip").write_text("fake FastQC archive\n")


def fake_bowtie2_build(args: list[str]) -> None:
    prefix = Path(args[-1])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("1", "2", "3", "4", "rev.1", "rev.2"):
        Path(f"{prefix}.{suffix}.bt2").write_text("fake bowtie2 index\n")


def fake_bowtie2(args: list[str]) -> None:
    r1 = argument_value(args, "-1")
    output = argument_value(args, "-S")
    names = read_fastq_names(r1)
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chrMini\tLN:4000"]
    for index, name in enumerate(names):
        first = 100 + index * 120
        second = first + 50
        lines.extend(
            [
                f"{name}\t99\tchrMini\t{first}\t60\t50M\t=\t{second}\t100\t"
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t"
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                f"{name}\t147\tchrMini\t{second}\t60\t50M\t=\t{first}\t-100\t"
                "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\t"
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ]
        )
    Path(output).write_text("\n".join(lines) + "\n")
    print(f"{len(names) * 2} reads; 100.00% overall alignment rate", file=sys.stderr)


def sam_records(text: str) -> list[list[str]]:
    return [
        line.split("\t")
        for line in text.splitlines()
        if line and not line.startswith("@")
    ]


def fake_samtools(args: list[str]) -> None:
    if not args:
        raise SystemExit("fake samtools: missing subcommand")
    command = args[0]
    rest = args[1:]
    if command == "quickcheck":
        return
    if command == "index":
        output = Path(rest[-1])
        output.write_text("fake BAM index\n")
        return
    if command in {"collate", "fixmate"}:
        source = rest[-1] if rest else "-"
        sys.stdout.write(read_text(source))
        return
    if command == "sort":
        output = argument_value(rest, "-o", "-")
        source = rest[-1] if rest and rest[-1] != output else "-"
        text = read_text(source) if source != "-" else sys.stdin.read()
        write_text(output, text)
        return
    if command == "view":
        output = argument_value(rest, "-o", "-")
        source = rest[-1] if rest and rest[-1] != output else "-"
        text = read_text(source) if source != "-" else sys.stdin.read()
        records = sam_records(text)
        required = int(argument_value(rest, "-f", "0"))
        excluded = int(argument_value(rest, "-F", "0"))
        selected = [
            record
            for record in records
            if int(record[1]) & required == required and int(record[1]) & excluded == 0
        ]
        if "-c" in rest:
            print(len(selected))
        elif "-b" in rest:
            headers = [line for line in text.splitlines() if line.startswith("@")]
            write_text(output, "\n".join(headers + ["\t".join(row) for row in selected]) + "\n")
        else:
            write_text(output, "\n".join("\t".join(row) for row in selected) + ("\n" if selected else ""))
        return
    if command == "flagstat":
        records = sam_records(read_text(rest[-1]))
        total = len(records)
        print(f"{total} + 0 in total (QC-passed reads + QC-failed reads)")
        print(f"{total} + 0 mapped (100.00% : N/A)")
        print(f"{total} + 0 properly paired (100.00% : N/A)")
        return
    if command == "stats":
        records = sam_records(read_text(rest[-1]))
        print(f"SN\traw total sequences:\t{len(records)}")
        print(f"SN\treads mapped:\t{len(records)}")
        print("IS\t100\t2\t2\t0\t0")
        return
    if command == "idxstats":
        records = sam_records(read_text(rest[-1]))
        print(f"chrMini\t4000\t{len(records)}\t0")
        print("*\t0\t0\t0")
        return
    if command == "markdup":
        metrics = Path(argument_value(rest, "-f"))
        metrics.write_text(
            json.dumps(
                {
                    "READ": 4,
                    "EXAMINED": 4,
                    "DUPLICATE TOTAL": 0,
                    "ESTIMATED LIBRARY SIZE": 2,
                }
            )
            + "\n"
        )
        return
    raise SystemExit(f"fake samtools: unsupported subcommand {command}")


def fake_bamcoverage(args: list[str]) -> None:
    Path(argument_value(args, "--outFileName")).write_text("fake bigWig\n")


def fake_compute_matrix(args: list[str]) -> None:
    matrix_path = Path(argument_value(args, "--outFileName"))
    matrix_table_path = Path(argument_value(args, "--outFileNameMatrix"))
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    matrix_table_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(matrix_path, "wt", encoding="utf-8") as handle:
        handle.write("fake deepTools matrix\n")
    matrix_table_path.write_text("fake deepTools matrix table\n", encoding="utf-8")


def fake_plot_profile(args: list[str]) -> None:
    profile_path = Path(argument_value(args, "--outFileName"))
    profile_table_path = Path(argument_value(args, "--outFileNameData"))
    profile_path.parent.mkdir(parents=True, exist_ok=True)
    profile_table_path.parent.mkdir(parents=True, exist_ok=True)
    values = ["4"] * 600
    values[:10] = ["2"] * 10
    values[-10:] = ["2"] * 10
    values[300] = "12"
    labels = [""] * 600
    labels[0] = "-3.0Kb"
    labels[299] = "TSS"
    labels[599] = "3.0Kb"
    bins = [str(index) for index in range(1, 601)]
    profile_path.write_text("fake PNG placeholder\n", encoding="utf-8")
    profile_table_path.write_text(
        "bin labels\t\t" + "\t".join(labels) + "\n"
        "bins\t\t" + "\t".join(bins) + "\n"
        "coverage.RPKM\tgenes\t" + "\t".join(values) + "\n",
        encoding="utf-8",
    )


def fake_macs2(args: list[str]) -> None:
    if not args or args[0] != "callpeak":
        raise SystemExit("fake macs2 supports only callpeak")
    name = argument_value(args, "-n")
    Path(f"{name}_peaks.xls").write_text("# fake MACS2 output\n")
    if "--broad" in args:
        peak = "chrMini\t80\t320\tpeak_1\t100\t.\t8\t12\t10\n"
        Path(f"{name}_peaks.broadPeak").write_text(peak)
        Path(f"{name}_peaks.gappedPeak").write_text(peak)
    else:
        Path(f"{name}_peaks.narrowPeak").write_text(
            "chrMini\t100\t300\tpeak_1\t100\t.\t8\t12\t10\t100\n"
        )
        Path(f"{name}_summits.bed").write_text("chrMini\t200\t201\tpeak_1\t10\n")


def fake_bedtools(args: list[str]) -> None:
    if not args:
        raise SystemExit("fake bedtools: missing subcommand")
    command = args[0]
    rest = args[1:]
    if command == "intersect":
        left = bed_rows(argument_value(rest, "-a"))
        right = bed_rows(argument_value(rest, "-b"))
        keep_overlap = "-u" in rest
        result = [
            row
            for row in left
            if any(overlap(row, other) for other in right) == keep_overlap
        ]
        write_text("-", "\n".join("\t".join(row) for row in result) + ("\n" if result else ""))
        return
    if command == "sort":
        rows = bed_rows(argument_value(rest, "-i"))
        rows.sort(key=lambda row: (row[0], int(row[1]), int(row[2])))
        write_text("-", "\n".join("\t".join(row) for row in rows) + ("\n" if rows else ""))
        return
    if command == "merge":
        rows = bed_rows(argument_value(rest, "-i"))
        write_text("-", "\n".join("\t".join(row[:3]) for row in rows) + ("\n" if rows else ""))
        return
    if command == "shuffle":
        templates = bed_rows(argument_value(rest, "-i"))
        result = []
        for index, row in enumerate(templates):
            width = int(row[2]) - int(row[1])
            start = 500 + index * (width + 10)
            result.append([row[0], str(start), str(start + width), row[3]])
        write_text("-", "\n".join("\t".join(row) for row in result) + "\n")
        return
    if command == "getfasta":
        reference = fasta_records(argument_value(rest, "-fi"))
        rows = bed_rows(argument_value(rest, "-bed"))
        tabular = "-tab" in rest
        output = []
        for row in rows:
            name = row[3] if len(row) > 3 else f"{row[0]}:{row[1]}-{row[2]}"
            sequence = reference[row[0]][int(row[1]) : int(row[2])]
            if tabular:
                output.append(f"{name}\t{sequence}")
            else:
                output.extend([f">{name}", sequence])
        write_text("-", "\n".join(output) + ("\n" if output else ""))
        return
    raise SystemExit(f"fake bedtools: unsupported subcommand {command}")


def fake_ame(args: list[str]) -> None:
    sys.stdout.write(
        "rank\tmotif_ID\tmotif_Alt_ID\tp-value\tadj_p-value\tpos\tneg\tenrichment\n"
        "1\tMA0139.1\tCTCF\t0.001\t0.01\t1\t1\t5.0\n"
    )


def fake_streme(args: list[str]) -> None:
    outdir = Path(argument_value(args, "--oc"))
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "streme.txt").write_text("fake de novo motif\n")


def fake_fimo(args: list[str]) -> None:
    outdir = Path(argument_value(args, "--oc"))
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "fimo.tsv").write_text(
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\t"
        "p-value\tq-value\tmatched_sequence\n"
        "MA0139.1\tCTCF\tpeak_000001\t96\t103\t+\t12\t0.0001\t0.001\tACGTACGT\n"
    )


def fake_multiqc(args: list[str]) -> None:
    outdir = Path(argument_value(args, "--outdir"))
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "multiqc_report.html").write_text("<html>fake MultiQC report</html>\n")
    data = outdir / "multiqc_data"
    data.mkdir(exist_ok=True)
    (data / "multiqc_data.json").write_text('{"fake": true}\n')


def main() -> None:
    tool = Path(sys.argv[0]).name
    args = sys.argv[1:]
    if args[:1] == ["--fake-tool"] and len(args) >= 2:
        tool = args[1]
        args = args[2:]
    if version(tool):
        return
    handlers = {
        "fastqc": fake_fastqc,
        "bowtie2-build": fake_bowtie2_build,
        "bowtie2": fake_bowtie2,
        "samtools": fake_samtools,
        "bamCoverage": fake_bamcoverage,
        "computeMatrix": fake_compute_matrix,
        "plotProfile": fake_plot_profile,
        "macs2": fake_macs2,
        "bedtools": fake_bedtools,
        "ame": fake_ame,
        "streme": fake_streme,
        "fimo": fake_fimo,
        "multiqc": fake_multiqc,
    }
    try:
        handler = handlers[tool]
    except KeyError as error:
        raise SystemExit(f"unsupported fake tool name: {tool}") from error
    handler(args)


if __name__ == "__main__":
    main()
