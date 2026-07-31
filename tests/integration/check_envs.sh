#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"

failures=0

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  failures=$((failures + 1))
}

check_environment() {
  local file=$1
  local package=$2
  local version=$3
  shift 3
  local expected_dependencies=("${package}=${version}" "$@")
  local path="envs/$file"

  if [[ ! -f "$path" ]]; then
    fail "missing $path"
    return
  fi

  mapfile -t channels < <(awk '
    /^channels:[[:space:]]*$/ { in_channels=1; next }
    in_channels && /^dependencies:[[:space:]]*$/ { exit }
    in_channels && /^[[:space:]]*-[[:space:]]+/ {
      sub(/^[[:space:]]*-[[:space:]]*/, "")
      print
    }
  ' "$path")
  if [[ ${#channels[@]} -lt 2 || ${channels[0]} != conda-forge || ${channels[1]} != bioconda ]]; then
    fail "$path must start channels with conda-forge then bioconda"
  fi

  if ! rg -qx -- "[[:space:]]*-[[:space:]]*${package}=${version}[[:space:]]*" "$path"; then
    fail "$path must pin ${package}=${version}"
  fi

  mapfile -t dependencies < <(awk '
    /^dependencies:[[:space:]]*$/ { in_dependencies=1; next }
    in_dependencies && /^[[:space:]]*-[[:space:]]+/ {
      sub(/^[[:space:]]*-[[:space:]]*/, "")
      print
    }
  ' "$path")
  if [[ "${dependencies[*]}" != "${expected_dependencies[*]}" ]]; then
    fail "$path must contain exactly these direct dependencies: ${expected_dependencies[*]}"
  fi
}

check_environment python.yml python 3.12.3 matplotlib=3.9.2
check_environment fastqc.yml fastqc 0.12.1
check_environment bowtie2.yml bowtie2 2.5.4
check_environment samtools.yml samtools 1.20
check_environment deeptools.yml deeptools 3.5.5
check_environment macs2.yml macs2 2.2.9.1
check_environment bedtools.yml bedtools 2.31.1
check_environment meme.yml meme 5.5.7
check_environment multiqc.yml multiqc 1.25.2

check_schema_contract() {
  local schema_errors
  if ! schema_errors=$(python3 -c '
import json
import sys

schema = json.load(open("nextflow_schema.json"))
properties = schema["properties"]
errors = []

for name in ("fasta", "bowtie2_index", "blacklist", "gtf", "tss_bed", "motif_db"):
    value = properties[name]
    if value.get("type") != "string":
        errors.append(f"{name} must have type string")
    if "default" in value and value["default"] is None:
        errors.append(f"{name} must omit a null default")

macs_genome_size = properties["macs_genome_size"]
variants = macs_genome_size.get("anyOf") or macs_genome_size.get("oneOf") or []
has_nonempty_string = any(
    variant.get("type") == "string" and variant.get("minLength", 0) >= 1
    for variant in variants
)
has_positive_integer = any(
    variant.get("type") == "integer" and variant.get("minimum", 0) >= 1
    for variant in variants
)
if not has_nonempty_string or not has_positive_integer:
    errors.append("macs_genome_size must accept a non-empty string or positive integer")

if errors:
    print("; ".join(errors), file=sys.stderr)
    sys.exit(1)
' 2>&1); then
    fail "$schema_errors"
  fi
}

check_project_qualified_conda_paths() {
  if ! rg -Fq 'conda "${projectDir}/envs/<tool>.yml"' nextflow.config; then
    fail 'nextflow.config must document conda "${projectDir}/envs/<tool>.yml" module paths'
  fi

  if [[ -d modules ]] && rg -n --glob '*.nf' "conda[[:space:]]+['\"](?:\./)?envs/" modules; then
    fail 'module conda directives must use ${projectDir}/envs/<tool>.yml, not relative env paths'
  fi
}

check_schema_contract
check_project_qualified_conda_paths

if (( failures > 0 )); then
  exit 1
fi

printf 'Environment definitions have the required channels and exact direct dependency pins.\n'
