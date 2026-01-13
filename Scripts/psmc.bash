#!/bin/bash
# Script Name: psmc.bash
# Description: Run PSMC analysis for a species
# Author: Aure Kylmänen
# ====================================
set -euo pipefail

# ============================
# USAGE
# ============================
# ./psmc.bash <species> <params.txt>

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <species> <params.txt>"
    exit 1
fi

SPECIES="$1"
PARAMS_FILE="$2"

# ============================
# LOAD PARAMETERS
# ============================
if [[ ! -f "$PARAMS_FILE" ]]; then
    echo "ERROR: Parameter file not found: $PARAMS_FILE"
    exit 1
fi

source "$PARAMS_FILE"

# ============================
# REQUIRED PARAMETERS CHECK
# ============================
: "${WORKING_DIR:?Missing WORKING_DIR in params file}"
: "${PSMC:?Missing PSMC executable in params file}"
: "${PSMC_UTILS:?Missing PSMC_UTILS path in params file}"
: "${DIPLOID_FASTQ:?Missing DIPLOID_FASTQ in params file}"
: "${PSMC_N_ITER:?Missing PSMC_N_ITER in params file}"
: "${PSMC_T:?Missing PSMC_T in params file}"
: "${PSMC_R:?Missing PSMC_R in params file}"
: "${PSMC_PATTERN:?Missing PSMC_PATTERN in params file}"

# ============================
# DIRECTORIES
# ============================
RESULTS_DIR="${WORKING_DIR}/${SPECIES}/results"
mkdir -p "$RESULTS_DIR"

cd "${WORKING_DIR}/${SPECIES}" || { echo "ERROR: Failed to cd into ${WORKING_DIR}/${SPECIES}"; exit 1; }

# ============================
# CONVERT FASTQ TO PSMCFA
# ============================

echo "Converting FASTQ to PSMCFA..."
${PSMC_UTILS}/fq2psmcfa -q20 "${DIPLOID_FASTQ}" > "${DIPLOID_PSMCFA}"

# ============================
# RUN PSMC
# ============================

echo "Running PSMC with parameters:"
echo "Iterations: $PSMC_N_ITER, t: $PSMC_T, r: $PSMC_R, pattern: $PSMC_PATTERN"

"${PSMC}" \
    -N "$PSMC_N_ITER" \
    -t "$PSMC_T" \
    -r "$PSMC_R" \
    -p "$PSMC_PATTERN" \
    -o "${DIPLOID_PSMC}" \
    "${DIPLOID_PSMCFA}"


# ============================
# GENERATE HISTORICAL INFERENCE
# ============================
echo "Generating historical inference scripts..."
"${PSMC_UTILS}/psmc2history.pl" "${DIPLOID_PSMC}" | "${PSMC_UTILS}/history2ms.pl" > "${SPECIES}_ms-cmd.sh"

"${PSMC_UTILS}/psmc_plot.pl" -u "${MUTATION_RATE}" -g "${GENERATION_TIME}" -p "${WORKING_DIR}/${SPECIES}/diploid" "${DIPLOID_PSMC}"


# ============================
# COPY RESULTS
# ============================
cp diploid.psmc "${RESULTS_DIR}/" 2>/dev/null || echo "⚠️ Failed to copy diploid.psmc"
cp "${SPECIES}_ms-cmd.sh" "${RESULTS_DIR}/" 2>/dev/null || echo "⚠️ Failed to copy ms command"
cp diploid*.png "${RESULTS_DIR}/" 2>/dev/null || echo "⚠️ Failed to copy PSMC plots"

echo "PSMC analysis completed successfully."
