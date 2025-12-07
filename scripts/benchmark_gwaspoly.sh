#!/bin/bash
# Benchmark: R/GWASpoly vs gwaspoly-rs (sequential) vs gwaspoly-rs (parallel)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

# Build release binary
echo "Building binx in release mode..."
cargo build --release -p binx-cli 2>/dev/null

BINX="$PROJECT_ROOT/target/release/binx"

# Test data paths
POTATO_GENO="tests/parity/data/potato/new_potato_geno.csv"
POTATO_PHENO="tests/parity/data/potato/new_potato_pheno.csv"
TOY_GENO="tests/parity/data/toy/toy.geno.tsv"
TOY_PHENO="tests/parity/data/toy/toy.pheno.tsv"

echo ""
echo "=========================================="
echo "GWASpoly Benchmark: R vs Rust"
echo "=========================================="

# Check if R and GWASpoly are available
if command -v Rscript &> /dev/null && Rscript -e "library(GWASpoly)" 2>/dev/null; then
    HAS_R=1
    echo "R/GWASpoly: Available"
else
    HAS_R=0
    echo "R/GWASpoly: Not available (skipping R benchmarks)"
fi

echo ""
echo "Test 1: Toy dataset (small)"
echo "---------------------------"
if [ -f "$TOY_GENO" ] && [ -f "$TOY_PHENO" ]; then
    # Rust sequential
    echo -n "  gwaspoly-rs (sequential): "
    START=$(python3 -c 'import time; print(time.time())')
    $BINX gwaspoly \
        --geno "$TOY_GENO" \
        --pheno "$TOY_PHENO" \
        --trait-name Trait1 \
        --ploidy 4 \
        --models additive,general \
        --out /tmp/bench_toy_rust.csv 2>/dev/null
    END=$(python3 -c 'import time; print(time.time())')
    TIME=$(python3 -c "print(f'{($END - $START)*1000:.2f} ms')")
    echo "$TIME"
    MARKERS=$(wc -l < /tmp/bench_toy_rust.csv)
    echo "    Results: $((MARKERS - 1)) marker-model combinations"

    # R GWASpoly
    if [ $HAS_R -eq 1 ]; then
        echo -n "  R/GWASpoly: "
        START=$(python3 -c 'import time; print(time.time())')
        Rscript scripts/gwaspoly/run_parity.R \
            --geno="$TOY_GENO" \
            --pheno="$TOY_PHENO" \
            --trait=Trait1 \
            --ploidy=4 \
            --models=additive,general \
            --out=/tmp/bench_toy_r.tsv 2>/dev/null
        END=$(python3 -c 'import time; print(time.time())')
        TIME=$(python3 -c "print(f'{($END - $START)*1000:.2f} ms')")
        echo "$TIME"
    fi
else
    echo "  Toy dataset not found"
fi

echo ""
echo "Test 2: Potato dataset (medium)"
echo "--------------------------------"
if [ -f "$POTATO_GENO" ] && [ -f "$POTATO_PHENO" ]; then
    # Rust sequential
    echo -n "  gwaspoly-rs (sequential): "
    START=$(python3 -c 'import time; print(time.time())')
    $BINX gwaspoly \
        --geno "$POTATO_GENO" \
        --pheno "$POTATO_PHENO" \
        --trait-name vine.maturity \
        --ploidy 4 \
        --models additive \
        --out /tmp/bench_potato_rust.csv 2>/dev/null
    END=$(python3 -c 'import time; print(time.time())')
    TIME=$(python3 -c "print(f'{($END - $START)*1000:.2f} ms')")
    echo "$TIME"
    MARKERS=$(wc -l < /tmp/bench_potato_rust.csv)
    echo "    Results: $((MARKERS - 1)) marker results"

    # R GWASpoly
    if [ $HAS_R -eq 1 ]; then
        echo -n "  R/GWASpoly: "
        START=$(python3 -c 'import time; print(time.time())')
        Rscript scripts/gwaspoly/run_parity.R \
            --geno="$POTATO_GENO" \
            --pheno="$POTATO_PHENO" \
            --trait=vine.maturity \
            --ploidy=4 \
            --models=additive \
            --out=/tmp/bench_potato_r.tsv 2>/dev/null
        END=$(python3 -c 'import time; print(time.time())')
        TIME=$(python3 -c "print(f'{($END - $START)*1000:.2f} ms')")
        echo "$TIME"
    fi
else
    echo "  Potato dataset not found"
fi

echo ""
echo "=========================================="
echo "Benchmark complete"
echo "=========================================="
