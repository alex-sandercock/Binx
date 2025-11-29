#!/usr/bin/env Rscript
# Generate reference values from R/rrBLUP::A.mat for Rust validation
# Run with: Rscript generate_amat_reference.R

library(rrBLUP)

cat("=" , rep("=", 70), "\n", sep="")
cat("Reference values from R/rrBLUP::A.mat v", as.character(packageVersion("rrBLUP")), "\n", sep="")
cat("=" , rep("=", 70), "\n\n", sep="")

# Helper to print matrix in Rust-friendly format
print_matrix <- function(name, mat, precision = 10) {
  cat(sprintf("// %s (%d x %d):\n", name, nrow(mat), ncol(mat)))
  cat("// [")
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      cat(sprintf("%.10f", mat[i,j]))
      if (i != nrow(mat) || j != ncol(mat)) cat(", ")
    }
    if (i != nrow(mat)) cat("\n//  ")
  }
  cat("]\n")
}

# ============================================================
# Test 1: Simple case (no missing, no filtering)
# ============================================================
cat("Test 1: Simple case (3 individuals x 5 markers)\n")
cat("-" , rep("-", 50), "\n", sep="")

set.seed(42)
X1 <- matrix(c(
  -1,  0,  1,  0, -1,
   0,  1, -1,  1,  0,
   1, -1,  0, -1,  1
), nrow = 3, byrow = TRUE)

cat("Input X:\n")
print(X1)

A1 <- A.mat(X1)
cat("\nOutput A:\n")
print(A1, digits = 10)
print_matrix("A", A1)
cat("\n")

# ============================================================
# Test 2: Larger case (5 individuals x 10 markers)
# ============================================================
cat("Test 2: Larger case (5 individuals x 10 markers)\n")
cat("-" , rep("-", 50), "\n", sep="")

set.seed(123)
X2 <- matrix(sample(c(-1, 0, 1), 50, replace = TRUE), nrow = 5, ncol = 10)

cat("Input X:\n")
print(X2)

A2 <- A.mat(X2)
cat("\nOutput A:\n")
print(A2, digits = 10)
print_matrix("A", A2)
cat("\n")

# ============================================================
# Test 3: With MAF filtering
# ============================================================
cat("Test 3: With MAF filtering (min.MAF = 0.2)\n")
cat("-" , rep("-", 50), "\n", sep="")

# Include a low-MAF marker
X3 <- matrix(c(
  -1, -1, -1,  0,  1,  # marker 1: low MAF
   0,  1, -1,  1,  0,
   1, -1,  0, -1,  1,
   0,  0,  1,  0, -1
), nrow = 4, byrow = TRUE)

cat("Input X:\n")
print(X3)

A3 <- A.mat(X3, min.MAF = 0.2)
cat("\nOutput A (min.MAF = 0.2):\n")
print(A3, digits = 10)
print_matrix("A", A3)
cat("\n")

# ============================================================
# Test 4: With missing values (mean imputation)
# ============================================================
cat("Test 4: With missing values (mean imputation)\n")
cat("-" , rep("-", 50), "\n", sep="")

X4 <- matrix(c(
  -1,  0,  1,  0,
   0, NA, -1,  1,
   1, -1, NA, -1,
   0,  1,  0, NA
), nrow = 4, byrow = TRUE)

cat("Input X (with NA):\n")
print(X4)

A4 <- A.mat(X4, impute.method = "mean")
cat("\nOutput A (mean imputation):\n")
print(A4, digits = 10)
print_matrix("A", A4)
cat("\n")

# ============================================================
# Test 5: With shrinkage (EJ method)
# ============================================================
cat("Test 5: With shrinkage (EJ method)\n")
cat("-" , rep("-", 50), "\n", sep="")

set.seed(456)
X5 <- matrix(sample(c(-1, 0, 1), 80, replace = TRUE), nrow = 8, ncol = 10)

cat("Input X (8 x 10):\n")
print(X5)

A5 <- A.mat(X5, shrink = TRUE)
cat("\nOutput A (EJ shrinkage):\n")
print(A5, digits = 10)
print_matrix("A", A5)
cat("\n")

# ============================================================
# Test 6: Return imputed matrix
# ============================================================
cat("Test 6: Return imputed matrix\n")
cat("-" , rep("-", 50), "\n", sep="")

X6 <- matrix(c(
  -1,  0,  1,
   0, NA, -1,
   1, -1,  0
), nrow = 3, byrow = TRUE)

cat("Input X (with NA):\n")
print(X6)

result6 <- A.mat(X6, impute.method = "mean", return.imputed = TRUE)
cat("\nOutput A:\n")
print(result6$A, digits = 10)
print_matrix("A", result6$A)
cat("\nImputed X:\n")
print(result6$imputed, digits = 10)
print_matrix("imputed", result6$imputed)
cat("\n")

# ============================================================
# Test 7: Realistic genomic data scenario
# ============================================================
cat("Test 7: Realistic scenario (20 individuals x 100 markers)\n")
cat("-" , rep("-", 50), "\n", sep="")

set.seed(789)
n7 <- 20
m7 <- 100

# Generate realistic genotype frequencies
freqs <- runif(m7, 0.1, 0.9)
X7 <- matrix(0, nrow = n7, ncol = m7)
for (j in 1:m7) {
  p <- freqs[j]
  probs <- c((1-p)^2, 2*p*(1-p), p^2)  # HWE
  X7[, j] <- sample(c(-1, 0, 1), n7, replace = TRUE, prob = probs)
}

A7 <- A.mat(X7)
cat("Output A diagonal (first 5):\n")
print(diag(A7)[1:5], digits = 10)
cat("\nOutput A off-diagonal (A[1,2:5]):\n")
print(A7[1, 2:5], digits = 10)

# Summary stats
cat("\nA matrix summary:\n")
cat(sprintf("  mean(diag(A)) = %.10f\n", mean(diag(A7))))
cat(sprintf("  range(A) = [%.10f, %.10f]\n", min(A7), max(A7)))
print_matrix("A[1:5, 1:5]", A7[1:5, 1:5])
cat("\n")

# ============================================================
# Test 8: All heterozygotes (edge case)
# ============================================================
cat("Test 8: All heterozygotes (edge case)\n")
cat("-" , rep("-", 50), "\n", sep="")

X8 <- matrix(0, nrow = 3, ncol = 4)  # All zeros = all heterozygotes

cat("Input X (all zeros):\n")
print(X8)

A8 <- A.mat(X8)
cat("\nOutput A:\n")
print(A8, digits = 10)
print_matrix("A", A8)
cat("\n")

# ============================================================
# Verification: Compare to VanRaden formula
# ============================================================
cat("Verification: Manual VanRaden calculation for Test 1\n")
cat("-" , rep("-", 50), "\n", sep="")

# VanRaden (2008) formula: A = WW'/(2*sum(p*(1-p)))
# where W = X - 2p (centered)
# But rrBLUP uses: W = X + 1 - 2*freq, A = WW'/(2*mean(p*(1-p))*m)

# For X1:
freq1 <- apply(X1 + 1, 2, mean) / 2
cat("Allele frequencies:\n")
print(freq1)

# Center: W = X + 1 - 2*freq
W1 <- X1 + 1 - 2*matrix(freq1, nrow = 3, ncol = 5, byrow = TRUE)
cat("\nCentered W:\n")
print(W1)

# var.A = 2 * mean(p*(1-p))
var.A <- 2 * mean(freq1 * (1 - freq1))
cat(sprintf("\nvar.A = 2 * mean(p*(1-p)) = %.10f\n", var.A))

# A = WW' / (var.A * m)
A_manual <- tcrossprod(W1) / (var.A * 5)
cat("\nManual A = WW'/(var.A * m):\n")
print(A_manual, digits = 10)

cat("\nDifference from A.mat:\n")
print(A1 - A_manual, digits = 10)

cat("\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("Reference generation complete.\n")
