#!/usr/bin/env Rscript
# Generate reference values from R/rrBLUP::mixed.solve for Rust validation
# Run with: Rscript generate_reference_values.R

library(rrBLUP)

cat("=" , rep("=", 70), "\n", sep="")
cat("Reference values from R/rrBLUP::mixed.solve v", as.character(packageVersion("rrBLUP")), "\n", sep="")
cat("=" , rep("=", 70), "\n\n", sep="")

# Helper to print results in Rust-friendly format
print_result <- function(name, result, precision = 10) {
  cat(sprintf("// Test: %s\n", name))
  cat(sprintf("// Vu = %.10f\n", result$Vu))
  cat(sprintf("// Ve = %.10f\n", result$Ve))
  cat(sprintf("// LL = %.10f\n", result$LL))
  cat("// beta = [", paste(sprintf("%.10f", result$beta), collapse = ", "), "]\n", sep="")
  cat("// u = [", paste(sprintf("%.10f", result$u), collapse = ", "), "]\n", sep="")
  if (!is.null(result$beta.SE)) {
    cat("// beta_se = [", paste(sprintf("%.10f", result$beta.SE), collapse = ", "), "]\n", sep="")
  }
  if (!is.null(result$u.SE)) {
    cat("// u_se = [", paste(sprintf("%.10f", result$u.SE), collapse = ", "), "]\n", sep="")
  }
  cat("\n")
}

# ============================================================
# Test 1: Simple intercept-only model (y only)
# ============================================================
cat("Test 1: Simple intercept-only model\n")
cat("-" , rep("-", 50), "\n", sep="")
set.seed(42)
y1 <- c(1.0, 2.0, 3.0, 4.0, 5.0)
result1 <- mixed.solve(y1)
print_result("test1_simple_intercept", result1)
cat("Input y = [", paste(y1, collapse = ", "), "]\n\n", sep="")

# ============================================================
# Test 2: Larger dataset with intercept only
# ============================================================
cat("Test 2: Larger dataset (n=10)\n")
cat("-" , rep("-", 50), "\n", sep="")
y2 <- c(1.2, 2.5, 3.1, 4.8, 5.2, 6.0, 7.3, 8.1, 9.5, 10.2)
result2 <- mixed.solve(y2)
print_result("test2_larger_intercept", result2)
cat("Input y = [", paste(y2, collapse = ", "), "]\n\n", sep="")

# ============================================================
# Test 3: With fixed effects matrix X (intercept + covariate)
# ============================================================
cat("Test 3: With fixed effects (intercept + covariate)\n")
cat("-" , rep("-", 50), "\n", sep="")
y3 <- c(1.0, 2.5, 4.2, 5.8, 7.1, 8.9)
X3 <- cbind(rep(1, 6), c(0, 1, 2, 3, 4, 5))  # Intercept + linear covariate
result3 <- mixed.solve(y3, X = X3)
print_result("test3_with_fixed_effects", result3)
cat("Input y = [", paste(y3, collapse = ", "), "]\n", sep="")
cat("Input X (6x2):\n")
print(X3)
cat("\n")

# ============================================================
# Test 4: With kinship matrix K
# ============================================================
cat("Test 4: With kinship matrix K\n")
cat("-" , rep("-", 50), "\n", sep="")
y4 <- c(2.1, 3.5, 4.2, 5.8, 6.1)
# Simple positive definite K matrix
K4 <- matrix(c(
  1.0, 0.5, 0.3, 0.2, 0.1,
  0.5, 1.0, 0.4, 0.3, 0.2,
  0.3, 0.4, 1.0, 0.5, 0.3,
  0.2, 0.3, 0.5, 1.0, 0.4,
  0.1, 0.2, 0.3, 0.4, 1.0
), nrow = 5, byrow = TRUE)
result4 <- mixed.solve(y4, K = K4)
print_result("test4_with_kinship", result4)
cat("Input y = [", paste(y4, collapse = ", "), "]\n", sep="")
cat("Input K (5x5):\n")
print(K4)
cat("\n")

# ============================================================
# Test 5: With Z matrix (random effects design)
# ============================================================
cat("Test 5: With Z matrix (n=6, m=3)\n")
cat("-" , rep("-", 50), "\n", sep="")
y5 <- c(1.5, 2.3, 3.1, 4.2, 5.0, 5.8)
# Z matrix mapping 6 observations to 3 random effects
Z5 <- matrix(c(
  1, 0, 0,
  1, 0, 0,
  0, 1, 0,
  0, 1, 0,
  0, 0, 1,
  0, 0, 1
), nrow = 6, byrow = TRUE)
result5 <- mixed.solve(y5, Z = Z5)
print_result("test5_with_z_matrix", result5)
cat("Input y = [", paste(y5, collapse = ", "), "]\n", sep="")
cat("Input Z (6x3):\n")
print(Z5)
cat("\n")

# ============================================================
# Test 6: With SE=TRUE (standard errors)
# ============================================================
cat("Test 6: With standard errors (SE=TRUE)\n")
cat("-" , rep("-", 50), "\n", sep="")
y6 <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
result6 <- mixed.solve(y6, SE = TRUE)
print_result("test6_with_se", result6)
cat("Input y = [", paste(y6, collapse = ", "), "]\n\n", sep="")

# ============================================================
# Test 7: ML method instead of REML
# ============================================================
cat("Test 7: ML method (instead of REML)\n")
cat("-" , rep("-", 50), "\n", sep="")
y7 <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
result7_reml <- mixed.solve(y7, method = "REML")
result7_ml <- mixed.solve(y7, method = "ML")
cat("REML results:\n")
print_result("test7_reml", result7_reml)
cat("ML results:\n")
print_result("test7_ml", result7_ml)
cat("Input y = [", paste(y7, collapse = ", "), "]\n\n", sep="")

# ============================================================
# Test 8: Complete model with X, Z, K, SE
# ============================================================
cat("Test 8: Complete model (X, Z, K, SE=TRUE)\n")
cat("-" , rep("-", 50), "\n", sep="")
y8 <- c(2.1, 3.2, 4.5, 5.1, 6.3, 7.2)
X8 <- cbind(rep(1, 6), c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0))
Z8 <- matrix(c(
  1, 0, 0,
  1, 0, 0,
  0, 1, 0,
  0, 1, 0,
  0, 0, 1,
  0, 0, 1
), nrow = 6, byrow = TRUE)
K8 <- matrix(c(
  1.0, 0.3, 0.1,
  0.3, 1.0, 0.2,
  0.1, 0.2, 1.0
), nrow = 3, byrow = TRUE)
result8 <- mixed.solve(y8, Z = Z8, K = K8, X = X8, SE = TRUE)
print_result("test8_complete", result8)
cat("Input y = [", paste(y8, collapse = ", "), "]\n", sep="")
cat("Input X (6x2):\n")
print(X8)
cat("Input Z (6x3):\n")
print(Z8)
cat("Input K (3x3):\n")
print(K8)
cat("\n")

# ============================================================
# Test 9: With NA values in y
# ============================================================
cat("Test 9: With NA values in y\n")
cat("-" , rep("-", 50), "\n", sep="")
y9 <- c(1.0, NA, 3.0, NA, 5.0, 6.0, 7.0, NA, 9.0, 10.0)
result9 <- mixed.solve(y9)
print_result("test9_with_na", result9)
cat("Input y = [", paste(ifelse(is.na(y9), "NA", y9), collapse = ", "), "]\n", sep="")
cat("Non-NA indices: ", which(!is.na(y9)), "\n\n", sep="")

# ============================================================
# Test 10: Genomic prediction scenario (larger, realistic)
# ============================================================
cat("Test 10: Genomic prediction scenario (n=20, m=20)\n")
cat("-" , rep("-", 50), "\n", sep="")
set.seed(123)
n10 <- 20
# Simulate a simple kinship matrix
G10 <- matrix(rnorm(n10 * 10), nrow = n10)  # 10 markers
K10 <- tcrossprod(G10) / 10  # Kinship from markers
# Make it positive definite
K10 <- K10 + diag(n10) * 0.1
# Simulate phenotypes
true_u <- rnorm(n10)
y10 <- 5 + K10 %*% true_u * 0.5 + rnorm(n10) * 0.3
y10 <- as.vector(y10)

result10 <- mixed.solve(y10, K = K10, SE = TRUE)
print_result("test10_genomic", result10)
cat("Input y (first 10): [", paste(sprintf("%.4f", y10[1:10]), collapse = ", "), ", ...]\n", sep="")
cat("K diagonal (first 5): [", paste(sprintf("%.4f", diag(K10)[1:5]), collapse = ", "), ", ...]\n\n", sep="")

cat("=" , rep("=", 70), "\n", sep="")
cat("Reference generation complete.\n")
cat("Copy the values above into Rust tests for validation.\n")
