#!/usr/bin/env Rscript
# Generate reference values from R/rrBLUP::kin.blup for Rust validation
# Run with: Rscript generate_kinblup_reference.R

library(rrBLUP)

cat("=" , rep("=", 70), "\n", sep="")
cat("Reference values from R/rrBLUP::kin.blup v", as.character(packageVersion("rrBLUP")), "\n", sep="")
cat("=" , rep("=", 70), "\n\n", sep="")

# ============================================================
# Test 1: Simple case (no kinship, replicated design)
# ============================================================
cat("Test 1: No kinship matrix (genotype as random factor)\n")
cat("-" , rep("-", 50), "\n", sep="")

set.seed(42)
data1 <- data.frame(
  gid = c("A", "B", "C", "A", "B", "C"),
  y = c(1.0, 2.0, 3.0, 1.5, 2.5, 3.5)
)

cat("Input data:\n")
print(data1)

result1 <- kin.blup(data1, geno = "gid", pheno = "y")
cat("\nResults:\n")
cat(sprintf("Vg = %.10f\n", result1$Vg))
cat(sprintf("Ve = %.10f\n", result1$Ve))
cat("g (genetic values):\n")
print(result1$g)
cat("pred (predictions):\n")
print(result1$pred)
cat("resid:\n")
print(result1$resid)
cat("\n")

# ============================================================
# Test 2: With kinship matrix
# ============================================================
cat("Test 2: With kinship matrix\n")
cat("-" , rep("-", 50), "\n", sep="")

data2 <- data.frame(
  gid = c("A", "B", "C", "D"),
  y = c(1.0, 2.5, 3.2, 4.1)
)

K2 <- matrix(c(
  1.0, 0.5, 0.3, 0.2,
  0.5, 1.0, 0.4, 0.3,
  0.3, 0.4, 1.0, 0.5,
  0.2, 0.3, 0.5, 1.0
), nrow = 4, byrow = TRUE)
rownames(K2) <- colnames(K2) <- c("A", "B", "C", "D")

cat("Input data:\n")
print(data2)
cat("\nKinship matrix K:\n")
print(K2)

result2 <- kin.blup(data2, geno = "gid", pheno = "y", K = K2)
cat("\nResults:\n")
cat(sprintf("Vg = %.10f\n", result2$Vg))
cat(sprintf("Ve = %.10f\n", result2$Ve))
cat("g (genetic values):\n")
print(result2$g)
cat("pred (predictions):\n")
print(result2$pred)
cat("\n")

# ============================================================
# Test 3: With missing phenotypes (prediction)
# ============================================================
cat("Test 3: With missing phenotypes (genomic prediction)\n")
cat("-" , rep("-", 50), "\n", sep="")

data3 <- data.frame(
  gid = c("A", "B", "C", "D", "E"),
  y = c(1.0, 2.5, NA, 4.0, NA)
)

K3 <- matrix(c(
  1.0, 0.5, 0.3, 0.2, 0.1,
  0.5, 1.0, 0.4, 0.3, 0.2,
  0.3, 0.4, 1.0, 0.5, 0.3,
  0.2, 0.3, 0.5, 1.0, 0.4,
  0.1, 0.2, 0.3, 0.4, 1.0
), nrow = 5, byrow = TRUE)
rownames(K3) <- colnames(K3) <- c("A", "B", "C", "D", "E")

cat("Input data (C and E have missing phenotypes):\n")
print(data3)

result3 <- kin.blup(data3, geno = "gid", pheno = "y", K = K3)
cat("\nResults:\n")
cat(sprintf("Vg = %.10f\n", result3$Vg))
cat(sprintf("Ve = %.10f\n", result3$Ve))
cat("g (genetic values) - includes predictions for C and E:\n")
print(result3$g)
cat("pred:\n")
print(result3$pred)
cat("\n")

# ============================================================
# Test 4: With PEV (prediction error variance)
# ============================================================
cat("Test 4: With PEV=TRUE\n")
cat("-" , rep("-", 50), "\n", sep="")

result4 <- kin.blup(data2, geno = "gid", pheno = "y", K = K2, PEV = TRUE)
cat("Results with PEV:\n")
cat(sprintf("Vg = %.10f\n", result4$Vg))
cat(sprintf("Ve = %.10f\n", result4$Ve))
cat("g:\n")
print(result4$g)
cat("PEV:\n")
print(result4$PEV)
cat("\n")

# ============================================================
# Test 5: With fixed effects
# ============================================================
cat("Test 5: With fixed effects (categorical factor)\n")
cat("-" , rep("-", 50), "\n", sep="")

data5 <- data.frame(
  gid = c("A", "B", "C", "A", "B", "C"),
  env = c("E1", "E1", "E1", "E2", "E2", "E2"),
  y = c(1.0, 2.0, 3.0, 2.0, 3.0, 4.0)
)

cat("Input data with environment fixed effect:\n")
print(data5)

result5 <- kin.blup(data5, geno = "gid", pheno = "y", fixed = "env")
cat("\nResults:\n")
cat(sprintf("Vg = %.10f\n", result5$Vg))
cat(sprintf("Ve = %.10f\n", result5$Ve))
cat("g:\n")
print(result5$g)
cat("\n")

# ============================================================
# Test 6: With covariates
# ============================================================
cat("Test 6: With covariates (continuous variable)\n")
cat("-" , rep("-", 50), "\n", sep="")

data6 <- data.frame(
  gid = c("A", "B", "C", "D"),
  cov = c(0.0, 1.0, 2.0, 3.0),
  y = c(1.0, 2.2, 3.5, 4.8)
)

K6 <- matrix(c(
  1.0, 0.5, 0.3, 0.2,
  0.5, 1.0, 0.4, 0.3,
  0.3, 0.4, 1.0, 0.5,
  0.2, 0.3, 0.5, 1.0
), nrow = 4, byrow = TRUE)
rownames(K6) <- colnames(K6) <- c("A", "B", "C", "D")

cat("Input data with covariate:\n")
print(data6)

result6 <- kin.blup(data6, geno = "gid", pheno = "y", K = K6, covariate = "cov")
cat("\nResults:\n")
cat(sprintf("Vg = %.10f\n", result6$Vg))
cat(sprintf("Ve = %.10f\n", result6$Ve))
cat("g:\n")
print(result6$g)
cat("\n")

# ============================================================
# Test 7: Larger realistic scenario
# ============================================================
cat("Test 7: Larger realistic scenario (20 genotypes, 50 observations)\n")
cat("-" , rep("-", 50), "\n", sep="")

set.seed(123)
n_geno <- 20
n_obs <- 50

# Generate random kinship matrix
G <- matrix(rnorm(n_geno * 30), nrow = n_geno)
K7 <- tcrossprod(G) / 30
K7 <- K7 + diag(n_geno) * 0.1  # Make positive definite
geno_names <- paste0("G", 1:n_geno)
rownames(K7) <- colnames(K7) <- geno_names

# Generate phenotypes
geno_effects <- rnorm(n_geno)
obs_geno <- sample(geno_names, n_obs, replace = TRUE)
y7 <- 5 + geno_effects[match(obs_geno, geno_names)] * 0.5 + rnorm(n_obs) * 0.3

data7 <- data.frame(gid = obs_geno, y = y7)

result7 <- kin.blup(data7, geno = "gid", pheno = "y", K = K7, PEV = TRUE)
cat("Results:\n")
cat(sprintf("Vg = %.10f\n", result7$Vg))
cat(sprintf("Ve = %.10f\n", result7$Ve))
cat("\nGenetic values (first 5):\n")
print(head(result7$g, 5))
cat("\nPEV (first 5):\n")
print(head(result7$PEV, 5))
cat("\n")

# ============================================================
# Test 8: With GAUSS=TRUE (Gaussian kernel)
# ============================================================
cat("Test 8: With Gaussian kernel (GAUSS=TRUE)\n")
cat("-" , rep("-", 50), "\n", sep="")

# Use a distance-like K matrix for Gaussian kernel
K8 <- matrix(c(
  0.0, 1.0, 2.0, 3.0,
  1.0, 0.0, 1.0, 2.0,
  2.0, 1.0, 0.0, 1.0,
  3.0, 2.0, 1.0, 0.0
), nrow = 4, byrow = TRUE)
rownames(K8) <- colnames(K8) <- c("A", "B", "C", "D")

data8 <- data.frame(
  gid = c("A", "B", "C", "D"),
  y = c(1.0, 1.8, 2.5, 3.5)
)

cat("Input data:\n")
print(data8)
cat("\nDistance matrix K:\n")
print(K8)

result8 <- kin.blup(data8, geno = "gid", pheno = "y", K = K8, GAUSS = TRUE)
cat("\nResults:\n")
cat(sprintf("Vg = %.10f\n", result8$Vg))
cat(sprintf("Ve = %.10f\n", result8$Ve))
cat("g:\n")
print(result8$g)
cat("Profile likelihood:\n")
print(result8$profile)
cat("\n")

cat("=" , rep("=", 70), "\n", sep="")
cat("Reference generation complete.\n")
