//! set_k.rs: Rust implementation of R/GWASpoly's set.K function
//!
//! This module computes kinship matrices for the polygenic effect in GWAS.
//!
//! ## R/GWASpoly set.K Reference
//! ```r
//! set.K <- function(data, K=NULL, n.core=1, LOCO=NULL) {
//!   # When LOCO = TRUE, K is computed for each chromosome as K=MM',
//!   # where M is the centered genotype matrix (lines x markers),
//!   # and scaled to have unit diagonal.
//!   #
//!   # When LOCO = FALSE, a single K matrix is computed for all markers.
//!
//!   if (LOCO) {
//!     markers <- split(data@map$Marker, data@map$Chrom)
//!   } else {
//!     markers <- list(all=data@map$Marker)
//!   }
//!
//!   f1 <- function(marks, data) {
//!     M <- scale(data@geno[,marks], center=T, scale=F)
//!     K <- tcrossprod(M)
//!     K <- K/mean(diag(K))
//!   }
//!
//!   K <- lapply(markers, f1, data=data)
//! }
//! ```
//!
//! ## makeLOCO Reference
//! ```r
//! makeLOCO <- function(K, exclude) {
//!   # K is list of per-chromosome kinship matrices
//!   # exclude is vector of indices to exclude
//!   n.chr <- length(K)
//!   tmp <- K[[1]] * 0
//!   keep <- setdiff(1:n.chr, exclude)
//!   for (i in keep) {
//!     tmp <- tmp + K[[i]]
//!   }
//!   return(tmp / length(keep))
//! }
//! ```

use anyhow::{anyhow, Result};
use binx_core::{GenotypeMatrixBiallelic, KinshipMatrix};
use ndarray::{Array2, Axis};
use std::collections::HashMap;

/// Compute kinship matrix for a subset of markers using GWASpoly method.
///
/// This matches R/GWASpoly's f1 function:
/// ```r
/// M <- scale(data@geno[,marks], center=T, scale=F)
/// K <- tcrossprod(M)
/// K <- K/mean(diag(K))
/// ```
///
/// # Arguments
/// * `geno` - Full genotype matrix
/// * `marker_indices` - Indices of markers to include
///
/// # Returns
/// Kinship matrix with mean(diag) = 1.0
fn compute_kinship_for_markers(
    geno: &GenotypeMatrixBiallelic,
    marker_indices: &[usize],
) -> Result<Array2<f64>> {
    let n_markers = marker_indices.len();
    let n_samples = geno.sample_ids.len();

    if n_markers == 0 {
        return Err(anyhow!("No markers provided"));
    }

    // Compute column means for selected markers
    // Matches R's: scale(..., center=T, scale=F) which centers by column mean
    let mut col_means = Vec::with_capacity(n_markers);
    for &midx in marker_indices {
        let row = geno.dosages.index_axis(Axis(0), midx);
        let mut sum = 0.0;
        let mut count = 0usize;
        for &v in row.iter() {
            if v.is_finite() {
                sum += v;
                count += 1;
            }
        }
        if count == 0 {
            return Err(anyhow!("Marker {} has all missing dosages", midx));
        }
        col_means.push(sum / (count as f64));
    }

    // Build centered matrix M (samples × markers) for M @ M'
    // Note: geno.dosages is (markers × samples), we need samples × markers
    // This matches R's: M <- scale(data@geno[,marks], center=T, scale=F)
    let mut centered = Array2::<f64>::zeros((n_samples, n_markers));
    for (new_idx, &midx) in marker_indices.iter().enumerate() {
        let mean = col_means[new_idx];
        for j in 0..n_samples {
            let val = geno.dosages[(midx, j)];
            centered[(j, new_idx)] = if val.is_finite() {
                val - mean
            } else {
                0.0 // Missing values set to mean (centered = 0)
            };
        }
    }

    // K = M @ M' (samples × samples)
    // This matches R's: K <- tcrossprod(M)
    let k = centered.dot(&centered.t());

    // Normalize: K = K / mean(diag(K))
    // This matches R's: K <- K/mean(diag(K))
    let diag_mean = (0..n_samples).map(|i| k[(i, i)]).sum::<f64>() / (n_samples as f64);
    if diag_mean <= 0.0 {
        return Err(anyhow!("Kinship diagonal mean is non-positive"));
    }
    let k_normalized = k / diag_mean;

    Ok(k_normalized)
}

/// Compute per-chromosome kinship matrices (GWASpoly set.K with LOCO=TRUE).
///
/// This matches R/GWASpoly's approach:
/// 1. Split markers by chromosome: `markers <- split(data@map$Marker, data@map$Chrom)`
/// 2. For each chromosome, compute K from ONLY that chromosome's markers
/// 3. Each K is normalized to have mean(diag) = 1.0
///
/// # Arguments
/// * `geno` - Genotype matrix with marker metadata containing chromosome info
///
/// # Returns
/// HashMap from chromosome name to kinship matrix (Array2)
pub fn compute_per_chromosome_kinship(
    geno: &GenotypeMatrixBiallelic,
) -> Result<HashMap<String, Array2<f64>>> {
    let meta = geno
        .marker_metadata
        .as_ref()
        .ok_or_else(|| anyhow!("Marker metadata required for per-chromosome kinship"))?;

    // Split markers by chromosome
    // Matches R's: markers <- split(data@map$Marker, data@map$Chrom)
    let mut chrom_markers: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, m) in meta.iter().enumerate() {
        chrom_markers.entry(m.chrom.clone()).or_default().push(idx);
    }

    // Compute K for each chromosome
    // Matches R's: K <- lapply(markers, f1, data=data)
    let mut per_chrom_k: HashMap<String, Array2<f64>> = HashMap::new();
    for (chrom, marker_indices) in chrom_markers {
        let k = compute_kinship_for_markers(geno, &marker_indices)?;
        per_chrom_k.insert(chrom, k);
    }

    Ok(per_chrom_k)
}

/// Combine per-chromosome kinship matrices for LOCO, excluding specified chromosome.
///
/// This matches R/GWASpoly's makeLOCO function:
/// ```r
/// makeLOCO <- function(K, exclude) {
///   n.chr <- length(K)
///   tmp <- K[[1]] * 0
///   keep <- setdiff(1:n.chr, exclude)
///   for (i in keep) {
///     tmp <- tmp + K[[i]]
///   }
///   return(tmp / length(keep))
/// }
/// ```
///
/// # Arguments
/// * `per_chrom_k` - HashMap of per-chromosome kinship matrices
/// * `exclude_chrom` - Chromosome to exclude
///
/// # Returns
/// LOCO kinship matrix (average of all other chromosomes)
pub fn make_loco_kinship(
    per_chrom_k: &HashMap<String, Array2<f64>>,
    exclude_chrom: &str,
) -> Result<Array2<f64>> {
    let chroms: Vec<&String> = per_chrom_k.keys().collect();
    let keep: Vec<&&String> = chroms.iter().filter(|c| c.as_str() != exclude_chrom).collect();

    if keep.is_empty() {
        return Err(anyhow!("No chromosomes remain after excluding {}", exclude_chrom));
    }

    // Get dimensions from first kept chromosome
    let first_k = per_chrom_k
        .get(*keep[0])
        .ok_or_else(|| anyhow!("Missing kinship for chromosome"))?;
    let n = first_k.nrows();

    // Sum all kept K matrices
    // Matches R's: for (i in keep) { tmp <- tmp + K[[i]] }
    let mut sum_k = Array2::<f64>::zeros((n, n));
    for chrom in &keep {
        let k = per_chrom_k
            .get(**chrom)
            .ok_or_else(|| anyhow!("Missing kinship for chromosome {}", chrom))?;
        sum_k = sum_k + k;
    }

    // Average
    // Matches R's: return(tmp / length(keep))
    let loco_k = sum_k / (keep.len() as f64);

    Ok(loco_k)
}

/// Compute LOCO kinship matrices for all chromosomes.
///
/// This is the main entry point that matches R/GWASpoly's set.K(LOCO=TRUE) + makeLOCO usage.
///
/// # Arguments
/// * `geno` - Genotype matrix with marker metadata
///
/// # Returns
/// HashMap from chromosome name to LOCO kinship matrix (for testing markers on that chromosome)
pub fn compute_loco_kinship_gwaspoly(
    geno: &GenotypeMatrixBiallelic,
) -> Result<HashMap<String, KinshipMatrix>> {
    // Step 1: Compute per-chromosome K's (matches set.K with LOCO=TRUE)
    let per_chrom_k = compute_per_chromosome_kinship(geno)?;

    // Step 2: For each chromosome, create LOCO kinship by averaging others (matches makeLOCO)
    let mut loco_kinships: HashMap<String, KinshipMatrix> = HashMap::new();

    for chrom in per_chrom_k.keys() {
        let loco_k = make_loco_kinship(&per_chrom_k, chrom)?;
        loco_kinships.insert(
            chrom.clone(),
            KinshipMatrix {
                sample_ids: geno.sample_ids.clone(),
                matrix: loco_k,
            },
        );
    }

    Ok(loco_kinships)
}

/// Compute single kinship matrix for all markers (set.K with LOCO=FALSE).
///
/// # Arguments
/// * `geno` - Genotype matrix
///
/// # Returns
/// Single kinship matrix
pub fn compute_kinship_all_markers(geno: &GenotypeMatrixBiallelic) -> Result<KinshipMatrix> {
    let n_markers = geno.marker_ids.len();
    let all_indices: Vec<usize> = (0..n_markers).collect();

    let k = compute_kinship_for_markers(geno, &all_indices)?;

    Ok(KinshipMatrix {
        sample_ids: geno.sample_ids.clone(),
        matrix: k,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use binx_core::MarkerMetadata;
    use ndarray::array;

    #[test]
    fn test_per_chromosome_kinship() {
        // Simple test with 2 chromosomes, 2 markers each
        let geno = GenotypeMatrixBiallelic {
            ploidy: 4,
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into()],
            marker_ids: vec!["M1".into(), "M2".into(), "M3".into(), "M4".into()],
            marker_metadata: Some(vec![
                MarkerMetadata { chrom: "chr1".into(), pos: 100.0 },
                MarkerMetadata { chrom: "chr1".into(), pos: 200.0 },
                MarkerMetadata { chrom: "chr2".into(), pos: 100.0 },
                MarkerMetadata { chrom: "chr2".into(), pos: 200.0 },
            ]),
            dosages: array![
                [0.0, 2.0, 4.0],  // M1 on chr1
                [1.0, 2.0, 3.0],  // M2 on chr1
                [4.0, 2.0, 0.0],  // M3 on chr2
                [3.0, 2.0, 1.0],  // M4 on chr2
            ],
        };

        let per_chrom = compute_per_chromosome_kinship(&geno).unwrap();

        assert!(per_chrom.contains_key("chr1"));
        assert!(per_chrom.contains_key("chr2"));

        // Each K should have mean(diag) ≈ 1.0
        let k1 = &per_chrom["chr1"];
        let diag_mean1: f64 = (0..3).map(|i| k1[(i, i)]).sum::<f64>() / 3.0;
        assert!((diag_mean1 - 1.0).abs() < 1e-10);

        let k2 = &per_chrom["chr2"];
        let diag_mean2: f64 = (0..3).map(|i| k2[(i, i)]).sum::<f64>() / 3.0;
        assert!((diag_mean2 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_make_loco() {
        let geno = GenotypeMatrixBiallelic {
            ploidy: 4,
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into()],
            marker_ids: vec!["M1".into(), "M2".into(), "M3".into(), "M4".into()],
            marker_metadata: Some(vec![
                MarkerMetadata { chrom: "chr1".into(), pos: 100.0 },
                MarkerMetadata { chrom: "chr1".into(), pos: 200.0 },
                MarkerMetadata { chrom: "chr2".into(), pos: 100.0 },
                MarkerMetadata { chrom: "chr2".into(), pos: 200.0 },
            ]),
            dosages: array![
                [0.0, 2.0, 4.0],
                [1.0, 2.0, 3.0],
                [4.0, 2.0, 0.0],
                [3.0, 2.0, 1.0],
            ],
        };

        let per_chrom = compute_per_chromosome_kinship(&geno).unwrap();

        // LOCO for chr1 should be just chr2's K (since only 2 chromosomes)
        let loco_chr1 = make_loco_kinship(&per_chrom, "chr1").unwrap();
        let k2 = &per_chrom["chr2"];

        for i in 0..3 {
            for j in 0..3 {
                assert!((loco_chr1[(i, j)] - k2[(i, j)]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_compute_loco_kinship_gwaspoly() {
        let geno = GenotypeMatrixBiallelic {
            ploidy: 4,
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into()],
            marker_ids: vec!["M1".into(), "M2".into(), "M3".into(), "M4".into()],
            marker_metadata: Some(vec![
                MarkerMetadata { chrom: "chr1".into(), pos: 100.0 },
                MarkerMetadata { chrom: "chr1".into(), pos: 200.0 },
                MarkerMetadata { chrom: "chr2".into(), pos: 100.0 },
                MarkerMetadata { chrom: "chr2".into(), pos: 200.0 },
            ]),
            dosages: array![
                [0.0, 2.0, 4.0],
                [1.0, 2.0, 3.0],
                [4.0, 2.0, 0.0],
                [3.0, 2.0, 1.0],
            ],
        };

        let loco_kinships = compute_loco_kinship_gwaspoly(&geno).unwrap();

        assert!(loco_kinships.contains_key("chr1"));
        assert!(loco_kinships.contains_key("chr2"));

        // LOCO for chr1 uses only chr2 markers
        // LOCO for chr2 uses only chr1 markers
        assert_eq!(loco_kinships["chr1"].sample_ids.len(), 3);
        assert_eq!(loco_kinships["chr2"].sample_ids.len(), 3);
    }
}
