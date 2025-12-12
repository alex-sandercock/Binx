//! Output format handling for plots
//!
//! This module provides utilities for saving plots in various formats.
//! The actual rendering is handled by the individual plot modules (manhattan, qq).

/// Supported output formats
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// Scalable Vector Graphics (default)
    Svg,
    /// Portable Network Graphics (requires `png` feature)
    Png,
}

impl OutputFormat {
    /// Detect format from file extension
    pub fn from_extension(ext: &str) -> Option<Self> {
        match ext.to_lowercase().as_str() {
            "svg" => Some(Self::Svg),
            "png" => Some(Self::Png),
            _ => None,
        }
    }

    /// Get the file extension for this format
    pub fn extension(&self) -> &'static str {
        match self {
            Self::Svg => "svg",
            Self::Png => "png",
        }
    }
}

impl Default for OutputFormat {
    fn default() -> Self {
        Self::Svg
    }
}
