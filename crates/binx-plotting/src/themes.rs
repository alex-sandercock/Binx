//! Color themes for GWAS plots

use plotters::style::RGBColor;

/// Color theme for plots
#[derive(Debug, Clone)]
pub struct Theme {
    /// Background color
    pub background: RGBColor,
    /// Text color
    pub text: RGBColor,
    /// Axis color
    pub axis: RGBColor,
    /// Significance threshold line color
    pub significance_line: RGBColor,
    /// Suggestive threshold line color
    pub suggestive_line: RGBColor,
    /// Reference line color (e.g., diagonal in QQ plot)
    pub reference_line: RGBColor,
    /// Confidence band color
    pub confidence_band: RGBColor,
    /// Alternating chromosome colors
    pub chromosome_colors: Vec<RGBColor>,
}

impl Default for Theme {
    fn default() -> Self {
        Self::classic()
    }
}

impl Theme {
    /// Classic Manhattan plot theme with blue/orange alternating chromosomes
    pub fn classic() -> Self {
        Self {
            background: RGBColor(255, 255, 255),
            text: RGBColor(0, 0, 0),
            axis: RGBColor(100, 100, 100),
            significance_line: RGBColor(255, 0, 0),
            suggestive_line: RGBColor(0, 0, 255),
            reference_line: RGBColor(100, 100, 100),
            confidence_band: RGBColor(200, 200, 200),
            chromosome_colors: vec![
                RGBColor(31, 119, 180),   // Blue
                RGBColor(255, 127, 14),   // Orange
            ],
        }
    }

    /// Nature-style theme with muted colors
    pub fn nature() -> Self {
        Self {
            background: RGBColor(255, 255, 255),
            text: RGBColor(50, 50, 50),
            axis: RGBColor(80, 80, 80),
            significance_line: RGBColor(178, 34, 34),  // Firebrick
            suggestive_line: RGBColor(100, 100, 100),
            reference_line: RGBColor(100, 100, 100),
            confidence_band: RGBColor(220, 220, 220),
            chromosome_colors: vec![
                RGBColor(77, 77, 77),      // Dark gray
                RGBColor(153, 153, 153),   // Light gray
            ],
        }
    }

    /// Colorful theme with distinct chromosome colors
    pub fn colorful() -> Self {
        Self {
            background: RGBColor(255, 255, 255),
            text: RGBColor(0, 0, 0),
            axis: RGBColor(100, 100, 100),
            significance_line: RGBColor(255, 0, 0),
            suggestive_line: RGBColor(128, 128, 128),
            reference_line: RGBColor(100, 100, 100),
            confidence_band: RGBColor(200, 200, 200),
            chromosome_colors: vec![
                RGBColor(228, 26, 28),    // Red
                RGBColor(55, 126, 184),   // Blue
                RGBColor(77, 175, 74),    // Green
                RGBColor(152, 78, 163),   // Purple
                RGBColor(255, 127, 0),    // Orange
                RGBColor(255, 255, 51),   // Yellow
                RGBColor(166, 86, 40),    // Brown
                RGBColor(247, 129, 191),  // Pink
            ],
        }
    }

    /// Dark theme for presentations
    pub fn dark() -> Self {
        Self {
            background: RGBColor(30, 30, 30),
            text: RGBColor(220, 220, 220),
            axis: RGBColor(150, 150, 150),
            significance_line: RGBColor(255, 100, 100),
            suggestive_line: RGBColor(100, 150, 255),
            reference_line: RGBColor(150, 150, 150),
            confidence_band: RGBColor(70, 70, 70),
            chromosome_colors: vec![
                RGBColor(102, 194, 165),  // Teal
                RGBColor(252, 141, 98),   // Coral
            ],
        }
    }

    /// High contrast theme for accessibility
    pub fn high_contrast() -> Self {
        Self {
            background: RGBColor(255, 255, 255),
            text: RGBColor(0, 0, 0),
            axis: RGBColor(0, 0, 0),
            significance_line: RGBColor(0, 0, 0),
            suggestive_line: RGBColor(100, 100, 100),
            reference_line: RGBColor(0, 0, 0),
            confidence_band: RGBColor(200, 200, 200),
            chromosome_colors: vec![
                RGBColor(0, 0, 0),         // Black
                RGBColor(150, 150, 150),   // Gray
            ],
        }
    }
}
