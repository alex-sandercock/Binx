# Installation

Binx can be installed via pre-built binaries (recommended) or built from source.

## Pre-built Binaries (Recommended)

Download the latest release for your platform from [GitHub Releases](https://github.com/alex-sandercock/Binx/releases):

| Platform | Download |
|----------|----------|
| Linux (x86_64) | `binx-linux-x86_64.tar.gz` |
| macOS (Intel) | `binx-macos-x86_64.tar.gz` |
| macOS (Apple Silicon) | `binx-macos-aarch64.tar.gz` |

### Linux Installation

```bash
# Download the latest release
curl -LO https://github.com/alex-sandercock/Binx/releases/latest/download/binx-linux-x86_64.tar.gz

# Extract the archive
tar -xzf binx-linux-x86_64.tar.gz

# Verify the installation
./binx --help

# (Optional) Move to a directory in your PATH
mkdir -p ~/bin
mv binx ~/bin/

# Add to PATH if not already present (add to ~/.bashrc or ~/.zshrc)
export PATH="$HOME/bin:$PATH"
```

### macOS Installation

```bash
# For Apple Silicon (M1/M2/M3)
curl -LO https://github.com/alex-sandercock/Binx/releases/latest/download/binx-macos-aarch64.tar.gz
tar -xzf binx-macos-aarch64.tar.gz

# For Intel Macs
curl -LO https://github.com/alex-sandercock/Binx/releases/latest/download/binx-macos-x86_64.tar.gz
tar -xzf binx-macos-x86_64.tar.gz

# Verify the installation
./binx --help
```

> **Note for macOS users**: You may need to allow the binary to run in System Preferences > Security & Privacy if you see a security warning.

## Building from Source

Building from source requires the Rust toolchain (`cargo` + `rustc`).

### Install Rust

If you don't have Rust installed, get it from [rustup.rs](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

### Clone and Build

```bash
# Clone the repository
git clone https://github.com/alex-sandercock/Binx.git
cd Binx

# Build in release mode (optimized)
cargo build --release

# The binary will be at target/release/binx
./target/release/binx --help
```

### Install to Cargo Bin Directory

```bash
# Install globally via cargo
cargo install --path binx-cli

# Now binx is available anywhere
binx --help
```

## Verifying Your Installation

After installation, verify that Binx is working correctly:

```bash
# Check version
binx --version

# View available commands
binx --help

# Check a specific command
binx gwas --help
```

You should see output similar to:

```
binx 0.1.0
Rust command-line genomics workbench for diploid and polyploid species

USAGE:
    binx <COMMAND>

COMMANDS:
    gwas       GWASpoly-style GWAS with multiple genetic models
    kinship    Compute kinship matrix (VanRaden method)
    dosage     Estimate genotype dosages from read counts
    convert    Convert VCF to other formats
    plot       Generate Manhattan, QQ, or LD decay plots
    qtl        Identify significant QTLs from GWAS results
    threshold  Calculate significance thresholds
    help       Print this message or the help of the given subcommand(s)
```

## System Requirements

- **OS**: Linux (x86_64), macOS (Intel or Apple Silicon)
- **RAM**: Depends on dataset size; typically 4-16 GB for standard GWAS
- **Disk**: Minimal for the binary; data storage depends on your datasets

## Troubleshooting

### "Permission denied" error

Make the binary executable:

```bash
chmod +x binx
```

### "Command not found" after installation

Ensure the binary location is in your PATH:

```bash
# Check where binx is located
which binx

# If not found, add the directory to your PATH
export PATH="$HOME/bin:$PATH"  # or wherever you placed the binary
```

### macOS security warning

If macOS blocks the binary:

1. Go to System Preferences > Security & Privacy
2. Click "Allow Anyway" next to the Binx message
3. Run `binx --help` again and click "Open" in the dialog

## Next Steps

Once installed, proceed to the [Quick Start](quickstart.md) guide to run your first analysis.
