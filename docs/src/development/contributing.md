# Contributing

We welcome contributions to Binx! Here's how to get started.

## Development Setup

### Prerequisites

- Rust toolchain (1.70+)
- Git

### Clone and Build

```bash
git clone https://github.com/alex-sandercock/Binx.git
cd Binx
cargo build
```

### Run Tests

```bash
cargo test
```

### Run with Debug Output

```bash
RUST_LOG=debug cargo run -- gwas --help
```

## Code Style

- Follow Rust conventions
- Run `cargo fmt` before committing
- Run `cargo clippy` to check for issues

## Making Changes

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make your changes
4. Add tests
5. Run tests: `cargo test`
6. Commit: `git commit -m "Add my feature"`
7. Push: `git push origin feature/my-feature`
8. Open a Pull Request

## Areas to Contribute

- Documentation improvements
- New genetic models
- Performance optimizations
- Additional output formats
- Bug fixes

## Reporting Issues

When reporting bugs, include:
- Binx version (`binx --version`)
- Operating system
- Minimal reproducible example
- Expected vs actual behavior

## License

Contributions are licensed under GPL-3.0.

## Contact

- GitHub Issues: [alex-sandercock/Binx](https://github.com/alex-sandercock/Binx/issues)
