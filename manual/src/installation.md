# Installation

You can install MerKurio in several ways, depending on your system and whether you have Rust installed.

[1. Precompiled Binaries (No Rust Needed)](#option-1-precompiled-binaries-no-rust-needed)\
[2. Install via Cargo (Requires Rust)](#option-2-install-via-cargo-requires-rust)\
[3. Build Manually Without Installing (Requires Rust)](#option-3-build-manually-without-installing-requires-rust)

After installation, verify if it works by running:

```bash
merkurio --help
```

Or, if you didn't add it to your PATH:

```bash
./path/to/merkurio --help
```

## Option 1: Precompiled Binaries (No Rust Needed)

Download a binary for Linux, Windows, or macOS from the [releases page](https://github.com/lschoenm/MerKurio/releases), then extract the archive:

```bash
tar -xzf path/to/release.tar.gz
```

On Linux/macOS, make it executable if needed:

```bash
chmod u+x path/to/merkurio
```

The `merkurio-x86_64-unknown-linux-musl` is compatible with a wider range of systems but can have worse performance.

## Option 2: Install via Cargo (Requires Rust)

If you have [Rust installed (edition 2024)](https://doc.rust-lang.org/cargo/getting-started/installation.html), the easiest way is:

```bash
cargo install merkurio
```

This pulls the latest version from [crates.io](https://crates.io/crates/merkurio).

To install a tagged release from GitHub:

```bash
cargo install --git https://github.com/lschoenm/MerKurio --tag vX.X.X

```

## Option 3: Build Manually Without Installing (Requires Rust)

```bash
git clone https://github.com/lschoenm/MerKurio
cd MerKurio
cargo build --release
```

The binary will be in `target/release/`.
