# Installation

Different installation methods are available for **MerKurio**, depending on whether the user has `cargo`/Rust installed and wants to build from source, or only download the executable. 

You can test if the installation was successful by executing MerKurio: 

```bash
# In case you downloaded the standalone executable
./path/to/merkurio --help
# If MerKurio is added to the PATH
merkurio --help
```

## Precompiled Binaries

The most straightforward way with no external dependencies. 

Precompiled binaries for Linux, Windows, and MacOS are available on the [releases page](https://github.com/lschoenm/MerKurio/releases). Download the appropriate binary for your system, decompress and untar it. Optionally, add it to your PATH. The "musl" Linux version can be slower on some systems than the "gnu" version, but should be compatible with a wider range of systems.

To decompress the gzipped tarball, run the following command:

```bash
tar -xzf path/to/release.tar.gz
```

Then, add user permissions (`chmod u+x path/to/merkurio`) to be able to execute the program from this path. 

## Install a Release using `cargo`

Alternatively, you can install a released version of **MerKurio** from source using `cargo`. This requires you to have Rust installed, [which is described here](https://doc.rust-lang.org/cargo/getting-started/installation.html). To install a specific version of MerKurio, run the following command (replace `X.X.X` with the version number):

```bash
cargo install --git https://github.com/lschoenm/MerKurio --tag vX.X.X
```

### Install from Source

To install from source, download the source code from the [releases page](https://github.com/lschoenm/MerKurio/releases): 

```bash
cargo install --path path/to/source/
```

Or clone the repo to get the latest version, decompress and untar it, and install MerKurio using `cargo`. This will automatically add it to your path:

```bash
git clone https://github.com/lschoenm/MerKurio
cd MerKurio
cargo install --path .
```

### Build from Source

In order to build the program from source without having it added to your PATH, clone the repository and run the following command in the root directory. The binary will be located in the `target/release` directory:

```bash
cargo build --release
```
