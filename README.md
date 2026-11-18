# defcell

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**defcell** is a Python command-line tool that uses the [PyMOL](https://pymol.org/) API to deform crystal unit cells and molecular coordinates loaded from a PyMOL session (`.pse`) file, then save the result as PDB.

## Features

- Deforms lattice vectors using a deformation gradient tensor
- Updates unit-cell parameters and the space group
- Moves molecular coordinates consistently with the deformation
- Supports normal and shear strains
- Supports pure-shear and simple-shear deformation modes

Supported strain components:

- Normal strain: `exx`, `eyy`, `ezz`
- Shear strain: `gyz`, `gzx` (`gxz` alias), `gxy`

## Requirements

- Python 3.10 or later
- [PyMOL](https://pymol.org/) Python API
- NumPy

## Installation

Clone the repository:

```bash
git clone https://github.com/wako-k/defcell.git
cd defcell
```

Install NumPy:

```bash
python -m pip install numpy
```

Make sure PyMOL is available in the same Python environment. Installation may vary depending on your platform or whether you use Conda.

## Usage

```bash
python defcell.py <PSE_FILE> --strain <COMPONENT> --amount <STRAIN> [options]
```

### Options

| Option | Description |
|---|---|
| `PSE_FILE` | Input PyMOL session file (`.pse`) |
| `-s`, `--strain` | Strain component |
| `-a`, `--amount` | Strain amount in percent |
| `--shear-mode` | `pure` (default) or `simple` |
| `-o`, `--output` | Output PDB filename (default: `out.pdb`) |

Available strain components:

```text
exx, eyy, ezz, gyz, gzx, gxz, gxy
```

### Examples

Apply 10% normal strain along the x-axis:

```bash
python defcell.py tests/anthracene.pse --strain exx --amount 10 --output anthracene_exx10.pdb
```

Apply 5% pure shear to `gxy`:

```bash
python defcell.py tests/anthracene.pse --strain gxy --amount 5 --shear-mode pure --output anthracene_gxy05.pdb
```

Apply 5% simple shear to `gzx`:

```bash
python defcell.py tests/anthracene.pse --strain gzx --amount 5 --shear-mode simple --output out.pdb
```

Display command help:

```bash
python defcell.py --help
```

## Test Files

The `tests/` directory contains PyMOL session files, PDB files, and shell scripts for testing.

Run the shell script from the repository root:

```bash
bash tests/test.sh
```

## License

This project is licensed under the [MIT License](LICENSE).
