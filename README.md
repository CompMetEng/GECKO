<picture>
  <source media="(prefers-color-scheme: dark)" srcset="img/gecko-logo-white.png">
  <img src="img/gecko-logo.png" alt="GECKO" width="400px">
</picture>

[![Documentation](https://readthedocs.org/projects/gecko-docs/badge/?version=latest)](https://gecko-docs.readthedocs.io/en/latest/)
[![GitHub Discussions](https://img.shields.io/github/discussions-search?query=repo%3ACompMetEng%2FGECKO&label=GitHub%20Discussions)](https://github.com/CompMetEng/GECKO/discussions)

**[Documentation](https://gecko-docs.readthedocs.io/en/latest/)** · [Installation](https://gecko-docs.readthedocs.io/en/latest/installation/) · [Guide](https://gecko-docs.readthedocs.io/en/latest/guide/) · [API reference](https://gecko-docs.readthedocs.io/en/latest/api/) · [Migrate](https://gecko-docs.readthedocs.io/en/latest/migrate.html) · [Citations](https://gecko-docs.readthedocs.io/en/latest/references.html)

## About GECKO

The **GECKO** toolbox enhances a **G**enome-scale model to account for **E**nzyme **C**onstraints, using **K**inetics and **O**mics. The resulting enzyme-constrained model (**ecModel**) can be used to perform simulations where enzyme allocation is either drawn from a total protein pool, or constrained by measured protein levels from proteomics data.

GECKO is written in MATLAB. Its Python port, [geckopy](https://github.com/CompMetEng/geckopy), reads and writes the same ecModel YAML files, and the two toolboxes are kept in step by a [parity ledger](https://github.com/CompMetEng/raven-gecko-parity). The documentation covers both.

This is GECKO **4.0.0b1**. For the changes from GECKO 3, see [GECKO 3 → GECKO 4](https://gecko-docs.readthedocs.io/en/latest/gecko3-to-gecko4.html); for the differences between GECKO and geckopy, see [GECKO vs. geckopy](https://gecko-docs.readthedocs.io/en/latest/gecko-to-geckopy.html).

## Get started

1. Install GECKO as described in the [installation instructions](https://gecko-docs.readthedocs.io/en/latest/installation/).
2. Follow [Getting started](https://gecko-docs.readthedocs.io/en/latest/guide/getting-started.html) and [Building an ecModel](https://gecko-docs.readthedocs.io/en/latest/guide/building-ec-model.html) in the guide, or open `tutorials/full_ecModel/protocol.m` or `tutorials/light_ecModel/protocol.m` in MATLAB. Each script reconstructs and analyses an ecModel step by step, in its _full_ or _light_ form ([GECKO light](https://gecko-docs.readthedocs.io/en/latest/guide/gecko-light.html)), and complements the [Nature Protocols](https://doi.org/10.1038/s41596-023-00931-7) paper.

## Documentation and support

The installation instructions, the protocol and the API reference for MATLAB and Python are at [gecko-docs](https://gecko-docs.readthedocs.io/en/latest/). Use [GitHub Discussions](https://github.com/CompMetEng/GECKO/discussions) for questions and comments.

## Cite us

If you have used GECKO in your work, please cite:

> Chen Y, Gustafsson J, Tafur Rangel A, Anton M, Domenzain I, Kittikunapong C, Li F, Yuan L, Nielsen J & Kerkhoven EJ. Reconstruction, simulation and analysis of enzyme-constrained metabolic models using GECKO Toolbox 3.0. *Nature Protocols* **19**, 629-667 (2024). doi:[10.1038/s41596-023-00931-7](https://doi.org/10.1038/s41596-023-00931-7).

## Contributing

Contributions are welcome. Read the [contributing guidelines](https://github.com/CompMetEng/GECKO/blob/develop/.github/CONTRIBUTING.md) first. A change to a function that geckopy also has needs the matching change there; the [parity ledger](https://github.com/CompMetEng/raven-gecko-parity) lists the pairs.
