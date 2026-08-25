# scqc-panels

[![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/fauserf/scqc-panels/test.yaml?branch=main
[badge-docs]: https://app.readthedocs.org/projects/scqc-panels/badge/

: Standardised QC diagnostic panels for single-cell data, including MALAT1 fraction and ribosomal protein (cytosol) score metrics.

## Getting started

Please refer to the [documentation][],
in particular, the [API documentation][].

## Installation

You need to have Python 3.12 or newer installed on your system.
If you don't have Python installed, we recommend installing [uv][].

We recommend managing dependencies in project-specific virtual environments to avoid dependency conflicts.
This is most convenient using package managers such as [uv][].
Choose from the options below to install scqc-panels:

<!--
1. Add the latest release of `scqc-panels` from [PyPI][] to your `uv` project:

   ```bash
   uv add scqc-panels
   ```

1. Install the latest release into a [standard virtual environment][venv]:

   ```bash
   (after activating your venv)
   pip install scqc-panels
   ```

-->

1. Install the latest development version:

   ```bash
   pip install git+https://github.com/fauserf/scqc-panels.git  # (or `uv add`)
   ```

## Release notes

See the [changelog][].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][].
If you found a bug, please use the [issue tracker][].

## Citation

> t.b.a

[uv]: https://github.com/astral-sh/uv
[scverse discourse]: https://discourse.scverse.org/
[issue tracker]: https://github.com/fauserf/scqc-panels/issues
[tests]: https://github.com/fauserf/scqc-panels/actions/workflows/test.yaml
[documentation]: https://scqc-panels.readthedocs.io
[changelog]: https://scqc-panels.readthedocs.io/page/changelog.html
[api documentation]: https://scqc-panels.readthedocs.io/page/api.html
[pypi]: https://pypi.org/project/scqc-panels
[venv]: https://docs.python.org/3/tutorial/venv.html
