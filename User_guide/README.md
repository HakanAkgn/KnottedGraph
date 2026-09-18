# Start using KnottedGraph

Begin with the [installation guide](../doc/installation.md) and
[Quick Start](../doc/quickstart.md). The Quick Start uses the base installation
and prints the same Yamada polynomial for an abstract theta graph and its
crossing-free spatial embedding.

| Your next step | Open | What to expect |
| --- | --- | --- |
| Choose a workflow | [Guide index](00_user_guide.ipynb) | A map of inputs, core tools and applications |
| Run the first examples | [Getting started](01_getting_started.ipynb) | Installation checks and a small graph workflow |
| Work with your own data | [Core workflows](02_core_workflows.ipynb) and [Input Handling](../doc/user_guide/input_adapters.md) | Supported input calls, graph checks and projection |
| Understand advanced options | [Advanced and reproduction](03_advanced_and_reproduction.ipynb) | Optional workflows and reproduction boundaries |
| Choose a scientific example | [Applications](applications/README.md) | Six application notebooks and material/TPMS saved results |
| Check correctness evidence | [Yamada sanity checks](benchmarks/01_yamada_sanity_checks.ipynb) and [application regressions](benchmarks/02_application_regression_checks.ipynb) | Focused checks; publication benchmarks have separate resource needs |

Read notebooks on GitHub to inspect their saved content. To execute locally,
install the extras named by the notebook and start JupyterLab from the repository
root. The `notebook` extra supplies JupyterLab; optional scientific workflows
need their corresponding extras too. A complete Python application environment
uses `uv sync --extra all`.

On a cluster, installations, execution and rendering belong on a compute node.
Read each notebook's run-mode and resource notes before **Run All**. Saved
output is evidence of a previous run, and does not show that your environment
has executed the notebook.

When a command fails, use [Troubleshooting](../doc/troubleshooting.md), retain
the error and check the selected package version and optional dependencies.
Write new results to a new output directory; keep the reference results intact.
