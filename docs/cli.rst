Command Line Interface
======================

The ``ember`` toolkit provides a command-line interface with multiple
subcommands for computing entropy metrics, generating p-values, plotting
summaries, and extracting highly specific or non-specific genes.

This page documents the full CLI using the same argparse parser that the
tool uses internally.

.. argparse::
   :module: ember_py.cli
   :func: create_parser
   :prog: ember
   :nodefault:
