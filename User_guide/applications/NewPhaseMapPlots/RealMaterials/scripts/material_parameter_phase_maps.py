#!/usr/bin/env python3
"""Compatibility entry for the material research scan.

New users: run ``python -m knotted_graph.applications.phase_map_examples --help``
for separate inspect, plot and bounded quick-scan commands. This entry keeps the
research engine's fine-grained options; large scans belong on a compute node.
"""

if __name__ == "__main__":
    from knotted_graph.applications.phase_map_examples._runtime import (
        require_scan_dependencies,
    )

    require_scan_dependencies()
    from knotted_graph.applications.phase_map_examples._materials import main

    main()
