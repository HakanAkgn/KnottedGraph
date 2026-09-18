#!/usr/bin/env python3
"""Compatibility entry for the compact-TPMS research scan.

For a bounded first run use the ``phase_map_examples scan tpms`` command.
"""

if __name__ == "__main__":
    from knotted_graph.applications.phase_map_examples._runtime import (
        require_scan_dependencies,
    )

    require_scan_dependencies()
    from knotted_graph.applications.phase_map_examples._tpms import main

    main()
