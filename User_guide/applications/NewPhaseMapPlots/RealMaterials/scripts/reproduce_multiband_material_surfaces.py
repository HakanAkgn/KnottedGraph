#!/usr/bin/env python3
"""Compatibility entry for the original multiband surface reproduction workflow."""

if __name__ == "__main__":
    from knotted_graph.applications.phase_map_examples._runtime import (
        require_scan_dependencies,
    )

    require_scan_dependencies()
    from knotted_graph.applications.phase_map_examples._material_models import main

    main()
