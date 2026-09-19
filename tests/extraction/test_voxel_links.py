import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('voxel_link_audit', ROOT/'dev/voxel_link_audit.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_octahedral_patterns_and_closed_cube_contacts():
    assert MODULE.pattern_kind(0) == 'empty'
    assert MODULE.pattern_kind(255) == 'sphere'
    assert MODULE.pattern_kind(1) == 'disk'
    assert MODULE.pattern_kind(3) == 'disk'
    assert MODULE.pattern_kind(0x09) == 'nonmanifold'
    assert MODULE.pattern_kind(0x81) == 'nonmanifold'
    for shift, expected in (((1,0,0),True), ((1,1,0),False), ((1,1,1),False)):
        mask = np.zeros((4,4,4),dtype=bool)
        mask[1,1,1] = True
        mask[tuple(1+s for s in shift)] = True
        assert MODULE.audit_mask(mask)['is_pl_3_manifold_with_boundary'] is expected


def test_solid_shell_and_disconnected_balls_are_manifolds():
    for kind in ('block','shell','separate'):
        mask = np.zeros((7,7,7),dtype=bool)
        mask[1:6,1:6,1:6] = True
        if kind == 'shell':
            mask[2:5,2:5,2:5] = False
        if kind == 'separate':
            mask[:] = False
            mask[1,1,1] = mask[5,5,5] = True
        assert MODULE.audit_mask(mask)['is_pl_3_manifold_with_boundary']
    empty = MODULE.audit_mask(np.zeros((3,3,3),dtype=bool))
    assert empty['empty'] and not empty['is_pl_3_manifold_with_boundary']
