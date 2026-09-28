import numpy as np

from fragmentation import chord_to_arc_km, ring_area_km2, to_xyz


def test_one_degree_cell_at_equator():
    # A 1° x 1° cell at the equator is about 12,364 km².
    ring = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]], dtype=float)
    assert abs(ring_area_km2(ring) - 12_364) < 15


def test_area_ignores_winding_direction():
    ring = np.array([[10, 40], [11, 40], [11, 41], [10, 41], [10, 40]], dtype=float)
    assert np.isclose(ring_area_km2(ring), ring_area_km2(ring[::-1]))


def test_great_circle_distance_one_degree_of_latitude():
    # One degree of latitude is about 111.2 km.
    a, b = to_xyz(np.array([[0.0, 0.0], [0.0, 1.0]]))
    assert abs(chord_to_arc_km(np.linalg.norm(a - b)) - 111.2) < 0.2
