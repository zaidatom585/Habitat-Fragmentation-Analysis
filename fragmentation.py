"""
Habitat fragmentation metrics from FWS critical habitat KML files.

Pipeline for one species:
  1. Stream the KML file and read every polygon (outer ring + holes).
  2. Compute geodesic polygon areas on a spherical Earth.
  3. Merge polygons whose edges lie within `merge_km` of each other into
     habitat patches (fragments).
  4. Group patches into clusters: patches within `cluster_km` of each other
     belong to the same cluster (single-linkage).
  5. Measure isolation (distance from each patch to its nearest neighbor)
     and extent (edge-to-edge distance between the two farthest patches).

Distances are computed on 3D points on a sphere with k-d trees, so they are
great-circle distances rather than distances on a flat map projection.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

EARTH_RADIUS_KM = 6371.0088      # mean Earth radius
KML_NS = "{http://www.opengis.net/kml/2.2}"


# ---------------------------------------------------------------------------
# Reading KML
# ---------------------------------------------------------------------------

@dataclass
class Polygon:
    outer: np.ndarray            # (n, 2) array of lon, lat in degrees
    holes: list[np.ndarray]      # list of (m, 2) arrays


def _parse_coordinates(text: str) -> np.ndarray:
    """Convert a KML <coordinates> string into an (n, 2) lon/lat array."""
    pairs = [token.split(",")[:2] for token in text.split()]
    return np.asarray(pairs, dtype=float)


def read_kml_polygons(path: str) -> list[Polygon]:
    """Stream a KML file and return all polygons it contains.

    Uses iterparse and clears processed elements, so large files
    (the bull trout layer is over 100 MB) do not need to fit in memory
    as a full XML tree.
    """
    polygons: list[Polygon] = []
    for _, elem in ET.iterparse(path, events=("end",)):
        if elem.tag != f"{KML_NS}Polygon":
            continue
        outer_el = elem.find(f"{KML_NS}outerBoundaryIs/{KML_NS}LinearRing/{KML_NS}coordinates")
        if outer_el is None or not (outer_el.text or "").strip():
            elem.clear()
            continue
        holes = [
            _parse_coordinates(h.text)
            for h in elem.findall(f"{KML_NS}innerBoundaryIs/{KML_NS}LinearRing/{KML_NS}coordinates")
            if (h.text or "").strip()
        ]
        polygons.append(Polygon(_parse_coordinates(outer_el.text), holes))
        elem.clear()
    return polygons


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def ring_area_km2(ring: np.ndarray) -> float:
    """Geodesic area of a closed lon/lat ring on a sphere, in km².

    Uses the spherical-excess approximation from Chamberlain & Duquette
    (2007, "Some Algorithms for Polygons on a Sphere"), which is accurate
    for polygons of this size.
    """
    lon = np.radians(ring[:, 0])
    lat = np.radians(ring[:, 1])
    lon2, lat2 = np.roll(lon, -1), np.roll(lat, -1)
    dlon = (lon2 - lon + np.pi) % (2 * np.pi) - np.pi
    total = np.sum(dlon * (2 + np.sin(lat) + np.sin(lat2)))
    return abs(total) * EARTH_RADIUS_KM ** 2 / 2


def polygon_area_km2(poly: Polygon) -> float:
    return max(ring_area_km2(poly.outer) - sum(ring_area_km2(h) for h in poly.holes), 0.0)


def to_xyz(lonlat: np.ndarray) -> np.ndarray:
    """Convert lon/lat degrees to 3D points (km) on the sphere."""
    lon = np.radians(lonlat[:, 0])
    lat = np.radians(lonlat[:, 1])
    cos_lat = np.cos(lat)
    return EARTH_RADIUS_KM * np.column_stack((cos_lat * np.cos(lon), cos_lat * np.sin(lon), np.sin(lat)))


def chord_to_arc_km(chord: np.ndarray | float) -> np.ndarray | float:
    """Convert straight-line (chord) distance through the Earth to great-circle distance."""
    ratio = np.clip(np.asarray(chord) / (2 * EARTH_RADIUS_KM), 0, 1)
    return 2 * EARTH_RADIUS_KM * np.arcsin(ratio)


def densify(ring_xyz: np.ndarray, spacing_km: float) -> np.ndarray:
    """Add points along long edges so no gap between boundary points exceeds spacing_km.

    Vertex-to-vertex distances then approximate true edge-to-edge distances
    to within about half of spacing_km.
    """
    start, end = ring_xyz[:-1], ring_xyz[1:]
    seg_len = np.linalg.norm(end - start, axis=1)
    steps = np.maximum(np.ceil(seg_len / spacing_km).astype(int), 1)
    pieces = [start[steps == 1]]
    for n in np.unique(steps[steps > 1]):
        idx = np.where(steps == n)[0]
        t = np.arange(n)[None, :, None] / n
        pts = start[idx, None, :] + t * (end[idx] - start[idx])[:, None, :]
        pieces.append(pts.reshape(-1, 3))
    pts = np.vstack(pieces + [ring_xyz[-1:]])
    # project back onto the sphere surface
    return pts * (EARTH_RADIUS_KM / np.linalg.norm(pts, axis=1))[:, None]


class _UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, i: int) -> int:
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, i: int, j: int) -> None:
        ri, rj = self.find(i), self.find(j)
        if ri != rj:
            self.parent[rj] = ri

    def groups(self) -> list[list[int]]:
        out: dict[int, list[int]] = {}
        for i in range(len(self.parent)):
            out.setdefault(self.find(i), []).append(i)
        return list(out.values())


def grid_downsample(points: np.ndarray, cell_km: float) -> np.ndarray:
    """Keep one point per 3D grid cell of size cell_km (error <= cell_km * sqrt(3) / 2)."""
    keys = np.floor(points / cell_km).astype(np.int64)
    _, idx = np.unique(keys, axis=0, return_index=True)
    return points[np.sort(idx)]


class _Shapes:
    """Boundary point sets at two resolutions, for fast and accurate distance queries.

    Three levels of precision keep large layers fast:
      1. bounding spheres give cheap lower/upper bounds for every pair;
      2. a coarse copy of each boundary refines the estimate for candidate pairs;
      3. the fine boundary gives the precise distance, computed only where needed.
    """

    def __init__(self, point_sets: list[np.ndarray], fine_km: float, coarse_km: float):
        self.fine = [grid_downsample(p, fine_km) for p in point_sets]
        self.coarse = [grid_downsample(p, coarse_km) for p in point_sets]
        self.fine_trees = [cKDTree(p) for p in self.fine]
        self.coarse_trees = [cKDTree(p) for p in self.coarse]
        self.coarse_err = coarse_km * np.sqrt(3)
        centers = np.array([p.mean(axis=0) for p in point_sets])
        radii = np.array([np.linalg.norm(p - c, axis=1).max() for p, c in zip(point_sets, centers)])
        sq = np.einsum("ij,ij->i", centers, centers)
        d = np.sqrt(np.maximum(sq[:, None] + sq[None, :] - 2 * centers @ centers.T, 0))
        r = radii[:, None] + radii[None, :]
        self.lower = np.maximum(d - r, 0)
        self.upper = d + r
        self._fine_cache: dict[tuple[int, int], float] = {}
        self._coarse_cache: dict[tuple[int, int], float] = {}

    def __len__(self) -> int:
        return len(self.fine)

    @staticmethod
    def _min_dist(sets, trees, cache, i, j) -> float:
        key = (i, j) if i < j else (j, i)
        if key not in cache:
            a, b = (i, j) if len(sets[i]) >= len(sets[j]) else (j, i)
            dist, _ = trees[a].query(sets[b], k=1)
            cache[key] = float(dist.min())
        return cache[key]

    def coarse_distance(self, i: int, j: int) -> float:
        return self._min_dist(self.coarse, self.coarse_trees, self._coarse_cache, i, j)

    def distance(self, i: int, j: int) -> float:
        """Precise minimum boundary-to-boundary chord distance between shapes i and j."""
        return self._min_dist(self.fine, self.fine_trees, self._fine_cache, i, j)

    def linked_pairs(self, threshold: float):
        """Yield pairs (i, j) whose precise distance is <= threshold."""
        ci, cj = np.where(np.triu(self.lower <= threshold, k=1))
        err = self.coarse_err
        for i, j in zip(ci, cj):
            c = self.coarse_distance(i, j)
            if c - err > threshold:
                continue
            if c + err <= threshold or self.distance(i, j) <= threshold:
                yield i, j

    def nearest_neighbor(self, i: int) -> float:
        order = np.argsort(self.lower[i])
        err = self.coarse_err
        best_coarse = np.inf
        candidates = []
        for j in order:
            if j == i:
                continue
            if self.lower[i, j] > best_coarse + err:
                break
            c = self.coarse_distance(i, j)
            candidates.append((c, j))
            best_coarse = min(best_coarse, c)
        return min(self.distance(i, j) for c, j in candidates if c - err <= best_coarse + err)

    def farthest_pair_distance(self) -> float:
        n = len(self)
        iu, ju = np.triu_indices(n, k=1)
        if len(iu) == 0:
            return 0.0
        order = np.argsort(-self.upper[iu, ju])
        err = self.coarse_err
        best_coarse = 0.0
        candidates = []
        for k in order:
            i, j = iu[k], ju[k]
            if self.upper[i, j] < best_coarse - err:
                break
            c = self.coarse_distance(i, j)
            candidates.append((c, i, j))
            best_coarse = max(best_coarse, c)
        return max(self.distance(i, j) for c, i, j in candidates if c + err >= best_coarse - err)


# ---------------------------------------------------------------------------
# Fragmentation metrics
# ---------------------------------------------------------------------------

def analyze_species(
    path: str,
    merge_km: float = 0.5,
    cluster_km: float = 10.0,
    cluster_thresholds: tuple[float, ...] = (5, 10, 25, 50),
    min_patch_km2: float = 0.01,
    fine_km: float = 0.25,
    coarse_km: float = 2.0,
) -> dict:
    """Compute fragmentation metrics for one species' critical habitat KML.

    merge_km:        polygons closer than this are merged into one patch (fragment)
    cluster_km:      patches closer than this belong to the same cluster
    min_patch_km2:   patches smaller than this are treated as digitizing slivers
    fine_km:         boundary resolution for precise distances (error ~0.2 km)
    coarse_km:       resolution of the first-pass distance estimate
    """
    polygons = read_kml_polygons(path)
    areas = np.array([polygon_area_km2(p) for p in polygons])

    # Boundary points include holes, so an "island" polygon sitting inside
    # another polygon's hole is recognized as touching it.
    rings = [
        np.vstack([densify(to_xyz(r), fine_km) for r in [p.outer, *p.holes]])
        for p in polygons
    ]

    # 1. Merge touching or nearly touching polygons into patches.
    poly_shapes = _Shapes(rings, fine_km, coarse_km)
    uf = _UnionFind(len(polygons))
    for i, j in poly_shapes.linked_pairs(merge_km):
        uf.union(i, j)
    patch_members = uf.groups()
    patch_areas = np.array([areas[m].sum() for m in patch_members])

    # Drop digitizing slivers.
    keep = patch_areas >= min_patch_km2
    slivers = int((~keep).sum())
    patch_members = [m for m, k in zip(patch_members, keep) if k]
    patch_areas = patch_areas[keep]

    patches = _Shapes([np.vstack([rings[k] for k in m]) for m in patch_members], fine_km, coarse_km)
    n = len(patches)

    # 2. Isolation: distance from each patch to its nearest neighbor.
    nn = np.array([patches.nearest_neighbor(i) for i in range(n)]) if n > 1 else np.array([np.nan])

    # 3. Clusters at several distance thresholds (single-linkage).
    thresholds = sorted(set(cluster_thresholds) | {cluster_km})
    clusters = {}
    cluster_labels = np.zeros(n, dtype=int)
    for t in thresholds:
        uf_c = _UnionFind(n)
        for i, j in patches.linked_pairs(t):
            uf_c.union(i, j)
        groups = uf_c.groups()
        clusters[t] = len(groups)
        if t == cluster_km:
            for label, members in enumerate(groups):
                cluster_labels[members] = label

    # 4. Extent: edge-to-edge distance between the two farthest-apart patches.
    extent = patches.farthest_pair_distance() if n > 1 else 0.0

    return {
        "polygons": len(polygons),
        "patches": n,
        "slivers_removed": slivers,
        **{f"clusters_{t:g}km": clusters[t] for t in thresholds},
        "total_area_km2": float(patch_areas.sum()),
        "largest_patch_km2": float(patch_areas.max()),
        "smallest_patch_km2": float(patch_areas.min()),
        "median_patch_km2": float(np.median(patch_areas)),
        "largest_patch_share": float(patch_areas.max() / patch_areas.sum()),
        "median_nn_km": float(chord_to_arc_km(np.nanmedian(nn))) if n > 1 else float("nan"),
        "max_nn_km": float(chord_to_arc_km(np.nanmax(nn))) if n > 1 else float("nan"),
        "extent_km": float(chord_to_arc_km(extent)),
        "patch_areas": patch_areas,
        "patch_members": patch_members,
        "cluster_labels": cluster_labels,
        "polygons_lonlat": [p.outer for p in polygons],
    }
