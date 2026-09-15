"""Indexed exact observation map on a full-real-plane interval partition.

All detection thresholds must be record edges. Bins in Z are (lo, hi].
A singleton at an edge has the same observation signature as the adjacent
left interval under Z <= C, so this represents the corresponding quotient
of the full singleton/open-interval construction. No true support is used.
Internal bin densities are displays of bin probabilities, NOT a fitted
continuous within-bin density. Empty possible records are retained.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass
class GridModel:
    z_edges: tuple[np.ndarray, np.ndarray]
    limits_z: np.ndarray
    selection: str = 'union'

    def __post_init__(self):
        self.z_edges = tuple(np.asarray(e, dtype=float) for e in self.z_edges)
        self.limits_z = np.asarray(self.limits_z, dtype=float)
        if self.selection not in ('union', 'both'):
            raise ValueError('selection must be union or both')
        if self.limits_z.ndim != 2 or self.limits_z.shape[1] != 2:
            raise ValueError('limits must have shape (J, 2)')
        if any(e.ndim != 1 or np.any(np.diff(e) <= 0) or not np.isfinite(e).all() for e in self.z_edges):
            raise ValueError('Strictly increasing finite record edges required')
        # An aligned-threshold specialization, never silently quantize a limit.
        cut = np.empty_like(self.limits_z, dtype=int)
        for j, c in enumerate(self.limits_z):
            for b in (0, 1):
                hit = np.flatnonzero(self.z_edges[b] == c[b])
                if len(hit) != 1:
                    raise ValueError('Every threshold must exactly equal a declared record edge')
                cut[j, b] = hit[0]
        self.shape = tuple(len(e) + 1 for e in self.z_edges)
        self.H = int(np.prod(self.shape))
        self.J = len(self.limits_z)
        self.cell_bins = np.stack(np.meshgrid(np.arange(self.shape[0]), np.arange(self.shape[1]), indexing='ij'), axis=-1).reshape(-1, 2)
        self.labels, self.record_index, self.visible = [], [], []
        for j in range(self.J):
            delta = self.cell_bins <= cut[j]
            included = delta.any(axis=1) if self.selection == 'union' else delta.all(axis=1)
            raw = np.column_stack((delta.astype(int), np.where(delta, self.cell_bins, -1)))
            labels = np.unique(raw[included], axis=0)
            lookup = {tuple(r): k for k, r in enumerate(labels)}
            idx = np.full(self.H, -1, dtype=int)
            for h in np.flatnonzero(included):
                idx[h] = lookup[tuple(raw[h])]
            self.labels.append(labels)
            self.record_index.append(idx)
            self.visible.append(included)
        self.visible = np.asarray(self.visible)
        self.invisible = ~self.visible.any(axis=0)
        self.common_detected = np.all(self.cell_bins <= cut[-1], axis=1)
        self._lookup = [{tuple(r): k for k, r in enumerate(rows)} for rows in self.labels]

    def probabilities(self, mass: np.ndarray):
        mass = self.validate_mass(mass)
        a, q = [], []
        for j, idx in enumerate(self.record_index):
            vis = self.visible[j]
            aj = np.bincount(idx[vis], weights=mass[vis], minlength=len(self.labels[j]))
            a.append(aj)
            q.append(float(aj.sum()))
        return a, np.asarray(q)

    def validate_mass(self, mass):
        mass = np.asarray(mass, dtype=float)
        if mass.shape != (self.H,) or not np.isfinite(mass).all() or mass.min() < 0:
            raise ValueError('A finite nonnegative mass for every declared cell is required')
        if abs(mass.sum()-1) > 5e-10:
            raise ValueError('Mass is not normalized')
        return mass

    def record_catalogue(self, z: np.ndarray, field: np.ndarray):
        """Simulation/observation boundary. Never called inside the estimator."""
        z = np.asarray(z, dtype=float)
        field = np.asarray(field, dtype=int)
        if z.ndim != 2 or z.shape[1] != 2 or field.shape != (len(z),):
            raise ValueError('Expected z (N,2), field (N,)')
        if not np.isfinite(z).all() or np.any(field < 0) or np.any(field >= self.J):
            raise ValueError('Invalid catalogue')
        delta = z <= self.limits_z[field]
        inc = delta.any(axis=1) if self.selection == 'union' else delta.all(axis=1)
        bins = np.column_stack([np.searchsorted(e, z[:, b], side='left') for b, e in enumerate(self.z_edges)])
        labels = np.column_stack((delta.astype(int), np.where(delta, bins, -1)))
        indices = np.full(len(z), -1, dtype=int)
        for j in range(self.J):
            for i in np.flatnonzero(inc & (field == j)):
                indices[i] = self._lookup[j][tuple(labels[i])]
        return inc, delta, labels, indices

    def counts_from_indices(self, field, indices, subset=None):
        field, indices = np.asarray(field), np.asarray(indices)
        subset = np.ones(len(field), bool) if subset is None else np.asarray(subset, bool)
        return [np.bincount(indices[subset & (field == j) & (indices >= 0)], minlength=len(self.labels[j])).astype(int) for j in range(self.J)]

    def dense_operators(self):
        return [(idx[None, :] == np.arange(len(rows))[:, None]).astype(int) for idx, rows in zip(self.record_index, self.labels)]

    def as_dict(self):
        return {'z_edges': [e.tolist() for e in self.z_edges], 'limits_z': self.limits_z.tolist(),
                'selection': self.selection, 'shape': list(self.shape),
                'labels': [a.tolist() for a in self.labels], 'record_index': [a.tolist() for a in self.record_index]}

    @classmethod
    def from_dict(cls, d):
        obj = cls(tuple(np.asarray(e) for e in d['z_edges']), np.asarray(d['limits_z']), d['selection'])
        if 'record_index' in d and any(not np.array_equal(x, y) for x, y in zip(obj.record_index, d['record_index'])):
            raise ValueError('Serialized observation map disagrees with declared geometry')
        return obj
