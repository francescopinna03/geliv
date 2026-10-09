import json
from pathlib import Path
import numpy as np


def build(seed, n, cells):
    rng = np.random.default_rng(seed)
    cell = np.repeat(np.arange(cells), n // cells)[rng.permutation(n)]
    zc = rng.normal(size=(cells, 2))
    wc = np.column_stack([rng.normal(size=cells), np.ones(cells)])
    v0 = rng.normal(size=n)
    x = zc[cell, 0] + 0.3 * zc[cell, 1] + v0
    y = 0.7 * x + 0.2 * wc[cell, 0] + 0.2 + 0.5 * v0 + rng.normal(size=n)
    return rng, cell, y, x[:, None], wc[cell], np.column_stack([zc[cell], wc[cell]])


class Rows:
    def __init__(self, y, E, W, Z):
        self.y, self.X, self.Z = y, np.column_stack([E, W]), Z
        self.n = len(y)
        self.sx = np.sqrt((self.X ** 2).mean(0))
        self.sz = np.sqrt((Z ** 2).mean(0))

    def stats(self, b, lam, method):
        G = (self.Z / self.sz) * (self.y - self.X @ (b / self.sx))[:, None]
        v = G @ lam
        if method == 0:
            shift = v.max()
            w = np.exp(v - shift)
            sw = w.sum()
            g = G.T @ w / sw
            return shift + np.log(sw / self.n), g, (G.T * w) @ G / sw - np.outer(g, g)
        if method == -1:
            w = 1 / (1 - v)
            return -np.log(1 - v).sum() / self.n, G.T @ w / self.n, (G.T * w ** 2) @ G / self.n
        return (v + 0.5 * v ** 2).sum() / self.n, G.T @ (1 + v) / self.n, G.T @ G / self.n

    def outer(self, b, lam, f, method):
        Xs, Zs = self.X / self.sx, self.Z / self.sz
        t = Zs @ lam
        v = (self.y - Xs @ b) * t
        w = {0: np.exp(v - f), -1: 1 / (1 - v), 1: 1 + v}[method] / self.n
        return Xs.T @ (w * t)

    def vcov(self, b):
        Zs, Xs = self.Z / self.sz, self.X / self.sx
        G = Zs * (self.y - Xs @ b)[:, None]
        D = -Zs.T @ Xs / self.n
        O = G.T @ G / self.n
        return np.linalg.inv(D.T @ np.linalg.solve(O, D)) / self.n / np.outer(self.sx, self.sx)

    def tsls(self):
        ZZ = self.Z.T @ self.Z / np.outer(self.sz, self.sz)
        ZX = self.Z.T @ self.X / np.outer(self.sz, self.sx)
        Zy = self.Z.T @ self.y / self.sz
        A = np.linalg.solve(ZZ, ZX)
        return np.linalg.solve(ZX.T @ A, ZX.T @ np.linalg.solve(ZZ, Zy)) / self.sx

    def ols(self):
        return np.linalg.lstsq(self.X, self.y, rcond=None)[0]


class Cells:
    def __init__(self, y, E, W, Z, cell):
        self.n, self.e = len(y), E.shape[1]
        self.ord = np.argsort(cell, kind="stable")
        cs = cell[self.ord]
        self.start = np.r_[0, np.flatnonzero(np.diff(cs)) + 1]
        self.cnt = np.diff(np.r_[self.start, self.n]).astype(float)
        self.cid = np.repeat(np.arange(len(self.start)), self.cnt.astype(int))
        first = self.ord[self.start]
        self.Zc, self.Xc = Z[first], W[first]
        back = np.empty(self.n, dtype=int)
        back[self.ord] = self.cid
        assert np.array_equal(Z, self.Zc[back]) and np.array_equal(W, self.Xc[back])
        self.y, self.E = y[self.ord], E[self.ord]
        self.sx = np.r_[np.sqrt((self.E ** 2).mean(0)), np.sqrt((self.cnt[:, None] * self.Xc ** 2).sum(0) / self.n)]
        self.sz = np.sqrt((self.cnt[:, None] * self.Zc ** 2).sum(0) / self.n)

    def psum(self, a):
        return np.add.reduceat(a, self.start, axis=0)

    def resid(self, theta):
        return self.y - self.E @ theta[: self.e] - (self.Xc @ theta[self.e:])[self.cid]

    def stats(self, b, lam, method):
        Zs = self.Zc / self.sz
        u = self.resid(b / self.sx)
        v = u * (Zs @ lam)[self.cid]
        if method == 0:
            shift = v.max()
            w = np.exp(v - shift)
            sw = w.sum()
            g = Zs.T @ self.psum(w * u) / sw
            return shift + np.log(sw / self.n), g, (Zs.T * self.psum(w * u ** 2)) @ Zs / sw - np.outer(g, g)
        if method == -1:
            w = 1 / (1 - v)
            return (-np.log(1 - v).sum() / self.n, Zs.T @ self.psum(w * u) / self.n,
                    (Zs.T * self.psum(w ** 2 * u ** 2)) @ Zs / self.n)
        return ((v + 0.5 * v ** 2).sum() / self.n, Zs.T @ self.psum((1 + v) * u) / self.n,
                (Zs.T * self.psum(u ** 2)) @ Zs / self.n)

    def outer(self, b, lam, f, method):
        u = self.resid(b / self.sx)
        t = (self.Zc / self.sz) @ lam
        v = u * t[self.cid]
        w = {0: np.exp(v - f), -1: 1 / (1 - v), 1: 1 + v}[method] / self.n
        gE = (self.E / self.sx[: self.e]).T @ (w * t[self.cid])
        gW = (self.Xc / self.sx[self.e:]).T @ (t * self.psum(w))
        return np.r_[gE, gW]

    def vcov(self, b):
        Zs = self.Zc / self.sz
        u = self.resid(b / self.sx)
        D = -np.column_stack([Zs.T @ self.psum(self.E / self.sx[: self.e]),
                              (Zs.T * self.cnt) @ (self.Xc / self.sx[self.e:])]) / self.n
        O = (Zs.T * self.psum(u ** 2)) @ Zs / self.n
        return np.linalg.inv(D.T @ np.linalg.solve(O, D)) / self.n / np.outer(self.sx, self.sx)

    def tsls(self):
        Zs = self.Zc / self.sz
        ZZ = (Zs.T * self.cnt) @ Zs
        ZX = np.column_stack([Zs.T @ self.psum(self.E / self.sx[: self.e]),
                              (Zs.T * self.cnt) @ (self.Xc / self.sx[self.e:])])
        Zy = Zs.T @ self.psum(self.y)
        A = np.linalg.solve(ZZ, ZX)
        return np.linalg.solve(ZX.T @ A, ZX.T @ np.linalg.solve(ZZ, Zy)) / self.sx

    def ols(self):
        EW = self.psum(self.E).T @ self.Xc
        XX = np.block([[self.E.T @ self.E, EW], [EW.T, (self.Xc.T * self.cnt) @ self.Xc]])
        Xy = np.r_[self.E.T @ self.y, self.Xc.T @ self.psum(self.y)]
        return np.linalg.solve(XX, Xy)

    def to_original(self, w_sorted):
        w = np.empty(self.n)
        w[self.ord] = w_sorted
        return w


def main():
    rng, cell, y, E, W, Z = build(20261009, 1200, 60)
    rows, cells = Rows(y, E, W, Z), Cells(y, E, W, Z, cell)
    report = {"n": 1200, "cells": 60, "checks": {}}
    report["checks"]["scaling"] = float(max(np.abs(rows.sx - cells.sx).max(), np.abs(rows.sz - cells.sz).max()))
    report["checks"]["tsls"] = float(np.abs(rows.tsls() - cells.tsls()).max())
    report["checks"]["ols"] = float(np.abs(rows.ols() - cells.ols()).max())
    b0 = rows.tsls() * rows.sx
    report["checks"]["vcov"] = float(np.abs(rows.vcov(b0) - cells.vcov(b0)).max())
    for method, name in ((0, "et"), (-1, "el"), (1, "cue")):
        worst_inner, worst_outer = 0.0, 0.0
        for _ in range(25):
            b = (np.array([0.71, 0.18, 0.21]) + 0.05 * rng.normal(size=3)) * rows.sx
            lam = 0.02 * rng.normal(size=Z.shape[1])
            v = ((rows.Z / rows.sz) * (rows.y - rows.X @ (b / rows.sx))[:, None]) @ lam
            if v.max() > 0.8:
                continue
            fr, gr, Hr = rows.stats(b, lam, method)
            fc, gc, Hc = cells.stats(b, lam, method)
            worst_inner = max(worst_inner, abs(fr - fc), np.abs(gr - gc).max(), np.abs(Hr - Hc).max())
            worst_outer = max(worst_outer, np.abs(rows.outer(b, lam, fr, method) - cells.outer(b, lam, fc, method)).max())
        report["checks"][name + "_inner_value_gradient_hessian"] = float(worst_inner)
        report["checks"][name + "_envelope_gradient"] = float(worst_outer)
    u_sorted = cells.resid(rows.tsls())
    u_rows = rows.y - rows.X @ rows.tsls()
    report["checks"]["weight_reordering"] = float(np.abs(cells.to_original(np.exp(0.01 * u_sorted)) - np.exp(0.01 * u_rows)).max())
    report["tolerance"] = 1e-10
    report["passed"] = all(v < report["tolerance"] for v in report["checks"].values())
    out = Path(__file__).resolve().parents[1] / "output" / "cells_validation.json"
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
