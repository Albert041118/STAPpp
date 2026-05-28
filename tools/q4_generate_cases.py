#!/usr/bin/env python3
"""Generate Q4 verification input decks for STAP++.

Generated cases:
- q4_patch_2x2.dat: constant strain patch test using equivalent nodal forces.
- q4_conv_tension_{1,2,4,8}.dat: left-fixed square plate with right-edge horizontal traction.
- q4_conv_shear_{1,2,4,8}.dat: left-fixed square plate with right-edge vertical traction.

The scripts intentionally keep STAP++ input format unchanged: non-zero displacement
boundary conditions are avoided by choosing zero-valued constrained target DOFs in
the patch test and using concentrated nodal loads in all cases.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

E = 210000.0
NU = 0.3
THICKNESS = 1.0
MODE = 1  # plane stress


@dataclass
class Node:
    nid: int
    x: float
    y: float
    fixed_x: bool = False
    fixed_y: bool = False


Element = Tuple[int, int, int, int, int]  # eid, n1, n2, n3, n4
Load = Tuple[int, int, float]  # node, dof(1=x,2=y), value


def elasticity_matrix(E: float = E, nu: float = NU, mode: int = MODE) -> List[List[float]]:
    if mode == 1:
        c = E / (1.0 - nu * nu)
        return [[c, c * nu, 0.0], [c * nu, c, 0.0], [0.0, 0.0, c * (1.0 - nu) / 2.0]]
    c = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return [[c * (1.0 - nu), c * nu, 0.0], [c * nu, c * (1.0 - nu), 0.0], [0.0, 0.0, c * (1.0 - 2.0 * nu) / 2.0]]


def shape_derivatives(xi: float, eta: float) -> Tuple[List[float], List[float]]:
    dxi = [
        -0.25 * (1.0 - eta),
         0.25 * (1.0 - eta),
         0.25 * (1.0 + eta),
        -0.25 * (1.0 + eta),
    ]
    deta = [
        -0.25 * (1.0 - xi),
        -0.25 * (1.0 + xi),
         0.25 * (1.0 + xi),
         0.25 * (1.0 - xi),
    ]
    return dxi, deta


def b_matrix(coords: List[Tuple[float, float]], xi: float, eta: float) -> Tuple[List[List[float]], float]:
    dxi, deta = shape_derivatives(xi, eta)
    dx_dxi = sum(dxi[i] * coords[i][0] for i in range(4))
    dy_dxi = sum(dxi[i] * coords[i][1] for i in range(4))
    dx_deta = sum(deta[i] * coords[i][0] for i in range(4))
    dy_deta = sum(deta[i] * coords[i][1] for i in range(4))
    detj = dx_dxi * dy_deta - dx_deta * dy_dxi
    if detj <= 0.0:
        raise ValueError(f"Non-positive detJ={detj}")
    B = [[0.0 for _ in range(8)] for _ in range(3)]
    for i in range(4):
        dndx = (dy_deta * dxi[i] - dy_dxi * deta[i]) / detj
        dndy = (-dx_deta * dxi[i] + dx_dxi * deta[i]) / detj
        ux = 2 * i
        uy = 2 * i + 1
        B[0][ux] = dndx
        B[1][uy] = dndy
        B[2][ux] = dndy
        B[2][uy] = dndx
    return B, detj


def q4_stiffness(coords: List[Tuple[float, float]]) -> List[List[float]]:
    D = elasticity_matrix()
    K = [[0.0 for _ in range(8)] for _ in range(8)]
    a = 1.0 / math.sqrt(3.0)
    for xi in (-a, a):
        for eta in (-a, a):
            B, detj = b_matrix(coords, xi, eta)
            DB = [[sum(D[i][k] * B[k][j] for k in range(3)) for j in range(8)] for i in range(3)]
            for i in range(8):
                for j in range(8):
                    K[i][j] += sum(B[k][i] * DB[k][j] for k in range(3)) * detj * THICKNESS
    return K


def structured_mesh(nx: int, ny: int, L: float = 1.0, H: float = 1.0) -> Tuple[List[Node], List[Element]]:
    nodes: List[Node] = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            nid = j * (nx + 1) + i + 1
            nodes.append(Node(nid, L * i / nx, H * j / ny))
    elems: List[Element] = []
    for j in range(ny):
        for i in range(nx):
            eid = j * nx + i + 1
            n1 = j * (nx + 1) + i + 1
            n2 = n1 + 1
            n4 = (j + 1) * (nx + 1) + i + 1
            n3 = n4 + 1
            elems.append((eid, n1, n2, n3, n4))
    return nodes, elems


def write_dat(path: Path, title: str, nodes: List[Node], elems: List[Element], loads: List[Load]) -> None:
    lines: List[str] = [title, f"{len(nodes)}\t1\t1\t1"]
    for n in nodes:
        bx = 1 if n.fixed_x else 0
        by = 1 if n.fixed_y else 0
        lines.append(f"{n.nid}\t{bx}\t{by}\t1\t{n.x:.10g}\t{n.y:.10g}\t0.0")
    lines.append("1")
    lines.append(str(len(loads)))
    for nid, dof, val in loads:
        lines.append(f"{nid}\t{dof}\t{val:.12g}")
    lines.append(f"2\t{len(elems)}\t1")
    lines.append(f"1\t{E:.12g}\t{NU:.12g}\t{THICKNESS:.12g}\t{MODE}")
    for eid, n1, n2, n3, n4 in elems:
        lines.append(f"{eid}\t{n1}\t{n2}\t{n3}\t{n4}\t1")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_patch(outdir: Path) -> Dict:
    nx = ny = 2
    nodes, elems = structured_mesh(nx, ny)
    # Linear displacement field. Constrained DOFs are chosen where target value is zero.
    a = 1.0e-3
    b = -5.0e-4
    for n in nodes:
        n.fixed_x = abs(n.x) < 1e-12      # u = a*x = 0 on left edge
        n.fixed_y = abs(n.y) < 1e-12      # v = b*y = 0 on bottom edge

    ndof = 2 * len(nodes)
    Kglob = [[0.0 for _ in range(ndof)] for _ in range(ndof)]
    for _, n1, n2, n3, n4 in elems:
        ids = [n1, n2, n3, n4]
        coords = [(nodes[i - 1].x, nodes[i - 1].y) for i in ids]
        Ke = q4_stiffness(coords)
        lm: List[int] = []
        for nid in ids:
            lm.extend([2 * (nid - 1), 2 * (nid - 1) + 1])
        for i in range(8):
            for j in range(8):
                Kglob[lm[i]][lm[j]] += Ke[i][j]

    u: List[float] = []
    expected_nodes: Dict[str, Dict[str, float]] = {}
    for n in nodes:
        ux = a * n.x
        uy = b * n.y
        u.extend([ux, uy])
        expected_nodes[str(n.nid)] = {"ux": ux, "uy": uy}

    f = [sum(Kglob[i][j] * u[j] for j in range(ndof)) for i in range(ndof)]
    loads: List[Load] = []
    for n in nodes:
        dofs = [(1, n.fixed_x), (2, n.fixed_y)]
        for dof, fixed in dofs:
            idx = 2 * (n.nid - 1) + (dof - 1)
            if not fixed and abs(f[idx]) > 1e-8:
                loads.append((n.nid, dof, f[idx]))

    name = "q4_patch_2x2"
    write_dat(outdir / f"{name}.dat", "Q4_patch_test_2x2_constant_strain", nodes, elems, loads)
    D = elasticity_matrix()
    strain = [a, b, 0.0]
    stress = [sum(D[i][j] * strain[j] for j in range(3)) for i in range(3)]
    return {
        "name": name,
        "type": "patch",
        "expected_nodes": expected_nodes,
        "expected_stress": {"sigma_x": stress[0], "sigma_y": stress[1], "tau_xy": stress[2]},
    }


def edge_loads(nodes: List[Node], nx: int, ny: int, component: int, total: float) -> List[Load]:
    # Consistent enough nodal distribution for a uniform edge traction: half weight at corners.
    right_nodes = [n for n in nodes if abs(n.x - 1.0) < 1e-12]
    loads: List[Load] = []
    for n in right_nodes:
        weight = 0.5 if abs(n.y) < 1e-12 or abs(n.y - 1.0) < 1e-12 else 1.0
        val = total * weight / ny
        loads.append((n.nid, component, val))
    return loads


def generate_convergence(outdir: Path, family: str, n: int) -> Dict:
    nodes, elems = structured_mesh(n, n)
    for node in nodes:
        if abs(node.x) < 1e-12:
            node.fixed_x = True
            node.fixed_y = True
    if family == "tension":
        loads = edge_loads(nodes, n, n, 1, 1000.0)
    elif family == "shear":
        loads = edge_loads(nodes, n, n, 2, 500.0)
    else:
        raise ValueError(family)
    name = f"q4_conv_{family}_{n}x{n}"
    write_dat(outdir / f"{name}.dat", f"Q4_convergence_{family}_{n}x{n}", nodes, elems, loads)
    return {"name": name, "type": "convergence", "family": family, "n": n}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="data/q4_generated")
    parser.add_argument("--manifest", default="data/q4_generated/manifest.json")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cases: List[Dict] = [generate_patch(outdir)]
    for family in ("tension", "shear"):
        for n in (1, 2, 4, 8):
            cases.append(generate_convergence(outdir, family, n))
    manifest = {"material": {"E": E, "nu": NU, "thickness": THICKNESS, "mode": MODE}, "cases": cases}
    Path(args.manifest).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Generated {len(cases)} cases in {outdir}")


if __name__ == "__main__":
    main()
