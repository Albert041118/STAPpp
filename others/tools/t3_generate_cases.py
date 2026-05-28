#!/usr/bin/env python3
"""Generate T3 verification input decks for STAP++."""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

E = 210000.0
NU = 0.3
THICKNESS = 1.0
MODE = 1


@dataclass
class Node:
    nid: int
    x: float
    y: float
    fixed_x: bool = False
    fixed_y: bool = False


Element = Tuple[int, int, int, int]  # eid, n1, n2, n3
Load = Tuple[int, int, float]  # node, dof(1=x,2=y), value


def elasticity_matrix() -> List[List[float]]:
    c = E / (1.0 - NU * NU)
    return [[c, c * NU, 0.0], [c * NU, c, 0.0], [0.0, 0.0, c * (1.0 - NU) / 2.0]]


def t3_b_matrix(coords: List[Tuple[float, float]]) -> Tuple[List[List[float]], float]:
    (x1, y1), (x2, y2), (x3, y3) = coords
    two_area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    area = 0.5 * two_area
    if area <= 0.0:
        raise ValueError(f"Non-positive T3 area={area}")
    b = [y2 - y3, y3 - y1, y1 - y2]
    c = [x3 - x2, x1 - x3, x2 - x1]
    B = [[0.0 for _ in range(6)] for _ in range(3)]
    for i in range(3):
        ux = 2 * i
        uy = 2 * i + 1
        B[0][ux] = b[i] / two_area
        B[1][uy] = c[i] / two_area
        B[2][ux] = c[i] / two_area
        B[2][uy] = b[i] / two_area
    return B, area


def t3_stiffness(coords: List[Tuple[float, float]]) -> List[List[float]]:
    D = elasticity_matrix()
    B, area = t3_b_matrix(coords)
    DB = [[sum(D[i][k] * B[k][j] for k in range(3)) for j in range(6)] for i in range(3)]
    return [
        [sum(B[k][i] * DB[k][j] for k in range(3)) * area * THICKNESS for j in range(6)]
        for i in range(6)
    ]


def structured_t3_mesh(nx: int, ny: int, L: float = 1.0, H: float = 1.0) -> Tuple[List[Node], List[Element]]:
    nodes: List[Node] = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            nodes.append(Node(j * (nx + 1) + i + 1, L * i / nx, H * j / ny))
    elems: List[Element] = []
    eid = 1
    for j in range(ny):
        for i in range(nx):
            n1 = j * (nx + 1) + i + 1
            n2 = n1 + 1
            n4 = (j + 1) * (nx + 1) + i + 1
            n3 = n4 + 1
            elems.append((eid, n1, n2, n3))
            eid += 1
            elems.append((eid, n1, n3, n4))
            eid += 1
    return nodes, elems


def write_dat(path: Path, title: str, nodes: List[Node], elems: List[Element], loads: List[Load]) -> None:
    lines = [title, f"{len(nodes)}\t1\t1\t1"]
    for n in nodes:
        bx = 1 if n.fixed_x else 0
        by = 1 if n.fixed_y else 0
        lines.append(f"{n.nid}\t{bx}\t{by}\t1\t{n.x:.10g}\t{n.y:.10g}\t0.0")
    lines.append("1")
    lines.append(str(len(loads)))
    for nid, dof, val in loads:
        lines.append(f"{nid}\t{dof}\t{val:.12g}")
    lines.append(f"3\t{len(elems)}\t1")
    lines.append(f"1\t{E:.12g}\t{NU:.12g}\t{THICKNESS:.12g}\t{MODE}")
    for eid, n1, n2, n3 in elems:
        lines.append(f"{eid}\t{n1}\t{n2}\t{n3}\t1")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def edge_loads(nodes: List[Node], ny: int, component: int, total: float) -> List[Load]:
    right_nodes = [n for n in nodes if abs(n.x - 1.0) < 1e-12]
    loads: List[Load] = []
    for n in right_nodes:
        weight = 0.5 if abs(n.y) < 1e-12 or abs(n.y - 1.0) < 1e-12 else 1.0
        loads.append((n.nid, component, total * weight / ny))
    return loads


def generate_basic(outdir: Path) -> List[Dict]:
    cases: List[Dict] = []

    nodes = [Node(1, 0.0, 0.0, True, True), Node(2, 1.0, 0.0, False, True), Node(3, 0.0, 1.0, True, False)]
    elems = [(1, 1, 2, 3)]
    write_dat(outdir / "t3-single.dat", "T3_single_element_plane_stress", nodes, elems, [(2, 1, 500.0)])
    cases.append({"name": "t3-single", "type": "basic"})

    nodes, elems = structured_t3_mesh(1, 1)
    for n in nodes:
        if abs(n.x) < 1e-12:
            n.fixed_x = True
            n.fixed_y = True
    write_dat(outdir / "t3-two-element.dat", "T3_two_element_square_plane_stress", nodes, elems, edge_loads(nodes, 1, 1, 1000.0))
    cases.append({"name": "t3-two-element", "type": "basic"})

    nodes = [Node(1, 0.0, 0.0, True, True), Node(2, 1.0, 0.0, False, True), Node(3, 0.0, 1.0, True, False)]
    write_dat(outdir / "t3-invalid-mode.dat", "T3_invalid_material_mode", nodes, [(1, 1, 2, 3)], [(2, 1, 500.0)])
    path = outdir / "t3-invalid-mode.dat"
    lines = path.read_text(encoding="utf-8").splitlines()
    lines = ["1\t210000\t0.3\t1\t9" if line.startswith("1\t210000\t0.3\t1\t") else line for line in lines]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    nodes = [Node(1, 0.0, 0.0, True, True), Node(2, 1.0, 0.0, False, True), Node(3, 0.0, 1.0, True, False)]
    write_dat(outdir / "t3-invalid-area.dat", "T3_invalid_clockwise_nodes", nodes, [(1, 1, 3, 2)], [(2, 1, 500.0)])
    return cases


def generate_patch(outdir: Path) -> Dict:
    nodes, elems = structured_t3_mesh(2, 2)
    a = 1.0e-3
    b = -5.0e-4
    for n in nodes:
        n.fixed_x = abs(n.x) < 1e-12
        n.fixed_y = abs(n.y) < 1e-12

    ndof = 2 * len(nodes)
    Kglob = [[0.0 for _ in range(ndof)] for _ in range(ndof)]
    for _, n1, n2, n3 in elems:
        ids = [n1, n2, n3]
        coords = [(nodes[i - 1].x, nodes[i - 1].y) for i in ids]
        Ke = t3_stiffness(coords)
        lm: List[int] = []
        for nid in ids:
            lm.extend([2 * (nid - 1), 2 * (nid - 1) + 1])
        for i in range(6):
            for j in range(6):
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
        for dof, fixed in ((1, n.fixed_x), (2, n.fixed_y)):
            idx = 2 * (n.nid - 1) + dof - 1
            if not fixed and abs(f[idx]) > 1e-8:
                loads.append((n.nid, dof, f[idx]))

    name = "t3_patch_2x2"
    write_dat(outdir / f"{name}.dat", "T3_patch_test_2x2_constant_strain", nodes, elems, loads)
    D = elasticity_matrix()
    strain = [a, b, 0.0]
    stress = [sum(D[i][j] * strain[j] for j in range(3)) for i in range(3)]
    return {
        "name": name,
        "type": "patch",
        "expected_nodes": expected_nodes,
        "expected_stress": {"sigma_x": stress[0], "sigma_y": stress[1], "tau_xy": stress[2]},
    }


def generate_convergence(outdir: Path, family: str, n: int) -> Dict:
    nodes, elems = structured_t3_mesh(n, n)
    for node in nodes:
        if abs(node.x) < 1e-12:
            node.fixed_x = True
            node.fixed_y = True
    loads = edge_loads(nodes, n, 1 if family == "tension" else 2, 1000.0 if family == "tension" else 500.0)
    name = f"t3_conv_{family}_{n}x{n}"
    write_dat(outdir / f"{name}.dat", f"T3_convergence_{family}_{n}x{n}", nodes, elems, loads)
    return {"name": name, "type": "convergence", "family": family, "n": n}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="data/t3_generated")
    parser.add_argument("--manifest", default="data/t3_generated/manifest.json")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    cases = generate_basic(outdir)
    cases.append(generate_patch(outdir))
    for family in ("tension", "shear"):
        for n in (1, 2, 4, 8):
            cases.append(generate_convergence(outdir, family, n))

    manifest = {"cases": cases}
    Path(args.manifest).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Generated {len(cases)} T3 cases in {outdir}")


if __name__ == "__main__":
    main()
