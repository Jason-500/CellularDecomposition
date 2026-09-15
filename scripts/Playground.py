import sys
sys.path.append("/home/ccsghk/CellularDecomposition")

import os
import sys

# print("cwd =", os.getcwd())
# print("file =", __file__)
# print("sys.path[0:3] =", sys.path[0:3])

from libs.DFV import DFV
from libs.Functions import *
from libs.PolyComplexSolver import PolyComplexSolver
from data.DFVdata import *

# taking out 6-cells.
cell_indices_dim6 = [i+1 for i in range(len(D_raw)) if (poly_i := poly(D_raw[i],F_raw[i])).dim()==6 and 
    # Check if the interior points of those cells satisfying the strict inequality constraint. If not, they are Euclidean.
    # Only need to check if there's coords equal to 0 or 1.
    all((xi!=0) & (xi!=1) for xi in poly_i.center())]# and
    # # Check if there's loop in Delaunay: such loop can't exist.
    # all(e[0]!=e[1] for e in D_raw[i].edges())]
cells_dim6 = [
    DFV(
        D_raw[i-1],
        F_raw[i-1],
        V_raw[i-1],
        is_interior=True,
        is_canonical=True,
        is_root_cell=True
    )
    for i in cell_indices_dim6
]

testpoly2 = PolyComplexSolver([cells_dim6[0]])
testpoly2.solve_boundary_till_dim(0)


from sage.all import QQ

# 取出所有真正的 6-dimensional interior top cells
cell_indices_dim6 = []

for i in range(len(D_raw)):
    P = poly(D_raw[i], F_raw[i])

    if (
        P.dim() == 6
        and all(x != 0 and x != 1 for x in P.center())
    ):
        cell_indices_dim6.append(i)

print("number of top cells =", len(cell_indices_dim6))
# 应该得到 10

top_cells = [
    DFV(
        D_raw[i],
        F_raw[i],
        V_raw[i],
        is_interior=True,
        is_canonical=True,
        is_root_cell=True,
    )
    for i in cell_indices_dim6
]

solver = PolyComplexSolver(top_cells)
solver.solve_boundary_till_dim(0)

chi_eff = QQ(0)

for d in range(7):
    cells = solver.interior_cell_dict[d]

    weighted_number = QQ(0)

    for cell in cells:
        A = cell.planar_graph_automorphism_group()
        weighted_number += QQ(1) / A.order()

    contribution = (-1)**d * weighted_number
    chi_eff += contribution

    print(
        "dim =", d,
        "number =", len(cells),
        "weighted =", weighted_number,
        "contribution =", contribution,
    )

print("effective orbifold Euler characteristic =", chi_eff)
