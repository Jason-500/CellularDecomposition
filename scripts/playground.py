import sys
sys.path.append("/home/ccsghk/CellularDecomposition")

import os
import sys

print("cwd =", os.getcwd())
print("file =", __file__)
print("sys.path[0:3] =", sys.path[0:3])

from libs.DFV import DFV
from libs.functions import *
from libs.PolyComplexSolver import PolyComplexSolver
from data.DFVdata import *

# taking out 6-cells.
cell_indices_dim6 = [i+1 for i in range(len(D_raw)) if (poly_i := poly(D_raw[i],F_raw[i])).dim()==6 and 
    # Check if the interior points of those cells satisfying the strict inequality constraint. If not, they are Euclidean.
    # Only need to check if there's coords equal to 0 or 1.
    all((xi!=0) & (xi!=1) for xi in poly_i.center())]# and
    # # Check if there's loop in Delaunay: such loop can't exist.
    # all(e[0]!=e[1] for e in D_raw[i].edges())]
cells_dim6 = [DFV(D_raw[i-1],F_raw[i-1],V_raw[i-1],is_interior=True,is_canonical=True) for i in cell_indices_dim6]

testpoly2 = PolyComplexSolver([cells_dim6[0]])
testpoly2.solve_boundary_till_dim(0)

