from sage.all import *
from libs.functions import *

class DimensionMismatchError(ValueError):
    """Dimension mismatch."""
    pass   

class NotWellDefinedError(ValueError):
    """Needs to be well defined."""
    pass  


class PolyComplexSolver:
    def __init__(self,initial_dfv_list):
        """
        class PolyComplexModuli: for managing the whole complex.
        """
        self.top_cell_dim = None
        for dfv in initial_dfv_list:
            if self.top_cell_dim is None:
                self.top_cell_dim = dfv.dim
            else:
                if dfv.dim != self.top_cell_dim:
                    raise DimensionMismatchError("Initialization error: inconsistent cell dimensions exist.")
                if not dfv.is_canonical:
                    raise NotWellDefinedError("Warning: cells are not standard representative.")
                if not dfv.is_interior:
                    raise NotWellDefinedError("Warning: cells are not in the interior.")
        self.interior_cell_dict = {self.top_cell_dim:initial_dfv_list}
        self.boundary_cell_dict = {self.top_cell_dim:[]}
        self.boundary_dim = self.top_cell_dim
        
    def solve_boundary(self):
        current_interior_cell = self.interior_cell_dict[self.boundary_dim]
        # current_boundary_cell = self.boundary_cell_dict[self.boundary_dim]
        # Look up all cells, take the canonical one out, set the permutations for faces.
        new_interior_cell_list = []
        new_boundary_cell_list = []
        for dfv in current_interior_cell:
            interior_faces, boundary_faces = dfv.faces

            for interior_face in interior_faces:
                is_canonical = True
                for interior_cell in new_interior_cell_list:
                    is_isomorphic, permutation = interior_cell.is_orientation_preserving_isomorphic_to(interior_face)
                    if is_isomorphic:
                        interior_face.is_canonical = False
                        interior_face.canonical_image = (interior_cell,permutation) # may need to fix the direction of map here.
                        is_canonical = False
                        break
                if is_canonical:
                    interior_face.set_canonical()
                    new_interior_cell_list.append(interior_face)

            for boundary_face in boundary_faces:
                is_canonical = True
                for boundary_cell in new_boundary_cell_list:
                    is_isomorphic, permutation = boundary_cell.is_orientation_preserving_isomorphic_to(boundary_face)
                    if is_isomorphic:
                        boundary_face.is_canonical = False
                        boundary_face.canonical_image = (boundary_cell,permutation) # may need to fix the direction of map here.
                        is_canonical = False
                        break
                if is_canonical:
                    boundary_face.set_canonical()
                    new_boundary_cell_list.append(boundary_face)

        # for dfv in current_boundary_cell:
        #     interior_faces, boundary_faces = dfv.faces

        #     for interior_face in interior_faces:
        #         is_canonical = True
        #         for interior_cell in new_interior_cell_list:
        #             is_isomorphic, permutation = interior_cell.is_orientation_preserving_isomorphic_to(interior_face)
        #             if is_isomorphic:
        #                 interior_face.is_canonical = False
        #                 interior_face.canonical_image = (interior_cell,permutation) # may need to fix the direction of map here.
        #                 is_canonical = False
        #                 break
        #         if is_canonical:
        #             interior_face.set_canonical()
        #             new_interior_cell_list.append(interior_face)

        #     for boundary_face in boundary_faces:
        #         is_canonical = True
        #         for boundary_cell in new_boundary_cell_list:
        #             is_isomorphic, permutation = boundary_cell.is_orientation_preserving_isomorphic_to(boundary_face)
        #             if is_isomorphic:
        #                 boundary_face.is_canonical = False
        #                 boundary_face.canonical_image = (boundary_cell,permutation) # may need to fix the direction of map here.
        #                 is_canonical = False
        #                 break
        #         if is_canonical:
        #             boundary_face.set_canonical()
        #             new_boundary_cell_list.append(boundary_face)
        
        self.boundary_dim -= 1
        self.interior_cell_dict[self.boundary_dim] = new_interior_cell_list
        self.boundary_cell_dict[self.boundary_dim] = new_boundary_cell_list

    def solve_boundary_till_dim(self,dim):
        while self.boundary_dim != dim:
            self.solve_boundary()