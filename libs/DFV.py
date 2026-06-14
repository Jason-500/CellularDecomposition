from sage.all import *
import copy as pycopy
from libs.functions import *


class DFV:
    """
    class DFV: for managing Delaunay with faces and dual(Voronoi).

    D: Delaunay graph;

    F: Faces of Delaunay;
    
    V: Voronoi dual of D.
    """
    id = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
    
    def __init__(self,D,F,V,is_interior=False,is_canonical=False):
        """
        class DFV: for managing Delaunay with faces and dual(Voronoi).
        D: Delaunay graph;
        F: Faces of Delaunay;
        V: Voronoi dual of D.
        """
        self.D = D
        self.F = F
        self.V = V
        self.card_E = len(D.edges())
        self.card_F = len(F)
        self._poly = None
        self._dim = None
        self._faces = None

        self.is_interior = None
        if is_interior: 
            self.is_interior = True

        self.is_canonical = None
        self.canonical_image = None
        if is_canonical: 
            self.is_canonical = True
            self.canonical_image = (self,DFV.id)

        # self.edge_dict = *under construction*
        
        
    def remove_single_edge(self,i):
        """
        Return the subgraphs obtained by removing ith edge from given Delaunay graph D, with modification of F and V.
        """
        d = self.D.copy()
        ei = d.edges()[i]
        d.delete_edge(ei)
        f,v = modified_faces_with_voronoi(ei,self.F,self.V)
        return DFV(d,f,v)

    def remove_multiple_edges(self,indices,is_interior=True):
        d, f, v = [pycopy.deepcopy(x) for x in (self.D, self.F, self.V)]
        e_list = self.D.edges()
        for i in indices:
            d.delete_edge(e_list[i])
            f,v = modified_faces_with_voronoi(e_list[i],f,v)
        return DFV(d,f,v,is_interior)
    
    @property
    def poly(self):
        if self._poly: return self._poly
        self._poly = poly(self.D,self.F)
        self._dim = self.poly.dim()
        return self._poly
    
    @property
    def dim(self):
        return self.poly.dim()

    @property
    def faces(self):
        """
        Return the faces: [interiors, boundaries]
        """
        if self._faces: return self._faces
        raw_interior_faces = []
        raw_boundary_faces = []
        voronoi_edge_dict = {label:(start,end) for start,end,label in self.V.edges()}
        for face in self.poly.facets():
            center = face.as_polyhedron().center()
            index_of_collapsed_edges = [i for i, x in enumerate(center) if x == 1]
            labels_of_collapsed_edges = [self.D.edges()[i][2] for i in index_of_collapsed_edges]

            if 0 in center: # non-separating degeneration
                raw_boundary_faces.append(self.remove_multiple_edges(index_of_collapsed_edges,is_interior=False))

            elif not Graph([voronoi_edge_dict[label] for label in labels_of_collapsed_edges],multiedges=True, loops=True).is_forest(): 
                # separating degeneration. I think it will always occur by pinching a single selfloop in V. But first we still detect the loop.
                raw_boundary_faces.append(self.remove_multiple_edges(index_of_collapsed_edges,is_interior=False))

            else:
                raw_interior_faces.append(self.remove_multiple_edges(index_of_collapsed_edges,is_interior=True))
        self._faces=[raw_interior_faces,raw_boundary_faces]
        return self._faces
        
    def return_DFV_tuple(self):
        return (self.D,self.F,self.V)

    def is_orientation_preserving_isomorphic_to(self, dfv):
        return is_orientation_preserving_isomorphic(self, dfv)

    def set_canonical_image(self,dfv,iso):
        self.is_canonical = False
        self.canonical_image = (dfv,iso)

    def set_canonical(self):
        self.is_canonical = True
        self.canonical_image = (self,DFV.id)

    def automorphism_group(self):
        original_grp = self.D.automorphism_group()
        canonical_faces = [canonical_face(f) for f in self.F]
        group_element_list = []
        for g in original_grp:
            current_map_dict = {v: g(v) for v in self.D.vertices()}
            mapping_func = lambda a: current_map_dict[a]
            is_valid = True
            for fi in self.F:
                if canonical_image(fi,mapping_func) not in canonical_faces:
                    is_valid = False
                    break
            if is_valid:
                group_element_list.append(g)

        return PermutationGroup(group_element_list)