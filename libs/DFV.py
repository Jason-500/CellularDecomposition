from sage.all import *
import copy as pycopy
from libs.functions import *

from itertools import combinations

class DFV:
    """
    class DFV: for managing Delaunay with faces and dual(Voronoi).

    D: Delaunay graph;

    F: Faces of Delaunay;
    
    V: Voronoi dual of D.
    """
    vertex_id = {i: i for i in range(1, 7)}
    edge_id = {f"x{i}": f"x{i}" for i in range(1, 13)}
    id = [vertex_id,edge_id]
    
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
        self._facets = None
        self._faces = {}
        self._edge_labels = [label for _,_,label in D.edges()]

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
    def edge_labels(self):
        return self._edge_labels
    
    # def facets(self):
    #     """
    #     Return the facets: [interior faces, boundary faces]
    #     """
    #     if self._facets: return self._facets
    #     raw_interior_faces = []
    #     raw_boundary_faces = []
    #     voronoi_edge_dict = {label:(start,end) for start,end,label in self.V.edges()}
    #     for face in self.poly.facets():
    #         center = face.as_polyhedron().center()
    #         index_of_collapsed_edges = [i for i, x in enumerate(center) if x == 1]
    #         labels_of_collapsed_edges = [self.D.edges()[i][2] for i in index_of_collapsed_edges]

    #         if 0 in center: # non-separating degeneration
    #             raw_boundary_faces.append(self.remove_multiple_edges(index_of_collapsed_edges,is_interior=False))

    #         elif not Graph([voronoi_edge_dict[label] for label in labels_of_collapsed_edges],multiedges=True, loops=True).is_forest(): 
    #             # separating degeneration. I think it will always occur by pinching a single selfloop in V. But first we still detect the loop.
    #             raw_boundary_faces.append(self.remove_multiple_edges(index_of_collapsed_edges,is_interior=False))

    #         else:
    #             raw_interior_faces.append(self.remove_multiple_edges(index_of_collapsed_edges,is_interior=True))
    #     self._facets=[raw_interior_faces,raw_boundary_faces]
    #     return self._facets
    # 
    # Reuse faces():
    def facets(self):
        return self.faces(1)
    
    def faces(self, dim):
        """
        Return the faces with given dimension: [interior faces, boundary faces]
        """
        if dim in self._faces: return self._faces[dim]
        raw_interior_faces = []
        raw_boundary_faces = []
        voronoi_edge_dict = {label:(start,end) for start,end,label in self.V.edges()}
        for face in self.poly.faces(dim):
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
        self._faces[dim]=[raw_interior_faces,raw_boundary_faces]
        return self._faces[dim]

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
    
# methods for determining if isomorphic
def canonical_face(face):
    n = len(face)
    face_nolabel = [(e[0],e[1]) for e in face]
    rotations = [tuple(face_nolabel[i:] + face_nolabel[:i]) for i in range(n)]
    # reversed_rotations = [tuple((v2, v1) for (v1,v2,_) in reversed(r)) for r in rotations]
    # print(reversed_rotations)
    # return min(rotations + reversed_rotations)
    return min(rotations)

def flipped_canonical_face(face):
    n = len(face)
    face_nolabel_flipped = [(e[1],e[0]) for e in face]
    rotations = [tuple(face_nolabel_flipped[i:] + face_nolabel_flipped[:i]) for i in range(n)]
    return min(rotations)

def canonical_image(face,map):
    image = canonical_face([
        (map(edge[0]),map(edge[1]))
        for edge in face
    ])
    return image

def is_orientation_preserving_isomorphic(DFV1:DFV, DFV2:DFV):
    if not DFV1.V.is_isomorphic(DFV2.V):
        return False, None
    
    # 检查 Delaunay 图是否同构
    is_iso, iso_map = DFV1.D.is_isomorphic(DFV2.D, certificate=True)

    # 一个特判: 把0维退化胞腔丢掉.
    if DFV1.poly.is_empty():
        return True, iso_map

    if not is_iso:
        # print("not iso") 
        return False, None
    else:
        canonical_f2 = [canonical_face(f2) for f2 in DFV2.F]
        
        for phi in DFV1.D.automorphism_group():
            # 构造顶点置换dict
            current_map_dict = {v: iso_map[phi(v)] for v in DFV1.D.vertices()}
            
            mapping_func = lambda a: current_map_dict[a]
            
            is_iso = True
            for f1 in DFV1.F: # 检查所有面是否保持定向
                if all(f2 != canonical_image(f1, mapping_func) for f2 in canonical_f2):
                    is_iso = False
                    break
            
            # 如果找到了保持定向的同构，返回 True 和该顶点映射字典
            if is_iso: 
                return True, current_map_dict
        
        return False, None

def cyclic_match_position(source_word, target_word):
    """
    找 source_word 与 target_word 的有向循环匹配位移.
    若不存在, 返回 None.

    返回 shift, 使得:
        source_word[k] == target_word[(k + shift) % n]
    """
    n = len(source_word)

    if n != len(target_word):
        return None

    for shift in range(n):
        if all(
            source_word[k] == target_word[(k + shift) % n]
            for k in range(n)
        ):
            return shift

    return None

def induced_edge_label_perm_from_faces(dfv:DFV,g):
    """
    给定 DFV 和置换 g 求解边置换的 dict. 
    若 perm 不保持面集 F, 返回None.

    返回 label_perm, 使得
        label_perm[old_label] == new_label.
    """
    target_words = [
        [(u, v) for u, v, _ in face]
        for face in dfv.F
    ]

    label_perm = {}

    for source_face in dfv.F:
        mapped_word = [
            (g(u), g(v))
            for u, v, _ in source_face
        ]

        match = None

        for j, target_word in enumerate(target_words):
            shift = cyclic_match_position(mapped_word, target_word)

            if shift is not None:
                match = (j, shift)
                break

        if match is None:
            return None
            
        j, shift = match
        target_face = dfv.F[j]
        n = len(source_face)

        for k, (_, _, source_label) in enumerate(source_face):
            target_label = target_face[(k + shift) % n][2]

            if source_label in label_perm:
                if label_perm[source_label] != target_label:
                    return None
            else:
                label_perm[source_label] = target_label

    return label_perm

def orbits_of_perm(label_perm, labels):
    """
    label_perm: dict, 表示单个 perm 的作用: label -> label
    labels: 所有 label 列表

    返回: 所有轨道
    """
    seen = []
    orbits = []

    for x in labels:
        if x in seen:
            continue

        orbit = []
        cur = x

        while cur not in orbit:
            orbit.append(cur)
            cur = label_perm[cur]

        orbits.append(orbit)
        seen += orbit

    return orbits

def delaunay_edge_orbit_of_g(dfv:DFV, g):
    """
    给定 DFV 和顶点置换 g 求解边置换的轨道.

    返回: 所有边的轨道
    """
    labels = dfv.edge_labels
    label_perm = induced_edge_label_perm_from_faces(dfv,g)

    return orbits_of_perm(label_perm,labels)

def constraints_from_orbits(orbits,labels):
    """
    从 orbits 求解 Polyhedron 构造函数所需的线性条件.
    """
    n = len(labels) # amount of parameters
    eq_constraints = [[0]*(n+1)]

    for orbit in orbits:
        for a, b in combinations(orbit, 2):
            eq_constraints.append(
                [0]+[1 if a == label else -1 if b == label else 0 for label in labels]
            )
    
    return eq_constraints