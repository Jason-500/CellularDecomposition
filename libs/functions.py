from sage.all import *
from sage.geometry.polyhedron.constructor import Polyhedron

def is_edge_on_face(edge,face):
    """
    Determine if an edge is on a given face.

    Modified from edge_on_faces()
    """
    return any([e[2]==edge[2] for e in face])

def is_edge_on_faces(edge,faces):
    """
    Determine if an edge is contained by a bunch of face.

    Modified from edge_on_faces()
    """
    return any([is_edge_on_face(edge,f) for f in faces])

def is_loop(edge):
    """
    Determine if an edge is a loop.
    """
    return edge[0]==edge[1]

def generate_voronoi(D,F):
    """
    Generate Voronoi graph (with label) with Delaunay and corresponding faces. the number representing vertex is the corresponding index of face.

        D: Delaunay graph;
        F: Faces of Delaunay;

    Modified from generate_voronoi()
    """
    V = []
    for e in D.edges():
        indices = [i for i in range(len(F)) if is_edge_on_face(e,F[i])]
        if len(indices)==1:
            V.append((indices[0],indices[0],e[2]))
        else:
            V.append((indices[0],indices[1],e[2]))
    
    return Graph(V, multiedges=True, loops=True)


def poly(D,F):
    """
    Calculate the dimension of polygon representing the cell labelled by given Delaunay graph with face information.

        D: Delaunay graph;
        
        F: Faces of Delaunay;

    Modified from poly_dimension()
    """

    n = len(D.edges()) # amount of parameters

    ieq_constraints = []
    eq_constraints = []

    # constraints that all \theta_e are between 0 and 1 for any edge e. (scaled from [0,\pi]) 24 constraints.
    for i in range(n):
        lb = [1 if j==i+1 else 0 for j in range(n+1)] # lower bound
        ub = [1 if j==0 else -1 if j==i+1 else 0 for j in range(n+1)] # upper bound
        ieq_constraints += [tuple(lb),tuple(ub)]

    # constraints that came from Rivin's Theorem: -|F'|+\sum \theta_e' >=0. 255 constraints.
    psf = list(powerset(range(len(F))))
    for s in psf: # generate constraints for all subset; Idea 1 applys here.
        if s:
            sf = [F[si] for si in s]
            tc = [-len(s)] + [1 if is_edge_on_faces(e,sf) else 0 for e in D.edges()] 
            ieq_constraints += [tuple(tc)]

    # constraints that cone angles at vertices of Delaunay are all pi. equations. 6 constraints.
    for v in D.vertices():
        vc = [D.degree(v)-1] + [(-2 if is_loop(e) else -1) if v in e else 0 for e in D.edges()]
        eq_constraints += [tuple(vc)]
    
    return Polyhedron(ieqs=ieq_constraints,eqns=eq_constraints,base_ring=QQ)


# methods for modification
def edge_removed_face(edge,face):
    """
    Return face removed given edge. If no such edge exist, return original face.

    Modified from face_modification()
    """
    if not is_edge_on_face(edge,face):
        return face
    else:
        return [e for e in face if e[2]!=edge[2]]

def merge_faces(edge,faces):
    """
    Merge two faces along single edge.

    Modified from face_modification()
    """
    if len(faces)>2 or not faces:
        return ('error')
    elif len(faces)==1:
        return edge_removed_face(edge,faces[0]) # This case will generate disconnected graph, which will be ruled out.
    else:
        face1, face2 = tuple(faces)
        face = []
        for i, e in enumerate(face1):
            if e[2] == edge[2]:
                face += face1[i+1:] + face1[:i]
                break
        for i, e in enumerate(face2):
            if e[2] == edge[2]:
                face += face2[i+1:] + face2[:i]
                break
        return face
        

def modified_faces_with_voronoi(edge,F,V):
    """
    Return faces removed given edge and the corresponding Voronoi graph. The faces separated by this edge will be fused together. 

    Note: the edge here is from Delaunay graph, only useful information is its label.

    Modified from face_modification()
    """
    for e in V.edges():
        if e[2]==edge[2]: 
            fusing_vertices = set(e[:2])
            v=V.copy()
            v.delete_edge(e)
            break
    
    nf = len(F)
    faces_to_merge = [F[vi] for vi in fusing_vertices]
    f = [F[i] for i in range(nf) if i not in fusing_vertices]
    f.append(merge_faces(edge,faces_to_merge)) # Adjust the order of edges from fused face. This fused face is at len(F)-len(f..._v...) position.

    # construct mapping of vertices for contracting Voronoi.
    mapping_v = [
        nf - len(fusing_vertices) if i in fusing_vertices
        else i - sum(1 for j in fusing_vertices if j < i)
        for i in range(nf)
    ]

    v = Graph([(mapping_v[e[0]],mapping_v[e[1]],e[2]) for e in v.edges()], multiedges=True, loops=True)
    return (f,v)

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

def is_orientation_preserving_isomorphic(DFV1, DFV2):
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

    