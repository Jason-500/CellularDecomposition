from sage.all import *
from sage.geometry.polyhedron.constructor import Polyhedron
# from itertools import combinations

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

def polyhedra_3d_html(poly_list, filename="polyhedra.html",
                      opacity=0.45, frame=True,
                      online=True, eps=1e-10,
                      label_vertices=True,
                      label_full_coords=False,
                      label_fontsize=10,
                      show_if_possible=False):
    """
    Export multiple high-dimensional but jointly 3D Sage Polyhedra
    to one interactive threejs HTML file, using one shared isometric
    coordinate system inherited from the original ambient Euclidean space.
    """
    if len(poly_list) == 0:
        raise ValueError("poly_list is empty.")

    ambient_dims = [P.ambient_dim() for P in poly_list]
    if len(set(ambient_dims)) != 1:
        raise ValueError("All polyhedra must have the same ambient dimension.")

    # Collect V-representation data.
    all_verts = []
    all_dirs = []

    per_poly_data = []

    for P in poly_list:
        verts = [vector(RDF, v) for v in P.vertex_generator()]
        rays  = [vector(RDF, r) for r in P.ray_generator()]
        lines = [vector(RDF, l) for l in P.line_generator()]

        if len(verts) == 0:
            raise ValueError("Each polyhedron needs at least one vertex.")

        all_verts += verts
        all_dirs += rays + lines

        per_poly_data.append((P, verts, rays, lines))

    if len(all_verts) == 0:
        raise ValueError("No vertices found.")

    # Shared affine base point.
    p0 = all_verts[0]

    # Directions spanning the joint affine hull.
    dirs = [v - p0 for v in all_verts[1:]] + all_dirs

    # Shared Gram-Schmidt orthonormal basis.
    E = []
    for b in dirs:
        u = vector(RDF, b)
        for e in E:
            u -= (u * e) * e
        n = u.norm()
        if n > eps:
            E.append(u / n)
        if len(E) == 3:
            break

    if len(E) != 3:
        raise ValueError(
            "The joint affine span does not appear to be 3D; found dimension {}.".format(len(E))
        )

    def proj_point(x):
        x = vector(RDF, x) - p0
        return tuple(x * e for e in E)

    def proj_dir(x):
        x = vector(RDF, x)
        return tuple(x * e for e in E)

    # Build one combined Graphics3d object.
    G = None
    projected_polyhedra = []

    from sage.plot.plot3d.shapes2 import text3d

    for k, (P, verts, rays, lines) in enumerate(per_poly_data):
        verts3 = [proj_point(v) for v in verts]
        rays3  = [proj_dir(r) for r in rays]
        lines3 = [proj_dir(l) for l in lines]

        Q = Polyhedron(vertices=verts3, rays=rays3, lines=lines3)
        projected_polyhedra.append(Q)

        colors = [
            "#0072B2",  # blue
            "#D55E00",  # vermillion
            "#009E73",  # green
            "#CC79A7",  # reddish purple
            "#E69F00",  # orange
            "#56B4E9",  # sky blue
            "#F0E442",  # yellow
            "#000000",  # black
        ]

        c = colors[k % len(colors)]

        if Q.dim() == 1:
            H = Q.plot(
                color=c,
                thickness=4,
                frame=frame
            )
        else:
            H = Q.plot(
                opacity=opacity,
                frame=frame,
                color=c
            )

        # # 非 1D 的对象再叠加 wireframe
        # if Q.dim() > 1:
        #     H += Q.plot(
        #         wireframe=True,
        #         color=c,
        #         thickness=3
        #     )

        if label_vertices:
            # Small offset so labels do not sit exactly on the vertices.
            offset = vector(RDF, (0.03, 0.03, 0.03))

            for i, (v_orig, v3) in enumerate(zip(verts, verts3)):
                pos = tuple(vector(RDF, v3) + offset)

                if label_full_coords:
                    label = "P{}v{}={}".format(k, i, tuple(v_orig))
                else:
                    label = "P{}v{}".format(k, i)

                H += text3d(label, pos, fontsize=label_fontsize)

        if G is None:
            G = H
        else:
            G += H

    G.save(filename, viewer="threejs", online=online)

    if show_if_possible:
        try:
            from IPython.display import IFrame, display
            display(IFrame(filename, width=900, height=700))
        except Exception as err:
            print("Saved to {}, but could not display inline: {}".format(filename, err))

    return projected_polyhedra, G, filename, p0, E

def is_subsequence(A, B):
    """
    Check whether list A is a subsequence of list B (order preserved, not necessarily contiguous).
    """
    it = iter(B)
    return all(any(a == b for b in it) for a in A)