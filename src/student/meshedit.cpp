#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <cmath>

#include "../geometry/halfedge.h"
#include "debug.h"

namespace {
using H = Halfedge_Mesh;

H::HalfedgeRef prev_in_face(H::HalfedgeRef h) {
    auto cur = h;
    while(cur->next() != h) cur = cur->next();
    return cur;
}

Vec3 face_centroid(H::FaceRef f) {
    Vec3 c;
    float n = 0.0f;
    auto h = f->halfedge();
    do {
        c += h->vertex()->pos;
        n += 1.0f;
        h = h->next();
    } while(h != f->halfedge());
    return (n > 0.0f) ? (c / n) : Vec3{};
}

std::vector<H::VertexRef> face_vertices(H::FaceRef f) {
    std::vector<H::VertexRef> vs;
    auto h = f->halfedge();
    do {
        vs.push_back(h->vertex());
        h = h->next();
    } while(h != f->halfedge());
    return vs;
}

bool is_triangle_face(H::FaceRef f) {
    return !f->is_boundary() && f->degree() == 3;
}

static unsigned int hash_u32(unsigned int x) {
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

static float rand01(unsigned int x) {
    return (hash_u32(x) & 0x00ffffffu) / float(0x01000000u);
}

std::vector<H::VertexRef> vertices_in_face_rings(H::FaceRef f, unsigned int rings) {

    std::vector<H::VertexRef> out;
    std::unordered_set<unsigned int> visited;
    std::queue<std::pair<H::VertexRef, unsigned int>> q;

    auto h = f->halfedge();
    do {
        auto v = h->vertex();
        if(visited.insert(v->id()).second) {
            out.push_back(v);
            q.push({v, 0u});
        }
        h = h->next();
    } while(h != f->halfedge());

    while(!q.empty()) {
        auto [v, d] = q.front();
        q.pop();
        if(d >= rings) continue;

        auto hv = v->halfedge();
        auto cur = hv;
        do {
            auto vn = cur->twin()->vertex();
            if(visited.insert(vn->id()).second) {
                out.push_back(vn);
                q.push({vn, d + 1u});
            }
            cur = cur->twin()->next();
        } while(cur != hv);
    }

    return out;
}

} // namespace


/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it does not want to perform the operation for
    whatever reason (e.g. you don't want to allow the user to erase the last vertex).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementaiton, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
    This method should replace the given vertex and all its neighboring
    edges and faces with a single face, returning the new face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_vertex(Halfedge_Mesh::VertexRef v) {

    (void)v;
    return std::nullopt;
}

/*
    This method should erase the given edge and return an iterator to the
    merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_edge(Halfedge_Mesh::EdgeRef e) {

    (void)e;
    return std::nullopt;
}

/*
    This method should collapse the given edge and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {

if(e->on_boundary()) return std::nullopt;

    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h3 = h0->twin();
    FaceRef f0 = h0->face();
    FaceRef f1 = h3->face();
    if(f0->is_boundary() || f1->is_boundary()) return std::nullopt;
    if(f0->degree() != 3 || f1->degree() != 3) return std::nullopt;

    // Neighborhood halfedges (two triangles).
    HalfedgeRef h1 = h0->next();
    HalfedgeRef h2 = h1->next();
    HalfedgeRef h4 = h3->next();
    HalfedgeRef h5 = h4->next();

    // Outside halfedges across the two "other" edges adjacent to v1.
    HalfedgeRef h6 = h1->twin(); // v2 -> v1
    HalfedgeRef h5t = h5->twin(); // v1 -> v3

    // Outside halfedges across the two "other" edges adjacent to v0.
    HalfedgeRef h2t = h2->twin(); // v0 -> v2
    HalfedgeRef h4t = h4->twin(); // v3 -> v0

    VertexRef v0 = h0->vertex();
    VertexRef v1 = h3->vertex();
    VertexRef v2 = h2->vertex();
    VertexRef v3 = h5->vertex();

    // Move the kept vertex to the midpoint.
    v0->pos = 0.5f * (v0->pos + v1->pos);

    // Re-point all halfedges originating at v1 to originate at v0 (except those that will be erased).
    {
        HalfedgeRef hv = v1->halfedge();
        do {
            if(hv != h3 && hv != h1) {
                hv->vertex() = v0;
            }
            hv = hv->twin()->next();
        } while(hv != v1->halfedge());
    }

    // Stitch edge (v0,v2): use the existing edge of h2, but make twins be (v0->v2) and (v2->v0)
    // from the surrounding mesh, i.e., (h2t) and (h6) after collapsing v1->v0.
    EdgeRef e_v0v2 = h2->edge();
    EdgeRef e_v1v2 = h1->edge();
    h6->vertex() = v2; // ensure start is v2 (already), explicit for clarity
    h6->twin() = h2t;
    h2t->twin() = h6;
    h6->edge() = e_v0v2;
    h2t->edge() = e_v0v2;
    e_v0v2->halfedge() = h2t;

    // Stitch edge (v0,v3): keep edge of h4, twins become (v0->v3) and (v3->v0)
    // using (h5t) (after collapsing v1->v0) and (h4t).
    EdgeRef e_v0v3 = h4->edge();
    EdgeRef e_v1v3 = h5->edge();
    h5t->vertex() = v0; // v1 becomes v0
    h5t->twin() = h4t;
    h4t->twin() = h5t;
    h5t->edge() = e_v0v3;
    h4t->edge() = e_v0v3;
    e_v0v3->halfedge() = h5t;

    // Ensure v0 has a valid outgoing halfedge.
    v0->halfedge() = h2t;
    v2->halfedge() = h6;
    v3->halfedge() = h4t;

    // Erase the two faces, collapsing edge, and the two redundant edges from v1.
    erase(f0);
    erase(f1);
    erase(e);
    erase(e_v1v2);
    erase(e_v1v3);

    // Erase halfedges that lived in the two removed faces.
    erase(h0);
    erase(h1);
    erase(h2);
    erase(h3);
    erase(h4);
    erase(h5);

    // The vertex v1 is removed.
    erase(v1);

    return v0;
}

/*
    This method should collapse the given face and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(Halfedge_Mesh::FaceRef f) {
    if(f->is_boundary()) return std::nullopt;

    std::vector<HalfedgeRef> inner_halfedges;
    std::vector<VertexRef> verts;
    std::vector<EdgeRef> inner_edges;

    HalfedgeRef h = f->halfedge();
    do {
        inner_halfedges.push_back(h);
        verts.push_back(h->vertex());
        inner_edges.push_back(h->edge());
        h = h->next();
    } while(h != f->halfedge());

    Vec3 centroid;
    for(VertexRef v : verts) centroid += v->pos;
    centroid /= (float)verts.size();

    VertexRef vm = new_vertex();
    vm->pos = centroid;

    // for each inner halfedge, redirect its twin's vertex to vm
    for(HalfedgeRef ih : inner_halfedges) {
        HalfedgeRef twin = ih->twin();
        // redirect all halfedges that originated from this vertex to vm
        HalfedgeRef cur = twin;
        do {
            cur->vertex() = vm;
            cur = cur->twin()->next();
        } while(cur != twin);
    }

    vm->halfedge() = inner_halfedges[0]->twin()->next();
    
    // stitch outer halfedges
    for(HalfedgeRef ih : inner_halfedges) {
        HalfedgeRef twin = ih->twin();
        HalfedgeRef prev_twin = twin;
        while(prev_twin->next() != twin) prev_twin = prev_twin->next();
        HalfedgeRef outer_prev = prev_in_face(twin);
        outer_prev->next() = twin->next();
    }
    for(HalfedgeRef ih : inner_halfedges) {
        HalfedgeRef twin = ih->twin();
        if(!twin->face()->is_boundary()) {
            twin->face()->halfedge() = twin->next();
        }
    }

    for(HalfedgeRef ih : inner_halfedges) {
        erase(ih->twin());
        erase(ih);
    }
    for(EdgeRef ie : inner_edges) erase(ie);
    for(VertexRef v : verts) erase(v);
    erase(f);

    return vm;
}


/*
    This method should flip the given edge and return an iterator to the
    flipped edge.
*/
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(Halfedge_Mesh::EdgeRef e) {
    if(e->on_boundary()) return std::nullopt;

    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h3 = h0->twin();
    FaceRef f0 = h0->face();
    FaceRef f1 = h3->face();
    if(f0->is_boundary() || f1->is_boundary()) return std::nullopt;
    if(f0->degree() != 3 || f1->degree() != 3) return std::nullopt;

    // collect all halfedges
    HalfedgeRef h1 = h0->next();
    HalfedgeRef h2 = h1->next();
    HalfedgeRef h4 = h3->next();
    HalfedgeRef h5 = h4->next();

    // collect vertices
    VertexRef v0 = h0->vertex();
    VertexRef v1 = h3->vertex();
    VertexRef v2 = h2->vertex();
    VertexRef v3 = h5->vertex();

    // reassign every pointer exhaustively
    h0->set_neighbors(h2, h3, v3, e, f0);
    h3->set_neighbors(h5, h0, v2, e, f1);
    h1->set_neighbors(h3, h1->twin(), v1, h1->edge(), f1);
    h2->set_neighbors(h4, h2->twin(), v2, h2->edge(), f0);
    h4->set_neighbors(h0, h4->twin(), v0, h4->edge(), f0);
    h5->set_neighbors(h1, h5->twin(), v3, h5->edge(), f1);

    // fix vertex halfedge pointers
    v0->halfedge() = h4;
    v1->halfedge() = h1;
    v2->halfedge() = h3;
    v3->halfedge() = h0;

    // fix face halfedge pointers
    f0->halfedge() = h0;
    f1->halfedge() = h3;

    return e;
}

/*
    This method should split the given edge and return an iterator to the
    newly inserted vertex. The halfedge of this vertex should point along
    the edge that was split, rather than the new edges.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(Halfedge_Mesh::EdgeRef e) {

    if(e->on_boundary()) return std::nullopt;

    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h3 = h0->twin();
    FaceRef f0 = h0->face();
    FaceRef f1 = h3->face();
    if(!is_triangle_face(f0) || !is_triangle_face(f1)) return std::nullopt;

    HalfedgeRef h1 = h0->next();
    HalfedgeRef h2 = h1->next();
    HalfedgeRef h4 = h3->next();
    HalfedgeRef h5 = h4->next();

    VertexRef v0 = h0->vertex();
    VertexRef v1 = h3->vertex();
    VertexRef v2 = h2->vertex();
    VertexRef v3 = h5->vertex();

    // Allocate new elements.
    VertexRef vm = new_vertex();
    vm->pos = 0.5f * (v0->pos + v1->pos);

    FaceRef f2 = new_face(false);
    FaceRef f3 = new_face(false);

    EdgeRef e_vm_v1 = new_edge();
    EdgeRef e_vm_v2 = new_edge();
    EdgeRef e_vm_v3 = new_edge();

    HalfedgeRef h6 = new_halfedge();  // vm -> v1
    HalfedgeRef h7 = new_halfedge();  // v1 -> vm
    HalfedgeRef h8 = new_halfedge();  // vm -> v2
    HalfedgeRef h9 = new_halfedge();  // v2 -> vm
    HalfedgeRef h10 = new_halfedge(); // vm -> v3
    HalfedgeRef h11 = new_halfedge(); // v3 -> vm

    // Set twins.
    h6->twin() = h7;
    h7->twin() = h6;
    h8->twin() = h9;
    h9->twin() = h8;
    h10->twin() = h11;
    h11->twin() = h10;

    // Reuse the original edge for (v0,vm).
    h0->vertex() = v0; // v0 -> vm
    h3->vertex() = vm; // vm -> v0

    // Assign edges.
    h0->edge() = e;
    h3->edge() = e;
    e->halfedge() = h0;

    h6->edge() = e_vm_v1;
    h7->edge() = e_vm_v1;
    e_vm_v1->halfedge() = h6;

    h8->edge() = e_vm_v2;
    h9->edge() = e_vm_v2;
    e_vm_v2->halfedge() = h8;

    h10->edge() = e_vm_v3;
    h11->edge() = e_vm_v3;
    e_vm_v3->halfedge() = h10;

    // Faces and next pointers.
    // f0 becomes (v0, vm, v2): h0 (v0->vm), h8 (vm->v2), h2 (v2->v0)
    h0->face() = f0;
    h8->face() = f0;
    h2->face() = f0;
    h0->next() = h8;
    h8->next() = h2;
    h2->next() = h0;
    f0->halfedge() = h0;

    // f2 is (vm, v1, v2): h6 (vm->v1), h1 (v1->v2), h9 (v2->vm)
    h6->face() = f2;
    h1->face() = f2;
    h9->face() = f2;
    h6->next() = h1;
    h1->next() = h9;
    h9->next() = h6;
    f2->halfedge() = h6;

    // f1 becomes (v1, vm, v3): h7 (v1->vm), h10 (vm->v3), h5 (v3->v1)
    h7->face() = f1;
    h10->face() = f1;
    h5->face() = f1;
    h7->next() = h10;
    h10->next() = h5;
    h5->next() = h7;
    f1->halfedge() = h7;

    // f3 is (vm, v0, v3): h3 (vm->v0), h4 (v0->v3), h11 (v3->vm)
    h3->face() = f3;
    h4->face() = f3;
    h11->face() = f3;
    h3->next() = h4;
    h4->next() = h11;
    h11->next() = h3;
    f3->halfedge() = h3;

    // Set remaining vertices for new halfedges.
    h6->vertex() = vm;
    h7->vertex() = v1;
    h8->vertex() = vm;
    h9->vertex() = v2;
    h10->vertex() = vm;
    h11->vertex() = v3;

    // Update vertex halfedge pointers.
    v0->halfedge() = h0;
    v1->halfedge() = h1;
    v2->halfedge() = h2;
    v3->halfedge() = h5;
    vm->halfedge() = h3; // along the (split) original edge

    return vm;
}


/* Note on the beveling process:

    Each of the bevel_vertex, bevel_edge, and bevel_face functions do not represent
    a full bevel operation. Instead, they should update the _connectivity_ of
    the mesh, _not_ the positions of newly created vertices. In fact, you should set
    the positions of new vertices to be exactly the same as wherever they "started from."

    When you click on a mesh element while in bevel mode, one of those three functions
    is called. But, because you may then adjust the distance/offset of the newly
    beveled face, we need another method of updating the positions of the new vertices.

    This is where bevel_vertex_positions, bevel_edge_positions, and
    bevel_face_positions come in: these functions are called repeatedly as you
    move your mouse, the position of which determins the normal and tangent offset
    parameters. These functions are also passed an array of the original vertex
    positions: for  bevel_vertex, it has one element, the original vertex position,
    for bevel_edge,  two for the two vertices, and for bevel_face, it has the original
    position of each vertex in halfedge order. You should use these positions, as well
    as the normal and tangent offset fields to assign positions to the new vertices.

    Finally, note that the normal and tangent offsets are not relative values - you
    should compute a particular new position from them, not a delta to apply.
*/

/*
    This method should replace the vertex v with a face, corresponding to
    a bevel operation. It should return the new face.  NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_vertex_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(Halfedge_Mesh::VertexRef v) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)v;
    return std::nullopt;
}

/*
    This method should replace the edge e with a face, corresponding to a
    bevel operation. It should return the new face. NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_edge_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(Halfedge_Mesh::EdgeRef e) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)e;
    return std::nullopt;
}

/*
    This method should replace the face f with an additional, inset face
    (and ring of faces around it), corresponding to a bevel operation. It
    should return the new face.  NOTE: This method is responsible for updating
    the *connectivity* of the mesh only---it does not need to update the vertex
    positions. These positions will be updated in
    Halfedge_Mesh::bevel_face_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_face(Halfedge_Mesh::FaceRef f) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."
    if(f->is_boundary()) return std::nullopt;

    // Gather the boundary halfedges/vertices of the face.
    std::vector<HalfedgeRef> hs;
    std::vector<VertexRef> vs;
    {
        HalfedgeRef h = f->halfedge();
        do {
            hs.push_back(h);
            vs.push_back(h->vertex());
            h = h->next();
        } while(h != f->halfedge());
    }
    const size_t n = hs.size();
    if(n < 3) return std::nullopt;

    // Reuse the original face as the inset face.
    FaceRef inset = f;

    // Allocate new side faces, inset vertices, and new edges/halfedges.
    std::vector<FaceRef> side(n);
    std::vector<VertexRef> u(n);
    std::vector<EdgeRef> e_inset(n);
    std::vector<EdgeRef> e_spoke(n);
    std::vector<HalfedgeRef> h_inset(n), h_inset_twin(n);
    std::vector<HalfedgeRef> h_spoke_out(n), h_spoke_in(n);

    for(size_t i = 0; i < n; i++) {
        side[i] = new_face(false);
        u[i] = new_vertex();
        u[i]->pos = vs[i]->pos;
        e_inset[i] = new_edge();
        e_spoke[i] = new_edge();

        h_inset[i] = new_halfedge();
        h_inset_twin[i] = new_halfedge();
        h_spoke_out[i] = new_halfedge();
        h_spoke_in[i] = new_halfedge();
    }

    // Set up inset face loop (u_i -> u_{i+1}).
    for(size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        h_inset[i]->set_neighbors(h_inset[j], h_inset_twin[i], u[i], e_inset[i], inset);
        h_inset_twin[i]->set_neighbors(h_spoke_out[i], h_inset[i], u[j], e_inset[i], side[i]);
        e_inset[i]->halfedge() = h_inset[i];
    }
    inset->halfedge() = h_inset[0];

    // Set up spoke halfedges and their twins across adjacent side faces.
    for(size_t i = 0; i < n; i++) {
        size_t prev = (i + n - 1) % n;
        // u_i -> v_i belongs to side_i
        h_spoke_out[i]->twin() = h_spoke_in[i];
        // v_i -> u_i belongs to side_{i-1}
        h_spoke_in[i]->twin() = h_spoke_out[i];

        h_spoke_out[i]->vertex() = u[i];
        h_spoke_in[i]->vertex() = vs[i];

        h_spoke_out[i]->edge() = e_spoke[i];
        h_spoke_in[i]->edge() = e_spoke[i];
        e_spoke[i]->halfedge() = h_spoke_out[i];

        // Assign faces later when wiring side face cycles.
        (void)prev;
    }

    // Wire each side face: (v_i -> v_{i+1}) [existing], (v_{i+1} -> u_{i+1}), (u_{i+1} -> u_i),
    // (u_i -> v_i).
    for(size_t i = 0; i < n; i++) {
        size_t next = (i + 1) % n;

        HalfedgeRef h_orig = hs[i]; // v_i -> v_{i+1}
        h_orig->face() = side[i];

        // v_{i+1} -> u_{i+1} is the "in" spoke halfedge for vertex (i+1), which belongs to side_i
        HalfedgeRef h_vnext_unext = h_spoke_in[next];
        // u_{i+1} -> u_i is the twin inset halfedge for edge i
        HalfedgeRef h_unext_ui = h_inset_twin[i];
        // u_i -> v_i is the "out" spoke halfedge for vertex i
        HalfedgeRef h_ui_vi = h_spoke_out[i];

        // Set faces.
        h_vnext_unext->face() = side[i];
        h_unext_ui->face() = side[i];
        h_ui_vi->face() = side[i];

        // Set next pointers.
        h_orig->next() = h_vnext_unext;
        h_vnext_unext->next() = h_unext_ui;
        h_unext_ui->next() = h_ui_vi;
        h_ui_vi->next() = h_orig;

        // Ensure vertices/edges are consistent (some were set above, but be explicit).
        h_vnext_unext->vertex() = vs[next];
        h_unext_ui->vertex() = u[next];
        h_ui_vi->vertex() = u[i];

        side[i]->halfedge() = h_orig;
    }

    // Ensure each inset vertex has an outgoing halfedge.
    for(size_t i = 0; i < n; i++) {
        u[i]->halfedge() = h_inset[i];
    }

    return inset;
}


/*
    Compute new vertex positions for the vertices of the beveled vertex.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the original vertex position and its associated outgoing edge
    to compute a new vertex position along the outgoing edge.
*/
void Halfedge_Mesh::bevel_vertex_positions(const std::vector<Vec3>& start_positions,
                                           Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled edge.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the orig array) to compute an offset vertex position.

    Note that there is a 1-to-1 correspondence between halfedges in
    newHalfedges and vertex positions
    in orig.  So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vector3D pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_edge_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled face.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the start_positions array) to compute an offset vertex
    position.

    Note that there is a 1-to-1 correspondence between halfedges in
    new_halfedges and vertex positions
    in orig. So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vec3 pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_face_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset,
                                         float normal_offset) {

    if(flip_orientation) normal_offset = -normal_offset;
    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
    (void)normal_offset;
}

/*
    Splits all non-triangular faces into triangles.
*/
void Halfedge_Mesh::triangulate() {

    std::unordered_map<VertexRef, Index> v_index;
    std::vector<Vec3> verts;
    verts.reserve(n_vertices());

    Index next = 0;
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        v_index[v] = next++;
        verts.push_back(v->pos);
    }

    std::vector<std::vector<Index>> tris;
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        if(f->is_boundary()) continue;
        auto vs = face_vertices(f);
        if(vs.size() < 3) continue;
        if(vs.size() == 3) {
            tris.push_back({v_index[vs[0]], v_index[vs[1]], v_index[vs[2]]});
        } else {
            Index i0 = v_index[vs[0]];
            for(size_t i = 1; i + 1 < vs.size(); i++) {
                tris.push_back({i0, v_index[vs[i]], v_index[vs[i + 1]]});
            }
        }
    }

    from_poly(tris, verts);
}


/* Note on the quad subdivision process:

        Unlike the local mesh operations (like bevel or edge flip), we will perform
        subdivision by splitting *all* faces into quads "simultaneously."  Rather
        than operating directly on the halfedge data structure (which as you've
        seen is quite difficult to maintain!) we are going to do something a bit nicer:
           1. Create a raw list of vertex positions and faces (rather than a full-
              blown halfedge mesh).
           2. Build a new halfedge mesh from these lists, replacing the old one.
        Sometimes rebuilding a data structure from scratch is simpler (and even
        more efficient) than incrementally modifying the existing one.  These steps are
        detailed below.

  Step I: Compute the vertex positions for the subdivided mesh.
        Here we're going to do something a little bit strange: since we will
        have one vertex in the subdivided mesh for each vertex, edge, and face in
        the original mesh, we can nicely store the new vertex *positions* as
        attributes on vertices, edges, and faces of the original mesh. These positions
        can then be conveniently copied into the new, subdivided mesh.
        This is what you will implement in linear_subdivide_positions() and
        catmullclark_subdivide_positions().

  Steps II-IV are provided (see Halfedge_Mesh::subdivide()), but are still detailed
  here:

  Step II: Assign a unique index (starting at 0) to each vertex, edge, and
        face in the original mesh. These indices will be the indices of the
        vertices in the new (subdivided mesh).  They do not have to be assigned
        in any particular order, so long as no index is shared by more than one
        mesh element, and the total number of indices is equal to V+E+F, i.e.,
        the total number of vertices plus edges plus faces in the original mesh.
        Basically we just need a one-to-one mapping between original mesh elements
        and subdivided mesh vertices.

  Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
        the element indices defined above. In other words, each new quad should be
        of the form (i,j,k,l), where i,j,k and l are four of the indices stored on
        our original mesh elements.  Note that it is essential to get the orientation
        right here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces
        should circulate in the same direction as old faces (think about the right-hand
        rule).

  Step IV: Pass the list of vertices and quads to a routine that clears
        the internal data for this halfedge mesh, and builds new halfedge data from
        scratch, using the two lists.
*/

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    simple linear interpolation, e.g., the edge midpoints and face
    centroids.
*/
void Halfedge_Mesh::linear_subdivide_positions() {

    // For each vertex, assign Vertex::new_pos to
    // its original position, Vertex::pos.

    // For each edge, assign the midpoint of the two original
    // positions to Edge::new_pos.

    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::new_pos. Note
    // that in general, NOT all faces will be triangles!
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        v->new_pos = v->pos;
    }

    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        HalfedgeRef h = e->halfedge();
        Vec3 p0 = h->vertex()->pos;
        Vec3 p1 = h->twin()->vertex()->pos;
        e->new_pos = 0.5f * (p0 + p1);
    }

    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        if(f->is_boundary()) continue;
        f->new_pos = face_centroid(f);
    }
    
}


/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    the Catmull-Clark rules for subdivision.

    Note: this will only be called on meshes without boundary
*/
void Halfedge_Mesh::catmullclark_subdivide_positions() {
    // The implementation for this routine should be
    // a lot like Halfedge_Mesh:linear_subdivide_positions:(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        if(f->is_boundary()) continue;
        f->new_pos = face_centroid(f);
    }


    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        HalfedgeRef h = e->halfedge();
        VertexRef v0 = h->vertex();
        VertexRef v1 = h->twin()->vertex();
        FaceRef f0 = h->face();
        FaceRef f1 = h->twin()->face();
        e->new_pos = (v0->pos + v1->pos + f0->new_pos + f1->new_pos) / 4.0f;
    }

    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        // Catmull-Clark vertex rule:
        // new = (F + 2R + (n-3)P) / n
        Vec3 F, R;
        float n = 0.0f;

        HalfedgeRef h = v->halfedge();
        do {
            F += h->face()->new_pos;
            Vec3 p0 = h->vertex()->pos;
            Vec3 p1 = h->twin()->vertex()->pos;
            R += 0.5f * (p0 + p1);
            n += 1.0f;
            h = h->twin()->next();
        } while(h != v->halfedge());

        F /= n;
        R /= n;
        v->new_pos = (F + 2.0f * R + (n - 3.0f) * v->pos) / n;
    }
}

/*
        This routine should increase the number of triangles in the mesh
        using Loop subdivision. Note: this is will only be called on triangle meshes.
*/
void Halfedge_Mesh::loop_subdivide() {

    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::new_pos.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::new_pos.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::is_new. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::pos.

    // Each vertex and edge of the original surface can be associated with a
    // vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions
    // using the connectivity of the original (coarse) mesh; navigating this mesh
    // will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We
    // will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using
    // the Loop subdivision rule.

    // Next, compute the updated vertex positions associated with edges.

    // Next, we're going to split every edge in the mesh, in any order. For
    // future reference, we're also going to store some information about which
    // subdivided edges come from splitting an edge in the original mesh, and
    // which edges are new.
    // In this loop, we only want to iterate over edges of the original
    // mesh---otherwise, we'll end up splitting edges that we just split (and
    // the loop will never end!)

    // Finally, flip any new edge that connects an old and new vertex.

    // Copy the updated vertex positions to the subdivided mesh.
}

/*
    Isotropic remeshing. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if this is not a triangle mesh)
*/
bool Halfedge_Mesh::isotropic_remesh() {

    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}

/* Helper type for quadric simplification */
struct Edge_Record {
    Edge_Record() {
    }
    Edge_Record(std::unordered_map<Halfedge_Mesh::VertexRef, Mat4>& vertex_quadrics,
                Halfedge_Mesh::EdgeRef e)
        : edge(e) {
        Halfedge_Mesh::HalfedgeRef h = e->halfedge();
        Halfedge_Mesh::VertexRef v0 = h->vertex();
        Halfedge_Mesh::VertexRef v1 = h->twin()->vertex();
        Mat4 K = vertex_quadrics[v0] + vertex_quadrics[v1];

        // Build 4x4 M: first 3 rows from K, 4th row (0,0,0,1). Solve M*(x,y,z,1)^T = (0,0,0,1)^T.
        Mat4 M;
        for(int j = 0; j < 4; j++) {
            M[j] = Vec4(K[j][0], K[j][1], K[j][2], (j == 3) ? 1.0f : 0.0f);
        }
        Mat4 M_inv = Mat4::inverse(M);
        optimal = Vec3(M_inv[3][0], M_inv[3][1], M_inv[3][2]);

        // If singular (e.g. coplanar), use edge midpoint.
        if(optimal.x != optimal.x || optimal.y != optimal.y || optimal.z != optimal.z) {
            optimal = 0.5f * (v0->pos + v1->pos);
        }
        Vec4 u(optimal.x, optimal.y, optimal.z, 1.0f);
        cost = dot(u, K * u);
    }
    Halfedge_Mesh::EdgeRef edge;
    Vec3 optimal;
    float cost;
};

/* Comparison operator for Edge_Records so std::set will properly order them */
bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
    if(r1.cost != r2.cost) {
        return r1.cost < r2.cost;
    }
    Halfedge_Mesh::EdgeRef e1 = r1.edge;
    Halfedge_Mesh::EdgeRef e2 = r2.edge;
    return &*e1 < &*e2;
}

/** Helper type for quadric simplification
 *
 * A PQueue is a minimum-priority queue that
 * allows elements to be both inserted and removed from the
 * queue.  Together, one can easily change the priority of
 * an item by removing it, and re-inserting the same item
 * but with a different priority.  A priority queue, for
 * those who don't remember or haven't seen it before, is a
 * data structure that always keeps track of the item with
 * the smallest priority or "score," even as new elements
 * are inserted and removed.  Priority queues are often an
 * essential component of greedy algorithms, where one wants
 * to iteratively operate on the current "best" element.
 *
 * PQueue is templated on the type T of the object
 * being queued.  For this reason, T must define a comparison
 * operator of the form
 *
 *    bool operator<( const T& t1, const T& t2 )
 *
 * which returns true if and only if t1 is considered to have a
 * lower priority than t2.
 *
 * Basic use of a PQueue might look
 * something like this:
 *
 *    // initialize an empty queue
 *    PQueue<myItemType> queue;
 *
 *    // add some items (which we assume have been created
 *    // elsewhere, each of which has its priority stored as
 *    // some kind of internal member variable)
 *    queue.insert( item1 );
 *    queue.insert( item2 );
 *    queue.insert( item3 );
 *
 *    // get the highest priority item currently in the queue
 *    myItemType highestPriorityItem = queue.top();
 *
 *    // remove the highest priority item, automatically
 *    // promoting the next-highest priority item to the top
 *    queue.pop();
 *
 *    myItemType nextHighestPriorityItem = queue.top();
 *
 *    // Etc.
 *
 *    // We can also remove an item, making sure it is no
 *    // longer in the queue (note that this item may already
 *    // have been removed, if it was the 1st or 2nd-highest
 *    // priority item!)
 *    queue.remove( item2 );
 *
 */
template<class T> struct PQueue {
    void insert(const T& item) {
        queue.insert(item);
    }
    void remove(const T& item) {
        if(queue.find(item) != queue.end()) {
            queue.erase(item);
        }
    }
    const T& top(void) const {
        return *(queue.begin());
    }
    void pop(void) {
        queue.erase(queue.begin());
    }
    size_t size() {
        return queue.size();
    }

    std::set<T> queue;
};

/*
    Mesh simplification. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if you can't simplify the mesh any
    further without destroying it.)
*/
bool Halfedge_Mesh::simplify() {

    std::unordered_map<VertexRef, Mat4> vertex_quadrics;
    std::unordered_map<FaceRef, Mat4> face_quadrics;
    std::unordered_map<EdgeRef, Edge_Record> edge_records;
    PQueue<Edge_Record> edge_queue;

    // 1. count non-boundary faces
    size_t n_tris = 0;
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        if(!f->is_boundary()) n_tris++;
    }
    size_t target_tris = (n_tris / 4 >= 1) ? (n_tris / 4) : 1;

    // 2. face quadrics
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        if(f->is_boundary()) continue;
        Vec3 n = f->normal();
        Vec3 p = f->halfedge()->vertex()->pos;
        float d = -dot(n, p);
        Vec4 v(n.x, n.y, n.z, d);
        face_quadrics[f] = outer(v, v);
    }

    // 3. vertex quadrics
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        Mat4 Q = Mat4::Zero;
        HalfedgeRef h = v->halfedge();
        do {
            FaceRef f = h->face();
            if(!f->is_boundary()) Q += face_quadrics[f];
            h = h->twin()->next();
        } while(h != v->halfedge());
        vertex_quadrics[v] = Q;
    }

    // 4. build edge records and priority queue
    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        Edge_Record rec(vertex_quadrics, e);
        edge_records[e] = rec;
        edge_queue.insert(rec);
    }

    while(n_tris > target_tris && edge_queue.size() > 0) {
        Edge_Record best = edge_queue.top();
        edge_queue.pop();
        EdgeRef e = best.edge;
        if(e->on_boundary()) continue;

        HalfedgeRef h = e->halfedge();
        VertexRef v0 = h->vertex();
        VertexRef v1 = h->twin()->vertex();
        

        // 4. collect edges touching either endpoint
        std::set<EdgeRef> touch;
        HalfedgeRef hv = v0->halfedge();
        do {
            touch.insert(hv->edge());
            hv = hv->twin()->next();
        } while(hv != v0->halfedge());
        hv = v1->halfedge();
        do {
            touch.insert(hv->edge());
            hv = hv->twin()->next();
        } while(hv != v1->halfedge());

        for(EdgeRef e_touch : touch) {
            auto it = edge_records.find(e_touch);
            if(it != edge_records.end()) edge_queue.remove(it->second);
        }

        // new quadric = sum of endpoints.
        Mat4 Q_new = vertex_quadrics[v0] + vertex_quadrics[v1];
        std::optional<VertexRef> opt_v = collapse_edge_erase(e);
        if(!opt_v) continue;

        VertexRef v_new = *opt_v;
        v_new->pos = best.optimal;
        vertex_quadrics[v_new] = Q_new;

        for(EdgeRef e_touch : touch) {
            edge_records.erase(e_touch);
        }

        hv = v_new->halfedge();
        do {
            EdgeRef e_inc = hv->edge();
            if(edge_records.find(e_inc) == edge_records.end()) {
                Edge_Record rec(vertex_quadrics, e_inc);
                edge_records[e_inc] = rec;
                edge_queue.insert(rec);
            }
            hv = hv->twin()->next();
        } while(hv != v_new->halfedge());

        n_tris -= 2; // each collapse removes 2 triangles
    }

    return true;
}


// Extensions: Spikes and Dimples: 
bool Halfedge_Mesh::spikes(float height_scale, unsigned int seed, float power) {

    if(height_scale <= 0.0f) return false;
    if(!(power > 0.0f)) return false;

    double sum = 0.0;
    unsigned int n = 0u;
    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        auto h = e->halfedge();
        Vec3 p0 = h->vertex()->pos;
        Vec3 p1 = h->twin()->vertex()->pos;
        sum += (p1 - p0).norm();
        n++;
    }
    if(n == 0u) return false;

    float mean_edge = float(sum / double(n));
    if(!(mean_edge > 0.0f)) return false;

    float height = height_scale * mean_edge;

    bool moved = false;
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        if(v->on_boundary()) continue;

        Vec3 nn = v->normal();
        if(nn.norm_squared() <= 0.0f) continue;
        nn = nn.unit();

        float r = rand01(seed ^ (v->id() * 0x9e3779b9u));
        float w = 1.0f;

        if(std::abs(power - 1.0f) < 1e-6f) {
            w = r;
        } else if(std::abs(power - 2.0f) < 1e-6f) {
            w = r * r;
        } else if(std::abs(power - 3.0f) < 1e-6f) {
            w = r * r * r;
        } else if(std::abs(power - 4.0f) < 1e-6f) {
            float r2 = r * r;
            w = r2 * r2;
        } else {
            w = std::pow(r, power);
        }

        v->pos += (height * w) * nn;
        moved = true;
    }

    return moved;
}

bool Halfedge_Mesh::dimples(float density, float radius_scale, float depth_scale, unsigned int seed,
                            unsigned int rings) {

    if(density <= 0.0f || radius_scale <= 0.0f || depth_scale <= 0.0f) return false;

    double sum = 0.0;
    unsigned int n = 0u;
    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        auto h = e->halfedge();
        Vec3 p0 = h->vertex()->pos;
        Vec3 p1 = h->twin()->vertex()->pos;
        sum += (p1 - p0).norm();
        n++;
    }
    if(n == 0u) return false;

    float mean_edge = float(sum / double(n));
    if(!(mean_edge > 0.0f)) return false;

    float base_radius = radius_scale * mean_edge;
    float base_depth = depth_scale * mean_edge;

    std::unordered_map<unsigned int, Vec3> pos0;
    std::unordered_map<unsigned int, Vec3> n0;
    std::unordered_map<unsigned int, Vec3> delta;
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        unsigned int id = v->id();
        pos0.emplace(id, v->pos);
        Vec3 nn = v->normal();
        if(nn.norm_squared() > 0.0f) nn = nn.unit();
        n0.emplace(id, nn);
        delta.emplace(id, Vec3{});
    }

    bool stamped = false;
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        if(f->is_boundary()) continue;

        unsigned int h = hash_u32(seed ^ (f->id() * 0x9e3779b9u));
        if(rand01(h) >= density) continue;

        float radius = base_radius * (0.75f + 0.75f * rand01(h ^ 0x243f6a88u));
        float depth = base_depth * (0.60f + 0.80f * rand01(h ^ 0xb7e15162u));
        Vec3 center = f->center();

        auto vs = vertices_in_face_rings(f, rings);
        for(auto v : vs) {
            unsigned int id = v->id();
            Vec3 p = pos0[id];
            float d = (p - center).norm();
            if(d >= radius) continue;

            float x = d / radius;
            float w = 1.0f - x * x;
            w = w * w;

            delta[id] += (-depth * w) * n0[id];
        }

        stamped = true;
    }

    if(!stamped) return false;

    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        unsigned int id = v->id();
        v->pos = pos0[id] + delta[id];
    }

    return true;
}