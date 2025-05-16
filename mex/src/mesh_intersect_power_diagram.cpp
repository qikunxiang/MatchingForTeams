#include "power_diagram_intersection.hpp"
#include "mesh_intersect_power_diagram.hpp"

/// @brief Compute the triangular cell complex formed by the intersection of a power diagram and a triangular mesh. 
/// Each non-triangular cell in the resulting intersection is triangulated. The function returns a std::vector containing 
/// the points in the triangulation, a std::vector containing the triangulation indices, a std::vector containing
/// the indices of the power cell each output triangle belongs to, and a std::vector containing the indices in the
/// triangular mesh each output triangle belongs to
/// @param mesh_vertices [input] std::vector containing the vertices in the triangular mesh
/// @param mesh_triangulation [input] std::vector containing the triangulation of the triangular mesh
/// @param circles [input] std::vector containing objects of std::pair which contains the centers and squared radii of circles
/// @param output_vertices [output] std::vector containing the vertices in the output triangulation (there might be irrelevant vertices that do not appear in the triangulation)
/// @param output_triangulation [output] std::vector containing the output triangulation
/// @param power_cell_indices [output] std::vector containing the indices of the power cell each output triangle belongs to
/// @param mesh_indices [output] std::vector containing the indices in the triangular mesh each output triangle belongs to
void compute_mesh_intersect_power_diagram(
    const std::vector<Point_2> & mesh_vertices, 
    const std::vector<tri_indices_t> & mesh_triangulation, \
    const std::vector<std::pair<Point_2, exact_field_t>> & circles, \
    std::vector<Point_2> & output_vertices, 
    std::vector<tri_indices_t> & output_triangulation, \
    std::vector<size_t> & power_cell_indices, 
    std::vector<size_t> & mesh_indices)
{
    // Overall procedure of the function:
    // [1] sweep through the mesh triangles and create Triangle_2 objects and Segment_2 objects
    // [2] compute the segments in the power diagram
    // [3] compute the arrangement
    // [4] sweep through all faces in the arrangement regardless of validity
    //  [4.1] each convex face is triangulated and the resulting triangles are stored with Vertex_2_handle
    //  [4.2] the correspondences between the triangles and the faces (centroids) are stored in a vector
    //  [4.3] non-convex faces are skipped
    //  [4.4] the centroid of each face is computed
    //  [4.5] determine which mesh triangle each centroid belongs to; if a centroid does not belong to any mesh triangle, then the entire face is invalid and should be skipped
    // [5] determine which power cell each centroid belongs to via a Kd-tree
    // [6] sweep through the triangles resulted from face triangulations
    //  [6.1] store vertex indices in an unordered_map
    //  [6.2] translate the Vertex_2_handle into vertex indices
    //  [6.3] store the mesh triangle membership and power cell membership indices using the correspondences between triangles and faces

    // [1] sweep through the mesh triangles and create Triangle_2 objects and Segment_2 objects
    // first construct Triangle_2 objects for each mesh triangle
    std::vector<Triangle_2> mesh_triangles;
    mesh_triangles.reserve(mesh_triangulation.size());

    // create an empty list and add all segments in the mesh to it
    std::forward_list<Segment_2> segments;

    for (tri_indices_t tri_indices : mesh_triangulation)
    {
        mesh_triangles.push_back(Triangle_2(mesh_vertices[std::get<0>(tri_indices)], mesh_vertices[std::get<1>(tri_indices)], \
            mesh_vertices[std::get<2>(tri_indices)]));
        
        // append the segments in the triangular mesh
        segments.push_front(Segment_2(mesh_vertices[std::get<0>(tri_indices)], mesh_vertices[std::get<1>(tri_indices)]));
        segments.push_front(Segment_2(mesh_vertices[std::get<1>(tri_indices)], mesh_vertices[std::get<2>(tri_indices)]));
        segments.push_front(Segment_2(mesh_vertices[std::get<2>(tri_indices)], mesh_vertices[std::get<0>(tri_indices)]));
    }

    // [2] compute the segments in the power diagram
    add_power_diagram_segments(circles, segments);

    // [3] compute the arrangement
    Arrangement arrangement;
    insert(arrangement, segments.begin(), segments.end());

    // [4] sweep through all faces in the arrangement regardless of validity
    // store the triangulation of the faces in a vector
    using tri_handle = std::tuple<Vertex_2_handle, Vertex_2_handle, Vertex_2_handle>;
    std::vector<tri_handle> face_triangulation;
    face_triangulation.reserve(arrangement.number_of_faces());

    // store the centroid of each face
    std::vector<Point_2> centroids;
    centroids.reserve(arrangement.number_of_faces());

    // store the correspondences between the triangles and the centroids
    std::vector<size_t> tri2centr_map;
    tri2centr_map.reserve(arrangement.number_of_faces());

    // store the correspondences between the centroids and the mesh triangles
    std::vector<size_t> centr2mesh_map;
    centr2mesh_map.reserve(arrangement.number_of_faces());
    
    // sweep through the faces
    size_t face_counter = 0;
    for (auto face_handle = arrangement.faces_begin(); face_handle != arrangement.faces_end(); ++face_handle)
    {
        // skip the unbounded face (there should be only 1)
        if (face_handle->is_unbounded())
            continue;

        // traverse around the outer CCB of the face
        auto outer_ccb_handle = face_handle->outer_ccb();
        auto first_vertex_handle = outer_ccb_handle->source();

        // for computing the centroid of the face (which is the arithmetic mean of the vertices)
        exact_field_t centroid_x = first_vertex_handle->point().x();
        exact_field_t centroid_y = first_vertex_handle->point().y();
        size_t face_tri_counter = 0;
        outer_ccb_handle = outer_ccb_handle->next();

        // for constructing a polygon corresponding to the face
        Polygon_2 face_polygon;
        face_polygon.push_back(first_vertex_handle->point());

        //  [4.1] each convex face is triangulated and the resulting triangles are stored with Vertex_2_handle
        while (true)
        {
            Vertex_2_handle edge_source_handle = outer_ccb_handle->source();
            Vertex_2_handle edge_target_handle = outer_ccb_handle->target();

            centroid_x += edge_source_handle->point().x();
            centroid_y += edge_source_handle->point().y();
            face_polygon.push_back(edge_source_handle->point());

            // assume first that the face is convex; if it is non-convex, it will be discarded later
            // since the face is a convex polygon, we adopt the fan triangulation
            if (edge_target_handle != first_vertex_handle)
            {
                // we add a triangle formed by the beginning vertex and the two ends of the edge
                face_triangulation.push_back(tri_handle{first_vertex_handle, edge_source_handle, edge_target_handle});

                //  [4.2] the correspondences between the triangles and the centroids are stored in vector
                // these might be removed later if the face is found out to be invalid
                tri2centr_map.push_back(face_counter);
                ++face_tri_counter;
            }
            else
            {
                // when the next vertex in the traversal equals the first one, the traversal of the face is complete
                break;
            }
            
            outer_ccb_handle = outer_ccb_handle->next();
        }

        //  [4.3] non-convex faces are skipped
        if (!face_polygon.is_convex())
        {
            // reset face_triangulation and tri2centr_map back to the states before encountering this non-convex face
            for (size_t reverse_id = 0; reverse_id < face_tri_counter; ++reverse_id)
            {
                face_triangulation.pop_back();
                tri2centr_map.pop_back();
            }
            
            continue;
        }

        //  [4.4] the centroid of each face is computed
        // the number of vertices in the face is equal to the number of triangles + 2
        Point_2 centroid = Point_2(centroid_x / exact_field_t(face_tri_counter + 2), centroid_y / exact_field_t(face_tri_counter + 2));

        //  [4.5] determine which mesh triangle each centroid belongs to; if a centroid does not belong to any mesh triangle, then the entire face is invalid and should be skipped
        bool is_face_valid = false;

        // sweep through all the mesh triangles
        for (size_t mesh_triangle_id = 0; mesh_triangle_id < mesh_triangles.size(); ++mesh_triangle_id)
        {
            if (mesh_triangles[mesh_triangle_id].has_on_bounded_side(centroid))
            {
                is_face_valid = true;
                centr2mesh_map[face_counter] = mesh_triangle_id;
                break;
            }
        }
        
        if (!is_face_valid)
        {
            // the centroid does not belong to any mesh triangle, and the entire face is invalid
            // reset face_triangulation and tri2centr_map back to the states before encountering this non-convex face
            for (size_t reverse_id = 0; reverse_id < face_tri_counter; ++reverse_id)
            {
                face_triangulation.pop_back();
                tri2centr_map.pop_back();
            }

            continue;
        }

        centroids.push_back(centroid);
        ++face_counter;
    }

    // [5] determine which power cell each centroid belongs to via a Kd-tree
    std::vector<size_t> centr2cell_map;
    compute_power_diagram_cell_indices(centroids, circles, centr2cell_map);

    
    // clear the output data structures before appending the outputs
    output_vertices.clear();
    output_triangulation.clear();
    power_cell_indices.clear();
    mesh_indices.clear();

    output_vertices.reserve(arrangement.number_of_vertices());
    output_triangulation.reserve(face_triangulation.size());
    power_cell_indices.reserve(face_triangulation.size());
    mesh_indices.reserve(face_triangulation.size());

    // [6] sweep through the triangles resulted from face triangulations

    // map for converting Vertex_2_handle to the corresponding index
    std::unordered_map<Vertex_2_handle, size_t> vert2id_map;

    size_t output_vertex_counter = 0;

    for (size_t face_triangle_id = 0; face_triangle_id < face_triangulation.size(); ++ face_triangle_id)
    {
        auto face_tri_handle = face_triangulation[face_triangle_id];

        //  [6.1] store vertex indices in an unordered_map
        if (!vert2id_map.contains(std::get<0>(face_tri_handle)))
        {
            vert2id_map.insert_or_assign(std::get<0>(face_tri_handle), output_vertex_counter);
            output_vertices.push_back(std::get<0>(face_tri_handle)->point());
            ++output_vertex_counter;
        }

        if (!vert2id_map.contains(std::get<1>(face_tri_handle)))
        {
            vert2id_map.insert_or_assign(std::get<1>(face_tri_handle), output_vertex_counter);
            output_vertices.push_back(std::get<1>(face_tri_handle)->point());
            ++output_vertex_counter;
        }

        if (!vert2id_map.contains(std::get<2>(face_tri_handle)))
        {
            vert2id_map.insert_or_assign(std::get<2>(face_tri_handle), output_vertex_counter);
            output_vertices.push_back(std::get<2>(face_tri_handle)->point());
            ++output_vertex_counter;
        }

        //  [6.2] translate the Vertex_2_handle into vertex indices
        output_triangulation.push_back(tri_indices_t(vert2id_map[std::get<0>(face_tri_handle)], \
            vert2id_map[std::get<1>(face_tri_handle)], vert2id_map[std::get<2>(face_tri_handle)]));

        //  [6.3] store the mesh triangle membership and power cell membership indices using the correspondences between triangles and faces
        power_cell_indices.push_back(centr2cell_map[tri2centr_map[face_triangle_id]]);
        mesh_indices.push_back(centr2mesh_map[tri2centr_map[face_triangle_id]]);
    }
}