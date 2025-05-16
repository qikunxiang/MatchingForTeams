#ifndef __MESH_INTERSECT_POWER_DIAGRAM_HPP__
#define __MESH_INTERSECT_POWER_DIAGRAM_HPP__

#include "power_diagram_intersection.hpp"

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
void compute_mesh_intersect_power_diagram( \
    const std::vector<Point_2> & mesh_vertices, \
    const std::vector<tri_indices_t> & mesh_triangulation, \
    const std::vector<std::pair<Point_2, exact_field_t>> & circles, \
    std::vector<Point_2> & output_vertices, \
    std::vector<tri_indices_t> & output_triangulation, \
    std::vector<size_t> & power_cell_indices, \
    std::vector<size_t> & mesh_indices);

#endif // __MESH_INTERSECT_POWER_DIAGRAM_HPP__