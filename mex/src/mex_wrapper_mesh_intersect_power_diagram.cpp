#include "power_diagram_intersection.hpp"
#include "mesh_intersect_power_diagram.hpp"

#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

/// @brief For generating a MATLAB function: [output_vertices, output_triangulation, power_cell_indices, mesh_indices] = mesh_intersect_power_diagram(mesh_vertices, mesh_triangulation, circle_centers, circle_radii, bbox_bottom_left, bbox_top_right)
/// This function takes a triangular mesh, a collection of circles, and a bounding box and computes the intersections of the triangles in the mesh with 
/// the cells in the power diagram. It returns triangles formed via triangulating the resulting intersections as well as which mesh triangle and which power
/// cell each of the output triangles belongs to. 
/// Reference: Altschuler, J.M. and Boix-Adserà, E., 2021. Wasserstein barycenters can be computed in polynomial time in fixed dimension. Journal of Machine Learning Research, 22(44), pp.1-19. See https://github.com/eboix/high_precision_barycenters/tree/master.
/// Inputs:
///     mesh_vertices: two-column matrix containing the vertices of the triangular mesh
///     mesh_triangulation: three-column matrix representing the mesh triangulation
///     circle_centers: two-column matrix containing the centers of the circles
///     circle_squared_radii: vector containing the squared radii of the circles
///     bbox_bottom_left: two-element vector representing the coordinate of the bottom-left corner of the bounding box; the bounding box must contain the entire triangular mesh
///     bbox_top_right: two-element vector representing the coordinate of the top-right corner of the bounding box; the bounding box must contain the entire triangular mesh
/// Outputs:
///     output_vertices: two-column matrix containing the vertices in the output triangular mesh
///     output_triangulation: three-column matrix representing the output triangulation
///     power_cell_indices: vector containing the power cell indices of the output triangles
///     mesh_indices: vector containing the mesh triangle indices of the output triangles
class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        // check the inputs
        if (outputs.size() > 4)
            throw matlab::data::InvalidNumberOfElementsProvidedException("Too many outputs");

        if (inputs.size() != 6)
            throw matlab::data::InvalidNumberOfElementsProvidedException("Inputs mis-specified: expecting 6 inputs");
        
        const matlab::data::TypedArray<double> & matlab_mesh_vertices = inputs[0];
        const matlab::data::TypedArray<double> & matlab_mesh_triangulation = inputs[1];
        const matlab::data::TypedArray<double> & circle_centers = inputs[2];
        const matlab::data::TypedArray<double> & circle_squared_radii = inputs[3];
        const matlab::data::TypedArray<double> & bbox_bottom_left = inputs[4];
        const matlab::data::TypedArray<double> & bbox_top_right = inputs[5];

        size_t mesh_vertex_num = matlab_mesh_vertices.getDimensions()[0];
        size_t mesh_triangle_num = matlab_mesh_triangulation.getDimensions()[0];
        size_t circle_num = circle_centers.getDimensions()[0];

        if (matlab_mesh_vertices.getDimensions()[1] != 2)
            throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: mesh_vertices needs to have 2 columns");
        
        if (matlab_mesh_triangulation.getDimensions()[1] != 3)
            throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: mesh_triangulation needs to have 3 columns");
        
        if (circle_centers.getDimensions()[1] != 2)
            throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: circle_centers needs to have 2 columns");

        if (circle_squared_radii.getNumberOfElements() != circle_num)
            throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: circle_squared_radii needs to be a vector with the \
                same number of rows as circle_centers");

        if (bbox_bottom_left.getNumberOfElements() != 2)
            throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: bbox_bottom_left needs to have 2 elements");

        if (bbox_top_right.getNumberOfElements() != 2)
            throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: bbox_top_right needs to have 2 elements");

        exact_field_t bbox_x_min = (double) bbox_bottom_left[0];
        exact_field_t bbox_x_max = (double) bbox_top_right[0];
        exact_field_t bbox_y_min = (double) bbox_bottom_left[1];
        exact_field_t bbox_y_max = (double) bbox_top_right[1];
        exact_field_t bbox_x_mid = (bbox_x_min + bbox_x_max) / exact_field_t(2);
        exact_field_t bbox_y_mid = (bbox_y_min + bbox_y_max) / exact_field_t(2);

        // this is an upper bound for the maximum distance within the bounding box; this is used to construct artificial points for truncating the rays in the power diagram to segments
        exact_field_t bbox_distance_bound = (bbox_x_max - bbox_x_min) + (bbox_y_max - bbox_y_min);

        // construct the data structures of the triangular mesh
        std::vector<Point_2> mesh_vertices;
        mesh_vertices.reserve(mesh_vertex_num);

        for (size_t vertex_id = 0; vertex_id < mesh_vertex_num; ++vertex_id)
        {
            mesh_vertices.push_back(Point_2(exact_field_t((double) (matlab_mesh_vertices[vertex_id][0])), \
                exact_field_t((double) (matlab_mesh_vertices[vertex_id][1]))));
        }

        std::vector<tri_indices_t> mesh_triangulation;
        mesh_triangulation.reserve(mesh_triangle_num);

        for (size_t triangle_id = 0; triangle_id < mesh_triangle_num; ++triangle_id)
        {
            mesh_triangulation.push_back(tri_indices_t((size_t) (matlab_mesh_triangulation[triangle_id][0]) - 1, \
                (size_t) (matlab_mesh_triangulation[triangle_id][1]) - 1, (size_t) (matlab_mesh_triangulation[triangle_id][2]) - 1));
        }

        // construct the collection of circles
        typedef std::pair<Point_2, exact_field_t> circle_t;
        typedef std::vector<circle_t> circle_set_t;
        circle_set_t circles;
        circles.reserve(circle_num);
        
        // add the circles to the list
        for (size_t circle_id = 0; circle_id < circle_num; ++circle_id)
        {
            circles.push_back(circle_t(Point_2(exact_field_t((double) (circle_centers[circle_id][0])), \
                exact_field_t((double) (circle_centers[circle_id][1]))), exact_field_t((double) (circle_squared_radii[circle_id]))));
        }

        // add 4 artificial points to the power diagram in order to truncate the relevant rays in the power diagram into segments; as a result rays in the power diagram can now be ignored
        circles.push_back(circle_t(Point_2(bbox_x_mid, bbox_y_max + bbox_distance_bound * exact_field_t(2)), exact_field_t(0)));
        circles.push_back(circle_t(Point_2(bbox_x_mid, bbox_y_min - bbox_distance_bound * exact_field_t(2)), exact_field_t(0)));
        circles.push_back(circle_t(Point_2(bbox_x_max + bbox_distance_bound * exact_field_t(2), bbox_y_mid), exact_field_t(0)));
        circles.push_back(circle_t(Point_2(bbox_x_min - bbox_distance_bound * exact_field_t(2), bbox_y_mid), exact_field_t(0)));

        // initialize the output data structures
        std::vector<Point_2> output_vertices;
        std::vector<tri_indices_t> output_triangulation;
        std::vector<size_t> power_cell_indices;
        std::vector<size_t> mesh_indices;

        // compute the output
        // note that since the bounding box is required to contain the entire triangular mesh and cells outside the triangular mesh are ignored, there will not be output triangles that belong to the 4 artificial power cells
        compute_mesh_intersect_power_diagram(mesh_vertices, mesh_triangulation, circles, output_vertices, output_triangulation, power_cell_indices, mesh_indices);

        // construct output
        matlab::data::ArrayFactory factory;

        auto & matlab_output_vertices = outputs[0];
        auto & matlab_output_triangulation = outputs[1];
        auto & matlab_power_cell_indices = outputs[2];
        auto & matlab_mesh_indices = outputs[3];

        matlab_output_vertices = factory.createArray<double>({output_vertices.size(), 2});
        
        for (size_t vertex_id = 0; vertex_id < output_vertices.size(); ++vertex_id)
        {
            matlab_output_vertices[vertex_id][0] = output_vertices[vertex_id].x().to_double();
            matlab_output_vertices[vertex_id][1] = output_vertices[vertex_id].y().to_double();
        }

        matlab_output_triangulation = factory.createArray<double>({output_triangulation.size(), 3});
        matlab_power_cell_indices = factory.createArray<double>({output_triangulation.size(), 1});
        matlab_mesh_indices = factory.createArray<double>({output_triangulation.size(), 1});

        for (size_t triangle_id = 0; triangle_id < output_triangulation.size(); ++triangle_id)
        {
            matlab_output_triangulation[triangle_id][0] = (double) (std::get<0>(output_triangulation[triangle_id]) + 1);
            matlab_output_triangulation[triangle_id][1] = (double) (std::get<1>(output_triangulation[triangle_id]) + 1);
            matlab_output_triangulation[triangle_id][2] = (double) (std::get<2>(output_triangulation[triangle_id]) + 1);

            matlab_power_cell_indices[triangle_id][0] = (double) (power_cell_indices[triangle_id] + 1);
            matlab_mesh_indices[triangle_id][0] = (double) (mesh_indices[triangle_id] + 1);
        }
    }
};