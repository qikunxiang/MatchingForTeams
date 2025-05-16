#include <CGAL/Gmpq.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/minkowski_sum_2.h>

#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

using FT = CGAL::Gmpq;
using Kernel = CGAL::Simple_cartesian<FT>;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

/// @brief For generating a MATLAB function: output_vertices = minkowski_sum_2d(input_vertices_cell).
/// The function computes the Minkowski sum of several two-dimensional polygons.
/// Inputs:
///     input_vertices_cell: cell array containing two-column matrices each representing the vertices of the input polygons
/// Outputs:
///     output_vertices: two-column matrix containing the vertices in the output polygon
class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        // check the inputs
        if (outputs.size() > 1)
            throw matlab::data::InvalidNumberOfElementsProvidedException("Too many outputs");

        if (inputs.size() != 1)
            throw matlab::data::InvalidNumberOfElementsProvidedException("Inputs mis-specified: expecting 1 input");
        
        const matlab::data::CellArray & input_vertices_cell = inputs[0];
        
        size_t polygon_num = input_vertices_cell.getDimensions()[0];

        Polygon_2 output_polygon;

        // sweep through the input polygons
        for (size_t p_id = 0; p_id < polygon_num; ++p_id)
        {
            const matlab::data::TypedArray<double> & matlab_input_vertices = input_vertices_cell[p_id][0];
            std::vector<Point_2> input_vertices;

            size_t vertex_num = matlab_input_vertices.getDimensions()[0];

            // build a Polygon_2 object representing the current input polygon
            for (size_t v_id = 0; v_id < vertex_num; ++v_id)
            {
                input_vertices.push_back(Point_2(FT((double) matlab_input_vertices[v_id][0]), \
                    FT((double) matlab_input_vertices[v_id][1])));
            }

            if (0 == p_id)
                CGAL::convex_hull_2(input_vertices.begin(), input_vertices.end(), std::back_insert_iterator(output_polygon));
            else
            {
                Polygon_2 input_polygon;
                CGAL::convex_hull_2(input_vertices.begin(), input_vertices.end(), std::back_insert_iterator(input_polygon));
                
                auto sum = CGAL::minkowski_sum_2(output_polygon, input_polygon);
                output_polygon = sum.outer_boundary();
            }
        }

        // output to MATLAB
        matlab::data::ArrayFactory factory;

        auto & output_vertices = outputs[0];
        output_vertices = factory.createArray<double>({output_polygon.size(), 2});
        
        for (size_t v_id = 0; v_id < output_polygon.size(); ++v_id)
        {
            output_vertices[v_id][0] = output_polygon[v_id].x().to_double();
            output_vertices[v_id][1] = output_polygon[v_id].y().to_double();
        }
    }
};