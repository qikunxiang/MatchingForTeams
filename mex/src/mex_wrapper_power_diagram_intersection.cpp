#include "power_diagram_intersection.hpp"

#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

/// @brief For generating a MATLAB function: cell_indices = power_diagram_intersection(circles_cell, bbox_bottom_left, bbox_top_right)
/// This function takes multiple collections of circles as well as a bounding box and computes the intersections of the cells in the corresponding power diagrams formed by these collections of circles. It returns a list of indices where the intersections of the corresponding power cells have a non-empty intersection with the bounding box.
/// Reference: Altschuler, J.M. and Boix-Adserà, E., 2021. Wasserstein barycenters can be computed in polynomial time in fixed dimension. Journal of Machine Learning Research, 22(44), pp.1-19. See https://github.com/eboix/high_precision_barycenters/tree/master.
/// Inputs:
///     circles_cell: cell array with two columns where the first column contains two-column matrices representing the centers of the circles and 
///     the second column contains vectors representing the squared radii of the circles
///     bbox_bottom_left: two-element vector representing the coordinate of the bottom-left corner of the bounding box
///     bbox_top_right: two-element vector representing the coordinate of the top-right corner of the bounding box
/// Output:
///     cell_indices: a matrix where each row corresponds to a combination of cell indices which results in a non-empty intersection and each column
///     corresponds to a collection of circles
class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        // check the inputs
        if (outputs.size() > 1)
            throw matlab::data::InvalidNumberOfElementsProvidedException("Too many outputs");

        if (inputs.size() != 3)
            throw matlab::data::InvalidNumberOfElementsProvidedException("Inputs mis-specified: expecting 3 inputs");
        
        if (inputs[0].getType() != matlab::data::ArrayType::CELL)
            throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: circles_cell needs to be a cell array");

        const matlab::data::CellArray & circles_cell = inputs[0];
        
        size_t set_num = circles_cell.getDimensions()[0];

        if (circles_cell.getDimensions()[1] != 2)
            throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: circles_cell needs to have 2 columns");
        
        const matlab::data::TypedArray<double> & bbox_bottom_left = inputs[1];
        const matlab::data::TypedArray<double> & bbox_top_right = inputs[2];

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

        // this is an upper bound for the maximum distance within the bounding box; this is used to construct artificial points for truncating the rays in the power 
        // diagram to segments
        exact_field_t bbox_distance_bound = (bbox_x_max - bbox_x_min) + (bbox_y_max - bbox_y_min);

        // construct the collection of circle sets
        typedef std::pair<Point_2, exact_field_t> circle_t;
        typedef std::vector<circle_t> circle_set_t;
        std::vector<circle_set_t> circle_sets(set_num);
        std::vector<size_t> circle_num_list;
        circle_num_list.reserve(set_num);
        
        // access data in circles_cell
        for (size_t set_id = 0; set_id < set_num; ++set_id)
        {
            const matlab::data::TypedArray<double> & circle_centers = circles_cell[set_id][0];
            const matlab::data::TypedArray<double> & circle_squared_radii = circles_cell[set_id][1];

            size_t circle_num = circle_centers.getDimensions()[0];
            circle_num_list.push_back(circle_num);

            if (circle_centers.getDimensions()[1] != 2)
                throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: matrices in the first column of circles_cell need to have 2 columns");
            
            if (circle_squared_radii.getDimensions()[1] != 1)
                throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: matrices in the second column of circles_cell need to have 1 column");
            
            if (circle_squared_radii.getDimensions()[0] != circle_num)
                throw matlab::data::InvalidArrayTypeException("Inputs mis-specified: matrices in the second column of circles_cell need to have the same \
                    number of rows as the matrices in the first column of circles_cell");
            
            // add the circles to the list
            for (size_t circle_id = 0; circle_id < circle_num; ++circle_id)
            {
                circle_sets[set_id].push_back(circle_t(Point_2(exact_field_t((double) (circle_centers[circle_id][0])), \
                    exact_field_t((double) (circle_centers[circle_id][1]))), exact_field_t((double) (circle_squared_radii[circle_id]))));
            }

            // add 4 artificial points to the power diagram in order to truncate the relevant rays in the power diagram into segments; as a result
            // rays in the power diagram can now be ignored
            circle_sets[set_id].push_back(circle_t(Point_2(bbox_x_mid, bbox_y_max + bbox_distance_bound * exact_field_t(2)), exact_field_t(0)));
            circle_sets[set_id].push_back(circle_t(Point_2(bbox_x_mid, bbox_y_min - bbox_distance_bound * exact_field_t(2)), exact_field_t(0)));
            circle_sets[set_id].push_back(circle_t(Point_2(bbox_x_max + bbox_distance_bound * exact_field_t(2), bbox_y_mid), exact_field_t(0)));
            circle_sets[set_id].push_back(circle_t(Point_2(bbox_x_min - bbox_distance_bound * exact_field_t(2), bbox_y_mid), exact_field_t(0)));
        }

        // create an empty list and add all segments in the power diagrams to it
        std::forward_list<Segment_2> segments;

        for (const circle_set_t & circle_set : circle_sets)
        {
            add_power_diagram_segments(circle_set, segments);
        }

        // compute a list of centroid of the non-empty power cells
        std::vector<Point_2> centroids;
        compute_power_diagram_intersection(segments, centroids);

        // store all the indices first
        std::vector<std::vector<size_t>> centroid_cell_indices(set_num);

        // find the indices of the cells that each centroid belongs to while eliminating irrelevant centroids;
        for (size_t set_id = 0; set_id < set_num; ++set_id)
        {
            // remove the artificial points for each set of circles
            circle_sets[set_id].pop_back();
            circle_sets[set_id].pop_back();
            circle_sets[set_id].pop_back();
            circle_sets[set_id].pop_back();

            compute_power_diagram_cell_indices(centroids, circle_sets[set_id], centroid_cell_indices[set_id]);
        }

        // construct output
        matlab::data::ArrayFactory factory;

        auto & cell_indices = outputs[0];
        cell_indices = factory.createArray<double>({centroids.size(), set_num});
        
        for (size_t centroid_id = 0; centroid_id < centroids.size(); ++centroid_id)
        {
            // assign the output matrix row by row
            for (size_t set_id = 0; set_id < set_num; ++set_id)
            {
                cell_indices[centroid_id][set_id] = (double) (centroid_cell_indices[set_id][centroid_id] + 1);
            }
        }
    }
};