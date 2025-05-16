#include "power_diagram_intersection.hpp"

#include <iostream>
#include <random>
#include <array>

int main(int argc, char* argv[])
{
    // number of sets of circles
    const size_t set_num = 3;

    // minimum and maximum number of circles per set
    const size_t circle_num_per_set_min = 10;
    const size_t circle_num_per_set_max = 20;

    // initialize the random number generator
    std::mt19937_64 rand_generator(1000);
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
    std::uniform_int_distribution<int> uniform_0_10_distribution(0, 10);
    std::uniform_int_distribution<int> circle_num_distribution(circle_num_per_set_min, circle_num_per_set_max);

    // construct the collection of circle sets
    typedef std::pair<Point_2, exact_field_t> circle_t;
    typedef std::vector<circle_t> circle_set_t;
    std::array<circle_set_t, set_num> circle_sets;
    std::vector<size_t> circle_num_list;

    // output a MATLAB .m script
    std::cout << "circles_cell = cell(" << set_num << ", 1);" << std::endl;

    for (size_t set_id = 0; set_id < set_num; ++set_id)
    {
        size_t current_circle_num = circle_num_distribution(rand_generator);
        circle_num_list.push_back(current_circle_num);
        circle_sets[set_id].reserve(current_circle_num);

        for (size_t circle_id = 0; circle_id < current_circle_num; ++circle_id)
        {
            exact_field_t circle_radius = uniform_distribution(rand_generator);
            circle_sets[set_id].push_back(circle_t(Point_2(exact_field_t(uniform_0_10_distribution(rand_generator)), \
                exact_field_t(uniform_0_10_distribution(rand_generator))), \
                circle_radius * circle_radius));
        }

        // add the 4 artificial points to truncate the relevant rays to segments
        circle_sets[set_id].push_back(circle_t(Point_2(exact_field_t(5), exact_field_t(-40)), exact_field_t(0)));
        circle_sets[set_id].push_back(circle_t(Point_2(exact_field_t(5), exact_field_t(50)), exact_field_t(0)));
        circle_sets[set_id].push_back(circle_t(Point_2(exact_field_t(-40), exact_field_t(5)), exact_field_t(0)));
        circle_sets[set_id].push_back(circle_t(Point_2(exact_field_t(50), exact_field_t(5)), exact_field_t(0)));

        std::cout << "circles_cell{" << (set_id + 1) << "} = [ ..." << std::endl;

        for (auto circle : circle_sets[set_id])
        {
            std::cout << "    " << circle.first.x().to_double() << ", " << circle.first.y().to_double() << ", " \
                << circle.second.to_double() << "; ..." << std::endl;
        }

        std::cout << "];" << std::endl;
    }

    std::cout << std::endl;

    // create an empty list
    std::forward_list<Segment_2> segments;
    auto segment_end_handle = segments.end();

    std::cout << "segments_cell = cell(" << set_num << ", 1);" << std::endl;
    size_t circle_counter = 0;

    for (circle_set_t circle_set : circle_sets)
    {
        add_power_diagram_segments(circle_set, segments);

        std::cout << "segments_cell{" << ++circle_counter << "} = [ ..." << std::endl;

        for (auto segment_handle = segments.begin(); segment_handle != segment_end_handle; ++segment_handle)
        {
            std::cout << "   " << segment_handle->source().x().to_double() << ", " << segment_handle->source().y().to_double() \
                << ", " << segment_handle->target().x().to_double() << ", " << segment_handle->target().y().to_double() << "; ..." << std::endl;
        }

        std::cout << "];" << std::endl;

        segment_end_handle = segments.begin();
    }

    std::cout << std::endl;

    std::vector<Point_2> centroids;
    compute_power_diagram_intersection(segments, centroids);

    std::cout << "centroids = [ ..." << std::endl;

    for (auto centroid : centroids)
    {
        std::cout << "    " << centroid.x().to_double() << ", " << centroid.y().to_double() << "; ..." << std::endl;
    }
    
    std::cout << "];" << std::endl;

    // store all the indices first
    std::vector<std::vector<size_t>> centroid_cell_indices(set_num);

    // find the indices of the cells that each centroid belongs to while eliminating irrelevant centroids;
    for (size_t set_id = 0; set_id < set_num; ++set_id)
    {
        circle_sets[set_id].pop_back();
        circle_sets[set_id].pop_back();
        circle_sets[set_id].pop_back();
        circle_sets[set_id].pop_back();

        compute_power_diagram_cell_indices(centroids, circle_sets[set_id], centroid_cell_indices[set_id]);
    }

    std::cout << "cell_indices = [ ..." << std::endl;

    for (size_t centroid_id = 0; centroid_id < centroids.size(); ++centroid_id)
    {
        std::cout << "    ";
        for (size_t set_id = 0; set_id < set_num; ++set_id)
        {
            std::cout << centroid_cell_indices[set_id][centroid_id] + 1;

            if (set_id != set_num - 1)
                std::cout << ", ";
        }

        std::cout << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl;

    return 0;
}