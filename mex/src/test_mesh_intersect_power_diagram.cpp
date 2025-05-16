#include "power_diagram_intersection.hpp"
#include "mesh_intersect_power_diagram.hpp"

#include <iostream>
#include <random>

int main(int argc, char* argv[])
{
    // initialize the random number generator
    std::mt19937_64 rand_generator(1000);
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
    std::uniform_int_distribution<int> uniform_0_10_distribution(0, 10);

    // construct a triangular mesh
    std::vector<Point_2> mesh_vertices;
    mesh_vertices.reserve(9);

    mesh_vertices.push_back(Point_2(exact_field_t(0), exact_field_t(0)));
    mesh_vertices.push_back(Point_2(exact_field_t(5), exact_field_t(0)));
    mesh_vertices.push_back(Point_2(exact_field_t(10), exact_field_t(0)));
    mesh_vertices.push_back(Point_2(exact_field_t(0), exact_field_t(5)));
    mesh_vertices.push_back(Point_2(exact_field_t(5), exact_field_t(5)));
    mesh_vertices.push_back(Point_2(exact_field_t(10), exact_field_t(5)));
    mesh_vertices.push_back(Point_2(exact_field_t(0), exact_field_t(10)));
    mesh_vertices.push_back(Point_2(exact_field_t(5), exact_field_t(10)));
    mesh_vertices.push_back(Point_2(exact_field_t(10), exact_field_t(10)));

    std::vector<tri_indices_t> mesh_triangulation;
    mesh_triangulation.reserve(8);

    mesh_triangulation.push_back(tri_indices_t(0, 1, 4));
    mesh_triangulation.push_back(tri_indices_t(0, 3, 4));
    mesh_triangulation.push_back(tri_indices_t(1, 2, 5));
    mesh_triangulation.push_back(tri_indices_t(1, 4, 5));
    mesh_triangulation.push_back(tri_indices_t(3, 4, 7));
    mesh_triangulation.push_back(tri_indices_t(3, 6, 7));
    mesh_triangulation.push_back(tri_indices_t(4, 5, 8));
    mesh_triangulation.push_back(tri_indices_t(4, 7, 8));

    // output to a MATLAB .m script
    std::cout << "mesh_vertices = [ ..." << std::endl;

    for (Point_2 mesh_vertex : mesh_vertices)
    {
        std::cout << "    " << mesh_vertex.x().to_double() << ", " << mesh_vertex.y().to_double() << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl << std::endl;

    std::cout << "mesh_triangulation = [ ..." << std::endl;

    for (tri_indices_t mesh_tri : mesh_triangulation)
    {
        std::cout << "    " << std::get<0>(mesh_tri) + 1 << ", " << std::get<1>(mesh_tri) + 1 << ", " \
            << std::get<2>(mesh_tri) + 1 << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl << std::endl;

    // construct the collection of circle sets
    size_t circle_num = 200;

    typedef std::pair<Point_2, exact_field_t> circle_t;
    typedef std::vector<circle_t> circle_set_t;

    circle_set_t circles;

    for (size_t circle_id = 0; circle_id < circle_num; ++circle_id)
    {
        exact_field_t circle_radius = uniform_distribution(rand_generator);
        circles.push_back(circle_t(Point_2(exact_field_t(uniform_distribution(rand_generator) * 10.0), \
            exact_field_t(uniform_distribution(rand_generator) * 10.0)), \
            circle_radius * circle_radius));
    }

    // add the 4 artificial points to truncate the relevant rays to segments
    circles.push_back(circle_t(Point_2(exact_field_t(5), exact_field_t(-40)), exact_field_t(0)));
    circles.push_back(circle_t(Point_2(exact_field_t(5), exact_field_t(50)), exact_field_t(0)));
    circles.push_back(circle_t(Point_2(exact_field_t(-40), exact_field_t(5)), exact_field_t(0)));
    circles.push_back(circle_t(Point_2(exact_field_t(50), exact_field_t(5)), exact_field_t(0)));

    // output the MATLAB .m script
    std::cout << "circles = [ ..." << std::endl;

    for (auto circle : circles)
    {
        std::cout << "    " << circle.first.x().to_double() << ", " << circle.first.y().to_double() << ", " \
            << circle.second.to_double() << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl << std::endl;

    std::forward_list<Segment_2> segments;

    add_power_diagram_segments(circles, segments);

    std::cout << "segments = [ ..." << std::endl;

    for (Segment_2 segment : segments)
    {
        std::cout << "    " << segment.source().x().to_double() << ", " << segment.source().y().to_double() \
            << ", " << segment.target().x().to_double() << ", " << segment.target().y().to_double() << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl << std::endl;

    std::vector<Point_2> output_vertices;
    std::vector<tri_indices_t> output_triangulation;
    std::vector<size_t> power_cell_indices;
    std::vector<size_t> mesh_indices;

    compute_mesh_intersect_power_diagram(mesh_vertices, mesh_triangulation, circles, \
        output_vertices, output_triangulation, power_cell_indices, mesh_indices);
    
    // print the output to the same MATLAB .m file
    std::cout << "output_vertices = [ ..." << std::endl;

    for (Point_2 output_vertex : output_vertices)
    {
        std::cout << "    " << output_vertex.x().to_double() << ", " << output_vertex.y().to_double() << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl << std::endl;

    std::cout << "output_triangulation = [ ..." << std::endl;

    for (tri_indices_t tri : output_triangulation)
    {
        std::cout << "    " << std::get<0>(tri) + 1 << ", " << std::get<1>(tri) + 1 << ", " << std::get<2>(tri) + 1 << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl << std::endl;

    std::cout << "power_cell_indices = [ ..." << std::endl;

    for (size_t power_cell_index : power_cell_indices)
    {
        std::cout << "    " << power_cell_index + 1 << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl << std::endl;

    std::cout << "mesh_indices = [ ..." << std::endl;


    for (size_t mesh_index : mesh_indices)
    {
        std::cout << "    " << mesh_index + 1 << "; ..." << std::endl;
    }

    std::cout << "];" << std::endl << std::endl;

    return 0;
}