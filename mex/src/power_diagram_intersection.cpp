#include "power_diagram_intersection.hpp"

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

/// @brief Compute the power diagram of a collection of circles with given centers and radii. 
/// Segments in the resulting power diagram are added to a given list and the rays are ignored.
/// Reference: Aurenhammer, F., 1987. Power diagrams: properties, algorithms and applications. SIAM Journal on Computing, 16(1), pp.78-96.
/// @param circles [input] std::vector containing objects of std::pair which contains the centers and squared radii of circles
/// @param segments [output] std::forward_list for appending the segments at the back
void add_power_diagram_segments(const std::vector<std::pair<Point_2, exact_field_t>> & circles, \
    std::forward_list<Segment_2> & segments)
{
    size_t circle_num = circles.size();

    // create the container of the lifted points; these points in 3D are dual to the hyperplanes corresponding to the circles
    std::vector<Point_3> lifted_points;
    lifted_points.reserve(circle_num);

    for (auto circle = circles.begin(); circle != circles.end(); ++circle)
    {
        Point_2 center = circle->first;
        exact_field_t squared_radius = circle->second;

        // the radius need to be non-negative
        CGAL_precondition(squared_radius >= 0);

        // suppose that the circle s has center (s1, s2) and radius r, then the circle corresponds to a hyperplane:
        // h = \Pi(s): x3 = 2 * s1 * x1 + 2 * s2 * x2 - s1^2 - s2^2 + r^2 and \Delta(h) = (s1, s2, s1^2 + s2^2 - r^2);
        // see Section 4.1 of Aurenhammer 1987
        exact_field_t lifted_x = center.x();
        exact_field_t lifted_y = center.y();
        exact_field_t lifted_z = lifted_x * lifted_x + lifted_y * lifted_y - squared_radius;

        // append the lifted point
        lifted_points.push_back(Point_3(lifted_x, lifted_y, lifted_z));
    }

    // create the 3D polyhedron to represent the convex hull of the lifted points
    Polyhedron_3 convhull;

    // compute the convex hull
    CGAL::convex_hull_3(lifted_points.begin(), lifted_points.end(), convhull);


    // create the map for storing the correspondence between facets of the convex hull and points in the power diagram
    std::unordered_map<Facet_3_handle, Point_2> facet_map;

    // traverse the facets of the convex hull
    for (Facet_3_handle facet_handle = convhull.facets_begin(); facet_handle != convhull.facets_end(); ++facet_handle)
    {
        const Facet_3 & facet = *facet_handle;

        // the facet should be a triangle
        CGAL_precondition(facet.is_triangle());

        auto halfedge_handle = facet.halfedge();
        
        // construct a triangle with the three vertices of the facet in order to compute the equation
        Triangle_3 facet_triangle = Triangle_3(halfedge_handle->vertex()->point(), \
            halfedge_handle->prev()->vertex()->point(), \
            halfedge_handle->prev()->prev()->vertex()->point());

        Plane_3 facet_plane = facet_triangle.supporting_plane();

        // skip all facets that are perpendicular to the xy-plane or facing downwards
        if (facet_plane.c() <= 0)
            continue;

        // insert the point in the power diagram that corresponds to the facet
        // suppose that the facet is contained in the hyperplane h: x3 = -a / c * x1 - b / c * x2 - d / c, then 
        // \Deta(h) = (-a / c / 2, -b / c / 2, d / c), and we will discard the third coordinate
        exact_field_t power_diagram_point_x = -facet_plane.a() / facet_plane.c() / 2;
        exact_field_t power_diagram_point_y = -facet_plane.b() / facet_plane.c() / 2;
        facet_map.insert_or_assign(facet_handle, Point_2(power_diagram_point_x, power_diagram_point_y));
        
        // traverse the neighboring facets and store the segments, in exist
        for (size_t halfedge_index = 0; halfedge_index < 3; ++halfedge_index)
        {
            // retrieve the neighboring facet through the opposite halfedge
            Facet_3_handle neighbor_facet_handle = halfedge_handle->opposite()->facet();

            try
            {
                const Point_2 & neighbor_point = facet_map.at(neighbor_facet_handle);
                Segment_2 new_segment = Segment_2(Point_2(power_diagram_point_x, power_diagram_point_y), neighbor_point);

                // add a segment that connects the power diagram point corresponding to the current facet with 
                // the power diagram point corresponding to the neighboring facet
                // only add the segment if it is non-degenerate
                if (!new_segment.is_degenerate())
                    segments.push_front(std::move(new_segment));
            }
            catch (const std::out_of_range e)
            {
                // the neighbor facet has not been traversed yet, or is not facing upwards
            }

            // move to the next halfedge
            halfedge_handle = halfedge_handle->next();
        }
    }
}


/// @brief Compute the polyhedral cell complex formed by the intersection of power diagrams. Rather than returning the resulting polyhedral 
/// cell complex, the centroid of each non-empty cell is returned in a std::vector. 
/// @param segments [input] a list containing all segments in the power diagrams; rays are ignored and thus artificial extra points need to be added to the construction of power diagrams to truncate the rays to segments
/// @param centroids [output] a std::vector containing the centroid points of the resulting polyhedral cell complex
void compute_power_diagram_intersection(const std::forward_list<Segment_2> & segments, \
    std::vector<Point_2> & centroids)
{
    // build the arrangement from the segments
    Arrangement arrangement;
    insert(arrangement, segments.begin(), segments.end());

    centroids.clear();
    centroids.reserve(arrangement.number_of_faces());

    for (auto face_handle = arrangement.faces_begin(); face_handle != arrangement.faces_end(); ++face_handle)
    {
        // skip the unbounded face (there should be only 1)
        if (face_handle->is_unbounded())
            continue;

        // start computing the centroid of the face (which is simply the arithmetic mean of the vertices)
        exact_field_t centroid_x = 0;
        exact_field_t centroid_y = 0;
        int face_vertex_counter = 0;

        // traverse around the outer CCB of the face
        auto outer_ccb_handle = face_handle->outer_ccb();
        auto first_vertex_handle = outer_ccb_handle->source();

        while (true)
        {
            centroid_x += outer_ccb_handle->source()->point().x();
            centroid_y += outer_ccb_handle->source()->point().y();
            ++face_vertex_counter;

            // when the next vertex in the traversal equals the first one, the traversal is complete
            if (outer_ccb_handle->target() == first_vertex_handle)
                break;
            
            outer_ccb_handle = outer_ccb_handle->next();
        }

        centroids.push_back(Point_2(centroid_x / face_vertex_counter, centroid_y / face_vertex_counter));
    }
}

/// @brief Given a collection of 2D points, query which cells in the power diagram they belong to.
/// The use of K-d tree (K-dimensional tree) was inspired by Altschuler, J.M. and Boix-Adserà, E., 2021. Wasserstein barycenters can be computed in polynomial time in fixed dimension. Journal of Machine Learning Research, 22(44), pp.1-19. See https://github.com/eboix/high_precision_barycenters/tree/master.
/// @param query_points [input] std::vector containing the points to be queried
/// @param circles [input] std::vector containing objects of std::pair which contains the centers and squared radii of circles
/// @param power_cell_indices [output] std::vector containing the indices of the power diagram cells that the query points belong to
void compute_power_diagram_cell_indices(const std::vector<Point_2> & query_points, \
    const std::vector<std::pair<Point_2, exact_field_t>> & circles, \
    std::vector<size_t> & power_cell_indices)
{
    using Kd_tree_K = CGAL::Simple_cartesian<double>;
    using Kd_tree_Point_3 = Kd_tree_K::Point_3;
    using TreeTraits_base = CGAL::Search_traits_3<Kd_tree_K>;
    using TreeTraits = CGAL::Search_traits_adapter<std::tuple<Kd_tree_Point_3, size_t>, 
        CGAL::Nth_of_tuple_property_map<0, std::tuple<Kd_tree_Point_3, size_t>>, TreeTraits_base>;
    using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<TreeTraits>;
    using Tree = Neighbor_search::Tree;

    // first compute the maximum radius of the circles
    double circle_squared_radius_max = -1;

    for (std::pair<Point_2, exact_field_t> circle : circles)
    {
        double current_circle_squared_radius = circle.second.to_double();
        if (current_circle_squared_radius > circle_squared_radius_max)
        {
            circle_squared_radius_max = current_circle_squared_radius;
        }
    }

    // construct Kd-tree
    std::vector<std::tuple<Kd_tree_Point_3, size_t>> tree_points_with_indices;
    tree_points_with_indices.reserve(circles.size());

    for (size_t circle_id = 0; circle_id < circles.size(); ++circle_id)
    {
        auto circle = circles[circle_id];

        // the x- and y-coordinates of the point in the K-d tree are identical to the center of the circle;
        // the z-coordinate of the point in the K-d tree is equal to sqrt(r_max^2 - r^2) where r is the radius of the circle;
        // subsequently, the squared Euclidean distance between (x, y, 0) and (sx, sy, sqrt(radius_max^2 - radius^2)) equals (x - sx)^2 + (y - sy)^2 + r_max^2 - r^2; since r_max is constant, this is equivalent to using the squared Euclidean distance (x - sx)^2 + (y - sy)^2 - r^2 for determining the nearest neighbor, which corresponds to finding the power diagram cell that a point belongs to
        double circle_center_x = circle.first.x().to_double();
        double circle_center_y = circle.first.y().to_double();
        double circle_squared_radius = circle.second.to_double();
        tree_points_with_indices.push_back(std::tuple<Kd_tree_Point_3, size_t>(Kd_tree_Point_3(circle_center_x, circle_center_y, \
            std::sqrt(circle_squared_radius_max - circle_squared_radius)), circle_id));
    }

    Tree tree(tree_points_with_indices.begin(), tree_points_with_indices.end());

    size_t point_num = query_points.size();

    power_cell_indices.clear();
    power_cell_indices.reserve(point_num);

    for (Point_2 point : query_points)
    {
        Kd_tree_Point_3 query(point.x().to_double(), point.y().to_double(), 0.0);
        Neighbor_search search(tree, query, 1);

        power_cell_indices.push_back(std::get<1>(search.begin()->first));
    }
}