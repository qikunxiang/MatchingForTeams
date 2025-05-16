#include <CGAL/convex_hull_2.h>
#include <CGAL/minkowski_sum_2.h>

#include <iostream>
#include <random>

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;

int main(int argc, char* argv[])
{
    Polygon_2 P, Q;

    std::vector<Point_2> P_points, Q_points;

    P_points.push_back(Point_2(0.0, 0.0));
    P_points.push_back(Point_2(1.0, 1.0));
    P_points.push_back(Point_2(1.0, 0.0));

    Q_points.push_back(Point_2(1.0, 0.0));
    Q_points.push_back(Point_2(0.0, -1.0));
    Q_points.push_back(Point_2(-1.0, 0.0));
    Q_points.push_back(Point_2(0.0, 1.0));

    CGAL::convex_hull_2(P_points.begin(), P_points.end(), std::back_inserter(P));
    CGAL::convex_hull_2(Q_points.begin(), Q_points.end(), std::back_inserter(Q));

    std::cout << "P = " << P << std::endl;
    std::cout << "Q = " << Q << std::endl;

    Polygon_with_holes_2 sum = CGAL::minkowski_sum_2(P, Q);

    std::cout << "sum = " << sum.outer_boundary() << std::endl;

    return 0;
}