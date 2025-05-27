classdef ProbMeas2DCPWADens < ProbMeas2DConvexPolytopeWithW2OT ...
        & HasTractableQuadraticIntegrals
    % Class for probability measures with continuous piece-wise affine (CPWA) density function. The support of the probability measure 
    % is formed by the union of triangles. The triangles need to satisfy the following condition: if two triangles A and B are not 
    % disjoint, then the intersection of A and B needs to be a face of both A and B, i.e., the intersection of A and B is either an 
    % edge of A or a vertex of A, and is also either an edge of B or a vertex of B. 

    methods(Access = public)
        function obj = ProbMeas2DCPWADens( ...
                vertices, ...
                triangles, ...
                dens_vertices)
            % Constructor function
            % Inputs:
            %   vertices: two-column matrix containing the vertices of all triangles
            %   triangles: three-column matrix indicating the triangulation dens_triangles: vector containing the density at each
            %   vertex (thus specifying the entire density function via interpolation)

            % call the superclass constructor to initialize the support
            obj@ProbMeas2DConvexPolytopeWithW2OT(vertices);

            obj.Dens = struct;

            if any(dens_vertices < 0)
                error('the density is not non-negative');
            end

            if any(all(dens_vertices(triangles') <= 1e-10, 1))
                error('there are triangles with 0 probability');
            end

            % check if there are duplicate vertices
            assert(size(unique(round(vertices, 6), 'rows'), 1) == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not appear in any triangle)
            assert(length(unique(triangles(:))) == size(vertices, 1), 'there are redundant vertices');

            % check that the triangles are non-empty (have non-zero area) and satisfy the condition
            [empty, is_simpcover] = check_triangulation_simplicial_cover(vertices, triangles);

            assert(all(~empty), 'there are empty triangles');
            assert(is_simpcover, 'the condition of simplicial cover is violated');

            obj.Supp.TriangleVertices = vertices;
            obj.Supp.Triangles = triangles;

            tri_num = size(triangles, 1);

            % evaluate the integral of the unnormalized density within each triangular region
            tri_integral = zeros(tri_num, 1);

            % evaluate the integrals of x1, x2, x1^2, x2^2, and x1*x2 with respect to the measure within each triangular region
            tri_lin_integral_x1 = zeros(tri_num, 1);
            tri_lin_integral_x2 = zeros(tri_num, 1);
            tri_quad_integral_x1sq = zeros(tri_num, 1);
            tri_quad_integral_x2sq = zeros(tri_num, 1);
            tri_quad_integral_x1x2 = zeros(tri_num, 1);

            for tri_id = 1:tri_num
                tri_verts = vertices(triangles(tri_id, :), :);
                dens_triangle = dens_vertices(triangles(tri_id, :));

                int_vals = ProbMeas2DAffDens.affineIntegralTriangle( ...
                        repmat(dens_triangle', 6, 1), ...
                        repmat(tri_verts(:, 1)', 6, 1), ...
                        repmat(tri_verts(:, 2)', 6, 1), ...
                        eye(6), ...
                        false(6, 1));

                tri_integral(tri_id) = int_vals(6);
                tri_lin_integral_x1(tri_id) = int_vals(4);
                tri_lin_integral_x2(tri_id) = int_vals(5);
                tri_quad_integral_x1sq(tri_id) = int_vals(1);
                tri_quad_integral_x2sq(tri_id) = int_vals(2);
                tri_quad_integral_x1x2(tri_id) = int_vals(3);
            end

            % calculate the normalizing constant and normalize the density
            norm_const = sum(tri_integral);
            dens_vertices = dens_vertices / norm_const;
            obj.Dens = struct('VertexDensities', dens_vertices, 'NormConst', norm_const);

            % calculate the expected value of x1, x2, x1^2, x2^2, x1*x2
            obj.FirstMomentVec = zeros(2, 1);
            obj.FirstMomentVec(1) = sum(tri_lin_integral_x1) / norm_const;
            obj.FirstMomentVec(2) = sum(tri_lin_integral_x2) / norm_const;
            obj.SecondMomentMat = zeros(2, 2);
            obj.SecondMomentMat(1, 1) = sum(tri_quad_integral_x1sq) / norm_const;
            obj.SecondMomentMat(2, 2) = sum(tri_quad_integral_x2sq) / norm_const;
            obj.SecondMomentMat(1, 2) = sum(tri_quad_integral_x1x2) / norm_const;
            obj.SecondMomentMat(2, 1) = obj.SecondMomentMat(1, 2);

            obj.Dens.InterpMat = zeros(tri_num, 3);
            obj.Dens.EdgeIneq = struct;
            
            ineq_w_cell = cell(tri_num, 1);
            ineq_b_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_verts = triangles(tri_id, :)';

                tri_vert_mat = vertices(tri_verts, :);

                % calculate the affine interpolation function for each triangle (ignoring the fact that interpolation only makes sense 
                % within the triangle)
                obj.Dens.InterpMat(tri_id, :) = [tri_vert_mat, ones(3, 1)] \ dens_vertices(tri_verts);

                % sort the three vertices counterclockwise
                ccorder = convhull(tri_vert_mat(:, 1), tri_vert_mat(:, 2), 'Simplify', true);

                % there are now 4 rows, the first and the last are identical
                tri_vert_cc = tri_vert_mat(ccorder, :);

                tri_vert1 = tri_vert_cc(1:3, :);
                tri_vert2 = tri_vert_cc(2:4, :);

                ineq_w_cell{tri_id} = [tri_vert2(:, 2) - tri_vert1(:, 2), tri_vert1(:, 1) - tri_vert2(:, 1)];
                ineq_b_cell{tri_id} = tri_vert2(:, 2) .* tri_vert1(:, 1) - tri_vert1(:, 2) .* tri_vert2(:, 1);

                % correct numerical inaccuracies and make sure that the vertices themselves are contained in the triangular region
                ineq_b_cell{tri_id} = max(ineq_b_cell{tri_id}, max(ineq_w_cell{tri_id} * tri_vert1', [], 2));
            end
            
            obj.Dens.EdgeIneq.w = vertcat(ineq_w_cell{:});
            obj.Dens.EdgeIneq.b = vertcat(ineq_b_cell{:});

            % prepare information needed for rejection sampling
            obj.RejSamp = struct;
            obj.RejSamp.Range = struct;
            obj.RejSamp.Range.x = [min(vertices(:, 1)), max(vertices(:, 1))];
            obj.RejSamp.Range.y = [min(vertices(:, 2)), max(vertices(:, 2))];

            obj.RejSamp.Multiplier = max(dens_vertices) * (obj.RejSamp.Range.x(2) - obj.RejSamp.Range.x(1)) ...
                * (obj.RejSamp.Range.y(2) - obj.RejSamp.Range.y(1));
        end

        function dens = densityFunction(obj, pts)
            % Compute the probability density function
            % Input: 
            %   pts: two-column matrix containing the input points
            % Output: 
            %   dens: vector containing the computed densities

            input_num = size(pts, 1);

            % evaluate the density by interpolation in each triangle (ignoring whether the point is in the triangle)
            dens_interp = obj.Dens.InterpMat * [pts'; ones(1, input_num)];

            tri_inside_mat = obj.checkIfInsideTriangularRegion(pts);

            % normalize the matrix such that each column sums to 1; if a point belongs to more than one triangles, then the densities
            % calculated with respect to each of such triangles should anyways be the same
            tri_inside_mat_sum = sum(tri_inside_mat, 2);
            tri_inside_mat_sum(tri_inside_mat_sum == 0) = 1;
            tri_inside_mat = tri_inside_mat ./ tri_inside_mat_sum;

            dens = sum(dens_interp' .* tri_inside_mat, 2);
        end
    
        function inside_mat = checkIfInsideTriangularRegion(obj, pts)
            % Check if each of the input points is inside each of the triangular regions
            % Input: 
            %   pts: two-column matrix containing the inputs
            % Output:
            %   inside_mat: logical matrix where each row represents an input point and each column represents a triangle region

            region_num = size(obj.Supp.Triangles, 1);
            input_num = size(pts, 1);

            region_ineq_mat = obj.Dens.EdgeIneq.w * pts' - obj.Dens.EdgeIneq.b <= ProbMeas2DConvexPolytope.INSIDE_TOLERANCE;
            inside_mat = reshape(all(reshape(region_ineq_mat(:), 3, region_num * input_num), 1)', region_num, input_num)';
        end

        function pts_san = sanitizePoints(obj, pts, tolerance)
            % Fix small numerical inaccuracies in the points such that points that are slightly outside the support are moved into the 
            % support
            % Inputs: 
            %   pts: two-column matrix containing the input points
            %   tolerance: numerical tolerance values for constraint violations (default is 1e-6)
            % Output:
            %   pts_san: two-column matrix containing the sanitized points

            if ~exist('tolerance', 'var') || isempty(tolerance)
                tolerance = 1e-6;
            end

            region_num = size(obj.Supp.Triangles, 1);
            input_num = size(pts, 1);
            constr_vio = obj.Dens.EdgeIneq.w * pts' - obj.Dens.EdgeIneq.b;

            vio_mat = reshape(max(reshape(constr_vio(:), 3, region_num * input_num), [], 1)', region_num, input_num)';

            if max(min(vio_mat, [], 2)) > tolerance
                error('constraint violations exceed the tolerance');
            end

            pts_san = pts;

            for input_id = 1:input_num
                % find the triangle that the point is closest to
                [min_vio, min_id] = min(vio_mat(input_id, :));

                if min_vio <= 0
                    continue;
                end

                pt = pts(input_id, :)';
                tri = obj.Supp.Triangles(min_id, :);
                tri_verts = obj.Supp.TriangleVertices(tri, :);

                % represent the point as a convex combination of the vertices of the triangle
                coefs = [tri_verts'; ones(1, 3)] \ [pt; 1];

                % make all the coefficients positive
                coefs = max(coefs, 1e-6);

                % normalize the coefficients to sum to 1
                coefs = coefs / sum(coefs);

                % update the point
                pts_san(input_id, :) = tri_verts' * coefs;
            end
        end

        function info = saveW2OptimalTransportInfo(obj)
            % Save optimal transport related information in a struct that can be used later
            % Output:
            %   info: struct that copies the fields from obj.W2OT

            info = saveW2OptimalTransportInfo@ProbMeas2DConvexPolytopeWithW2OT(obj);
            info.condrejsamp = obj.W2OT.CondRejSamp;
        end

        function loadW2OptimalTransportInfo(obj, info)
            % Load optimal transport related information from a previously saved struct
            % Input:
            %   info: struct that will be copied to the fields of obj.W2OT

            loadW2OptimalTransportInfo@ProbMeas2DConvexPolytopeWithW2OT(obj, info);
            obj.W2OT.CondRejSamp = info.condrejsamp;
        end
    end

    methods(Access = protected)
        function inside = doCheckIfInsideSupport(obj, pts)
            % Check if the input points are inside the support of the probability measure
            % Input:
            %   pts: two-column matrix containing the input points
            % Output: 
            %   inside: boolean vector indicating whether each point is inside the support

            tri_num = size(obj.Supp.Triangles, 1);
            input_num = size(pts, 1);

            tri_ineq_mat = obj.Dens.EdgeIneq.w * pts' - obj.Dens.EdgeIneq.b <= ProbMeas2DConvexPolytope.INSIDE_TOLERANCE;
            tri_inside_mat = reshape(all(reshape(tri_ineq_mat(:), 3, tri_num * input_num), 1)', tri_num, input_num);
            inside = any(tri_inside_mat, 1)';
        end

        function raw_samps = sampleFromProposalDistribution(obj, samp_num, rand_stream)
            % Generate random samples from a proposal distribution for rejection sampling. For this probability measure, the proposal 
            % distribution is uniform on a square.
            % Inputs:
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling

            width_x = diff(obj.RejSamp.Range.x);
            width_y = diff(obj.RejSamp.Range.y);

            raw_samps = rand(rand_stream, samp_num, 2) .* [width_x, width_y] + [obj.RejSamp.Range.x(1), obj.RejSamp.Range.y(1)];
        end

        function accept_probs = rejSampAcceptProbFunction(obj, raw_samps)
            % Computes the acceptance probability when using rejection sampling to generate independent samples from the probability
            % measure
            % Input:
            %   raw_samps: two-column matrix containing the raw samples generated from the proposal distribution before rejection
            % Output:
            %   accept_probs: vector containing the computed acceptance probability

            accept_probs = obj.densityFunction(raw_samps) / (obj.RejSamp.Multiplier / (diff(obj.RejSamp.Range.x) ...
                * diff(obj.RejSamp.Range.y)));
        end

        function prep = prepareOptimalTransport(obj, atoms, probs, angle_num)
            % Prepare quantities for the computation of semi-discrete optimal transport.
            % Inputs:
            %   atoms: two-column matrix containing the locations of the atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in the discrete measure
            %   angle_num: number of angles used in the integration

            prep = prepareOptimalTransport@ProbMeas2DConvexPolytope(obj, atoms, probs, angle_num);

            angle_list = prep.angle_list;

            atom_num = size(atoms, 1);
            tri_num = size(obj.Supp.Triangles, 1);
            prep.radial = cell(atom_num, 1);

            % coefficients in the left-hand side of the inequalities characterizing the triangles with repect to the distance away from
            % the atom; 
            % if this value is positive, it means that the ray is "leaving" the closed half-space, i.e., when the radial distance is
            % large enough the points on the ray will eventually be outside the closed half-space; 
            % if this value is non-positive, it means that the ray is "entering" or "remaining in" the closed half-space, i.e., when
            % the radial distance is large enough the points on the ray will eventually be inside the closed half-space
            ineq_coef = obj.Dens.EdgeIneq.w * [cos(angle_list)'; sin(angle_list)'];

            for atom_id = 1:atom_num
                % differences between the left-hand side of the inequalities characterizing the triangles and the right-hand side 
                % evaluated at the atom; 
                % if this value is positive it means that the atom is outside the closed half-space;
                % if this value is non-positive it means that the atom is inside the closed half-space
                ineq_diff = obj.Dens.EdgeIneq.w * prep.atoms(atom_id, :)' - obj.Dens.EdgeIneq.b;
                
                % radial distance of the point of intersection of the line going through the atom and a boundary of a triangle
                ineq_isec = -ineq_diff ./ ineq_coef;

                % 0/0 may occur, which means that the atom is on the boundary of the closed half-space and the ray is parallel to the
                % boundary; in this case, we set the intersection point to be at infinity
                ineq_isec_nan = isnan(ineq_isec);
                ineq_isec(ineq_isec_nan) = inf;
                
                % lower boundaries of the intervals
                itvl_lb = zeros(tri_num * 3, angle_num);

                % upper boundaries of the intervals
                itvl_ub = inf(tri_num * 3, angle_num);

                % boolean matrix indicating whether each atom is outside each closed half-space
                ineq_diff_pos = ineq_diff > 0;

                % boolean matrix indicating whether each ray is "leaving" each closed half-space
                ineq_coef_pos = ineq_coef > 0;

                % boolean matrix indicating if an atom is inside a closed half-space and its ray is "leaving" the closed half-space
                mat_npos_pos = (~ineq_diff_pos) & ineq_coef_pos;

                % boolean matrix indicating if an atom is outside a closed half-space and its ray is "entering" the closed half-space
                mat_pos_npos = ineq_diff_pos & (~ineq_coef_pos);

                % boolean matrix indicating if an atom is outside a closed half-space and its ray is not "entering" the closed
                % half-space
                mat_pos_nneg = ineq_diff_pos & (ineq_coef >= 0);

                % if an atom is outside a closed half-space and its ray is "entering" the closed half-space, then the intersection of 
                % the ray and the closed half-space is a ray emitting from the point of intersection
                itvl_lb(mat_pos_npos) = ineq_isec(mat_pos_npos);

                % if an atom is inside a closed half-space and its ray is "leaving" the closed half-space, then the intersection of the
                % ray and the closed half-space is a line segment between the atom and the point of intersection
                itvl_ub(mat_npos_pos) = ineq_isec(mat_npos_pos);

                % if an atom is outside a closed half-space and its ray is not "entering" the closed half-space, then the intersection
                % of the ray and the closed half-space is empty
                itvl_ub(mat_pos_nneg) = 0;

                % taking the intersection of the three intervals to get the intersection of the ray and the triangle
                tri_itvl_lb = reshape(max(reshape(itvl_lb(:), 3, tri_num * angle_num), [], 1), tri_num, angle_num);
                tri_itvl_ub = reshape(min(reshape(itvl_ub(:), 3, tri_num * angle_num), [], 1), tri_num, angle_num);

                % for intervals that are empty (or have empty interior), we set both their lower bounds and upper bounds to 0
                tri_itvl_empty = tri_itvl_lb >= tri_itvl_ub;
                tri_itvl_lb(tri_itvl_empty) = 0;
                tri_itvl_ub(tri_itvl_empty) = 0;

                % the density function evaluated on this interval is an affine function: aff_coefficient * r + aff_intercept
                aff_intercept = obj.Dens.InterpMat * [prep.atoms(atom_id, :)'; 1];
                aff_coefficient = obj.Dens.InterpMat * [cos(angle_list'); sin(angle_list'); zeros(1, angle_num)];

                prep.radial{atom_id} = struct( ...
                    'interval_lb', tri_itvl_lb, ...
                    'interval_ub', tri_itvl_ub, ...
                    'intercept', aff_intercept, ...
                    'coefficient', aff_coefficient);
            end
        end

        function [val1, val2] = computeInnerIntegral(~, ...
                atom_id, ...
                prep, ...
                batch_indices, ...
                upper_limits)
            % Evaluates inner integrals along a list of angles in the polar coordinate system. The integrands are r -> r * p(r) and 
            % r -> r^2 * p(r) where p(r) denotes the probability density function.
            % Inputs:
            %   atom_id: the index of the atom (which is treated as the origin)
            %   prep: struct returned by the function obj.prepareOptimalTransport
            %   batch_indices: vector indicating the indices of the pre-specified angles where the inner integral is evaluated
            %   upper_limits: the upper integral limits specified in a vector (the lower integral limits is always 0)
            % Outputs:
            %   val1: the integral of r -> r * p(r)
            %   val2: the integral of r -> r^2 * p(r)

            radial = prep.radial{atom_id};
            coef = radial.coefficient(:, batch_indices);
            
            % intersect the intervals (representing line segments in every direction) with (0, upper_limits), and making sure that 
            % empty intervals are represented with equal lower and upper limits
            itg_ub = min(radial.interval_ub(:, batch_indices), upper_limits');
            itg_lb = min(radial.interval_lb(:, batch_indices), itg_ub);

            upper_limits_p2 = itg_ub .^ 2 - itg_lb .^ 2;
            upper_limits_p3 = itg_ub .^ 3 - itg_lb .^ 3;
            upper_limits_p4 = itg_ub .^ 4 - itg_lb .^ 4;

            val1 = sum(coef / 3 .* upper_limits_p3 + radial.intercept / 2 .* upper_limits_p2, 1)';
            val2 = sum(coef / 4 .* upper_limits_p4 + radial.intercept / 3 .* upper_limits_p3, 1)';
        end
     
        function postProcessOptimalTransport(obj, pp_angle_num)
            % Post-processing after computing semi-discrete optimal transport including computing the contour of each cell in the
            % Laguerre diagram and preparing for conditional rejection sampling
            % Input:
            %   pp_angle_num: number of angles used in post-processing

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);
            
            % cell array containing information about rejection sampling from the conditional distributions (conditional on each atom)
            obj.OT.CondRejSamp = cell(atom_num, 1);

            hyp_a = (obj.OT.Weights - obj.OT.Weights') / 2;

            % take an evenly-spaced 10000 angles for this step
            prep = obj.prepareOptimalTransport(obj.OT.DiscMeas.Atoms, obj.OT.DiscMeas.Probs, pp_angle_num);
            angle_samp = prep.angle_list;

            for atom_id = 1:atom_num
                [rho_temp, hyp_full] = ProbMeas2DConvexPolytope.hyperbolic_Dirichlet_tessellation_polar( ...
                    angle_samp, hyp_a(:, atom_id), ...
                    prep.hyp.c(:, atom_id), ...
                    prep.hyp.intercept(:, atom_id), ...
                    prep.lin{atom_id}.c, ...
                    prep.lin{atom_id}.intercept, ...
                    true);

                radial = prep.radial{atom_id};
                rad_coef = radial.coefficient;
                obj.OT.CondRejSamp{atom_id} = struct;
                itvl_ub_max = max(radial.interval_ub, [], 1)';

                % take the distance to be the smaller of the distance calculated from the Laguerre diagram and the intersection with 
                % the triangles
                rho_i = min(rho_temp, itvl_ub_max);

                % add 10% of the standard deviation to make sure that the circle covers the whole region
                std_rho = std(rho_i);
                radius = max(rho_i) + std_rho * 0.1;

                if radius <= 1e-12
                    error('unexpected error causing the radius to be non-positive');
                end

                obj.OT.CondRejSamp{atom_id}.Radius = radius;
                obj.OT.CondRejSamp{atom_id}.NeighborAtoms = [find(any(hyp_full < radius, 1))'; atom_id];
                
                % conservative estimate of the maximum density in the cell in each direction
                tri_density_lb = radial.interval_lb .* rad_coef;
                tri_density_ub = radial.interval_ub .* rad_coef;
                interval_nonempty = radial.interval_lb < radial.interval_ub;
                tri_max_density = (radial.intercept + max(tri_density_lb, tri_density_ub)) .* interval_nonempty;
                max_density = max(max(tri_max_density)) * 1.01 / obj.OT.DiscMeas.Probs(atom_id);

                obj.OT.CondRejSamp{atom_id}.MaxDensity = max_density;
                obj.OT.CondRejSamp{atom_id}.Multiplier = max_density * pi * radius ^ 2;
            end
        end

        function raw_samps = sampleFromCondProposalDistribution(obj, ...
                atom_id, ...
                samp_num, ...
                rand_stream)
            % Generate random samples from a proposal distribution for conditional rejection sampling. For this procedure, the proposal
            % distribution is uniform on a circle.
            % Inputs:
            %   atom_id: the index of the corresponding atom in the discrete measure
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling

            % randomly generate points uniformly in the circle
            dists_from_center = sqrt(rand(rand_stream, samp_num, 1)) * obj.OT.CondRejSamp{atom_id}.Radius;
            angles = rand(rand_stream, samp_num, 1) * 2 * pi;

            % raw samples before the rejection step
            raw_samps = obj.OT.DiscMeas.Atoms(atom_id, :) + dists_from_center .* [cos(angles), sin(angles)];
        end

        function accept_probs = condRejSampAcceptProbFunction(obj, ...
                atom_id, ...
                raw_samps)
            % Compute the acceptance probability in conditional rejection sampling from a cell
            % Inputs: 
            %   atom_id: the index of the atom (corresponding to the cell that is being sampled)
            %   raw_samps: two-column matrix containing the raw samples before the rejection step 
            % Output:
            %   accept_probs: vector containing the computed acceptance probabilities

            raw_density = obj.densityFunction(raw_samps);

            condrejsamp = obj.OT.CondRejSamp{atom_id};
            
            radius = condrejsamp.Radius;
            multiplier = condrejsamp.Multiplier;
            neighbors = condrejsamp.NeighborAtoms;

            neighbor_atoms = obj.OT.DiscMeas.Atoms(neighbors, :);

            weighted_dist = sqrt((raw_samps(:, 1) - neighbor_atoms(:, 1)') .^ 2 + (raw_samps(:, 2) - neighbor_atoms(:, 2)') .^ 2) ...
                - obj.OT.Weights(neighbors)';

            accept_probs = raw_density .* all(weighted_dist(:, neighbors == atom_id) <= weighted_dist, 2) ...
                / (obj.OT.DiscMeas.Probs(atom_id) * multiplier / (pi * radius ^ 2));
        end

        function [int1, int2] = computeW2InnerIntegral(obj, ...
                vertices, ...
                triangles, ...
                centers, ...
                mesh_indices)
            % Compute integrals of constant 1 functions as well as the squared distance functions (from given centers) within triangles
            % with respect to the probability measure
            % Inputs:
            %   vertices: two-column matrix containing the vertices in the triangulation
            %   triangles: three-column matrix encoding the triangulation
            %   centers: two-column matrix containing centers with respect to which the second integral is evaluated
            %   mesh_indices: vector containing the indices of the triangles in the triangulation of the support; used to speed up the 
            %   computation of the density
            % Outputs:
            %   int1: vector containing the computed integrals of the constant 1 function within each triangle; each component
            %   corresponds to a triangle
            %   int2: vector containing the computed integrals of the quadratic function (x, y) -> (x - c_x)^2 + (y - c_y)^2, where 
            %   (c_x, c_y) are the coordinates of the center; each component corresponds to a triangle

            triangle_num = size(triangles, 1);

            vertices_x = vertices(:, 1);
            vertices_y = vertices(:, 2);

            % the x- and y-coordinates of the three vertices of the triangles
            tri_vertices_x = vertices_x(triangles')';
            tri_vertices_y = vertices_y(triangles')';

            % since the density is piece-wise affine, the weights and intercepts depend on the mesh indices
            interp_mat = obj.Dens.InterpMat(mesh_indices, :);
            dens_vertices = tri_vertices_x .* interp_mat(:, 1) + tri_vertices_y .* interp_mat(:, 2) + interp_mat(:, 3);

            int_vals = ProbMeas2DAffDens.affineIntegralTriangle( ...
                repmat(dens_vertices, 2, 1), ...
                repmat(tri_vertices_x, 2, 1), ...
                repmat(tri_vertices_y, 2, 1), ...
                [repmat([0, 0, 0, 0, 0, 1], triangle_num, 1); ...
                [ones(triangle_num, 2), zeros(triangle_num, 1), -2 * centers, sum(centers.^ 2, 2)]], ...
                false(2 * triangle_num, 1));

            int1 = int_vals(1:triangle_num);
            int2 = int_vals(triangle_num + 1:end);
        end

        function postProcessW2OptimalTransport(obj)
            % Post-processing after computing semi-discrete Wasserstein-2 optimal transport including computing the contour of each
            % cell in the Laguerre diagram and preparing for conditional rejection sampling 

            atom_num = size(obj.W2OT.DiscMeas.Atoms, 1);
            mesh_triangle_num = size(obj.Supp.Triangles, 1);

            % shift the weights to make them all positive so that they can be interpreted as squared radii of circles in a power 
            % diagram
            squared_radii = obj.W2OT.Weights - min(obj.W2OT.Weights) + 1;

            % invoke the mex function to compute a triangular mesh formed by the intersections of power cells and 
            [output_vertices, ...
                output_triangulation, ...
                power_cell_indices, ...
                mesh_indices] ...
                = mesh_intersect_power_diagram( ...
                obj.Supp.TriangleVertices, ...
                obj.Supp.Triangles, ...
                obj.W2OT.DiscMeas.Atoms, ...
                squared_radii, ...
                obj.W2OT.BoundingBox.BottomLeft, ...
                obj.W2OT.BoundingBox.TopRight);

            triangle_num = size(output_triangulation, 1);

            vertices_x = output_vertices(:, 1);
            vertices_y = output_vertices(:, 2);

            % the x- and y-coordinates of the three vertices of the triangles
            tri_vertices_x = vertices_x(output_triangulation')';
            tri_vertices_y = vertices_y(output_triangulation')';

            % since the density is piece-wise affine, the weights and intercepts depend on the mesh indices
            interp_mat = obj.Dens.InterpMat(mesh_indices, :);
            dens_vertices = tri_vertices_x .* interp_mat(:, 1) + tri_vertices_y .* interp_mat(:, 2) + interp_mat(:, 3);

            % compute the probabilities contained in each triangle in the output mesh
            tri_probs = ProbMeas2DAffDens.affineIntegralTriangle( ...
                dens_vertices, ...
                tri_vertices_x, ...
                tri_vertices_y, ...
                repmat([0, 0, 0, 0, 0, 1], triangle_num, 1), ...
                false(triangle_num, 1));

            obj.W2OT.CondRejSamp = cell(atom_num, 1);

            for atom_id = 1:atom_num
                % first create as many sub-cells as the number of triangles in the density mesh
                subcell_probs = zeros(mesh_triangle_num, 1);
                subcells = cell(mesh_triangle_num, 1);

                % keep track of the number of triangles with non-empty intersection with the power cell
                mesh_triangle_counter = 0;

                for mesh_triangle_id = 1:mesh_triangle_num
                    tri_intersect_list = power_cell_indices == atom_id & mesh_indices == mesh_triangle_id;

                    if all(~tri_intersect_list)
                        continue;
                    end

                    mesh_triangle_counter = mesh_triangle_counter + 1;

                    % gather the triangles that are in the intersection of the power cell and the triangle and take their union (which 
                    % is convex) to form the sub-cell
                    intersect_indices = output_triangulation(tri_intersect_list, :);
                    intersect_vertices = output_vertices(intersect_indices(:), :);
                    ch = convhull(intersect_vertices, 'Simplify', true);
                    ps = polyshape(intersect_vertices(ch, :), 'Simplify', false, 'KeepCollinearPoints', false);

                    % compute the total probability within this sub-cell
                    subcell_probs(mesh_triangle_counter) = sum(tri_probs(tri_intersect_list));

                    info = struct;

                    % store information about the shape of this sub-cell
                    info.BoxBottomLeft = min(intersect_vertices, [], 1)';
                    info.BoxWidthHeight = max(intersect_vertices, [], 1)' - info.BoxBottomLeft;
                    info.PolyShape = ps;

                    % store information about the density inside the sub-cell
                    info.DensityWeightVec = obj.Dens.InterpMat(mesh_triangle_id, 1:2)';
                    info.DensityIntercept = obj.Dens.InterpMat(mesh_triangle_id, 3);

                    % store information about the multiplier in the rejection sampling
                    max_density = max(intersect_vertices * info.DensityWeightVec + info.DensityIntercept) * 1.01 ...
                        / subcell_probs(mesh_triangle_counter);

                    info.MaxDensity = max_density;
                    info.Multiplier = max_density * info.BoxWidthHeight(1) * info.BoxWidthHeight(2);

                    subcells{mesh_triangle_counter} = info;
                end

                obj.W2OT.CondRejSamp{atom_id} = struct;
                obj.W2OT.CondRejSamp{atom_id}.SubCells = subcells(1:mesh_triangle_counter, :);
                obj.W2OT.CondRejSamp{atom_id}.SubCellProbabilities = subcell_probs(1:mesh_triangle_counter);
            end
        end

        function samps = condRejSampFromW2Coupling(obj, atom_id, samp_num, rand_stream)
            % Perform rejection sampling of the conditional distribution given an atom in the discrete measure; i.e., perform rejection
            % sampling within a power cell.
            % Inputs:
            %   atom_id: atom index
            %   samp_num: number of samples required
            %   rand_stream: RandStream object used for sampling
            % Output:
            %   samps: two-column matrix containing the generated samples

            subcells = obj.W2OT.CondRejSamp{atom_id}.SubCells;
            subcell_probs = obj.W2OT.CondRejSamp{atom_id}.SubCellProbabilities;
            subcell_num = length(subcells);

            % first sample the sub-cell indices; since we are sampling conditional to the power cell, the sub-cell probabilities need 
            % to be normalized
            subcell_index_samps = randsample(rand_stream, subcell_num, samp_num, true, subcell_probs / sum(subcell_probs));

            % count the number of samples for each cell
            subcell_samp_num_list = accumarray(subcell_index_samps, ones(samp_num, 1), [subcell_num, 1]);

            samp_cell = cell(subcell_num, 1);
            
            for subcell_id = 1:subcell_num
                subcell_info = subcells{subcell_id};
                subcell_samp_num = subcell_samp_num_list(subcell_id);
                samp_cell{subcell_id} = zeros(subcell_samp_num, 2);

                % the number of samples generated in each iteration is equal to 150% of samp_num times the expected number of 
                % iterations for each sample, capped at 1e4
                batch_samp_num = min(ceil(subcell_samp_num * subcell_info.Multiplier * 1.5), 1e4);
    
                samp_counter = 0;
    
                while samp_counter < subcell_samp_num
                    % raw samples before the rejection step
                    subcell_raw_samps = subcell_info.BoxWidthHeight' .* rand(rand_stream, batch_samp_num, 2) ...
                        + subcell_info.BoxBottomLeft';
    
                    % compute the acceptance probabilities
                    raw_density = (subcell_raw_samps * subcell_info.DensityWeightVec + subcell_info.DensityIntercept) ...
                        .* isinterior(subcell_info.PolyShape, subcell_raw_samps);
    
                    accept_probs = raw_density / (subcell_probs(subcell_id) * subcell_info.Multiplier ...
                        / (subcell_info.BoxWidthHeight(1) * subcell_info.BoxWidthHeight(2)));
    
                    if any(accept_probs > 1)
                        warning('acceptance probability larger than 1 (%.4f) has been detected', max(accept_probs));
                    end
    
                    % the list of samples accepted
                    accepted = rand(rand_stream, batch_samp_num, 1) <= accept_probs;
    
                    acc_samps = subcell_raw_samps(accepted, :);
                    acc_num = size(acc_samps, 1);
    
                    % the number of new samples to be added to the list
                    new_samp_num = min(acc_num, subcell_samp_num - samp_counter);
    
                    samp_cell{subcell_id}(samp_counter + (1:new_samp_num), :) = acc_samps(1:new_samp_num, :);
    
                    % update the sample counter
                    samp_counter = samp_counter + new_samp_num;
                end
            end

            samps = vertcat(samp_cell{:});

            % order the samples according the sampled sub-cell indices
            [~, ordering] = sort(subcell_index_samps, 'ascend');
            samps(ordering, :) = samps;
        end
    
        function integrals = doSetSimplicialTestFuncs(obj, ...
                vertices, triangles_cell, check_simpcover)
            % Initialize the test functions with respect to a simplicial cover and compute the integrals of the test functions with
            % respect to the probability measure
            % Inputs: 
            %   vertices: two-column matrix containing the vertices used in the triangulation; the first rows must be identical to the
            %   vertices of the support
            %   triangles_cell: cell array where the number of cells is equal to the number of triangular regions; each cell contains 
            %   a three-column matrix representing the triangulation within a triangular region
            %   check_simpcover: check whether the given triangulation forms a simplicial cover (default is true)
            % Output:
            %   integrals: vector containing integrals of the test functions with respect to the probability measure; the number of 
            %   test functions is equal to the number of vertices in the triangulation

            if ~exist('check_simpcover', 'var') || isempty(check_simpcover)
                check_simpcover = true;
            end

            region_num = size(obj.Supp.Triangles, 1);
            assert(length(triangles_cell) == region_num, 'triangulation mis-specified');

            % check if the first vertices are identical to the vertices in the triangular regions
            assert(max(max(abs(vertices(1:size(obj.Supp.TriangleVertices, 1), :) - obj.Supp.TriangleVertices))) < 1e-14, ...
                'the vertices do not begin with the vertices of the triangular regions');

            % check if there are duplicate vertices
            assert(size(unique(round(vertices, 6), 'rows'), 1) == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not appear in any triangle)
            triangles_agg = vertcat(triangles_cell{:});
            assert(length(unique(triangles_agg(:))) == size(vertices, 1), 'there are redundant vertices');

            region_inside_mat = obj.checkIfInsideTriangularRegion(vertices);

            for region_id = 1:region_num
                % check if all vertices are contained in the triangular region
                region_vert_indices = unique(triangles_cell{region_id}(:));
                assert(all(region_inside_mat(region_vert_indices, region_id)), ...
                    'there are vertices used in the triangulation of this region that is outside of the region');
            end

            if check_simpcover
                % check if the triangulation forms a simplicial cover
                [empty, is_simpcover] = check_triangulation_simplicial_cover(vertices, triangles_agg);
    
                if any(empty)
                    error('some triangles are empty');
                end
    
                if ~is_simpcover
                    error('the triangulation does not form a simplicial cover');
                end
            end
            
            obj.SimplicialTestFuncs = struct;
            obj.SimplicialTestFuncs.Vertices = vertices;
            obj.SimplicialTestFuncs.Triangles = triangles_agg;

            vert_num = size(vertices, 1);
            tri_num = size(triangles_agg, 1);

            dens_vertices = obj.densityFunction(vertices);

            % evaluate the integral of three interpolation functions with respect to the probability measure within each triangular
            % region via an affine transformation which places the vertices on (0, 0), (1, 0), and (0, 1)
            tri_integrals = zeros(tri_num, 3);

            % compute the inverse transformation matrix which transform a coordinate into weights in the affine combination
            invtrans_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert = vertices(triangles_agg(tri_id, :), :);
                tri_dens = dens_vertices(triangles_agg(tri_id, :));
                trans_mat = tri_vert(2:3, :)' - tri_vert(1, :)';

                tri_integrals(tri_id, :) = (sum(tri_dens) + tri_dens) / 24 * abs(det(trans_mat));

                % compute the inverse transformation matrix
                invtrans_cell{tri_id} = [tri_vert'; ones(1, 3)] \ eye(3);
            end

            integrals = accumarray(triangles_agg(:), tri_integrals(:), [vert_num, 1]);
            obj.SimplicialTestFuncs.Integrals = integrals;

            if abs(sum(integrals) - 1) > 1e-10
                error('integrals are incorrectly computed');
            end

            obj.SimplicialTestFuncs.InvTransMat = vertcat(invtrans_cell{:});
        end
    end
end