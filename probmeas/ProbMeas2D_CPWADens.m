classdef ProbMeas2D_CPWADens < ProbMeas2D_ConvexPolytope
    % Class for probability measures with continuous piece-wise affine
    % (CPWA) density function. The support of the probability measure is
    % formed by the union of triangles. The triangles need to satisfy the
    % following condition: if two triangles A and B are not disjoint, then
    % the intersection of A and B needs to be a face of both A and B, i.e.,
    % the intersection of A and B is either an edge of A or a vertex of A,
    % and is also either an edge of B or a vertex of B. 

    methods(Access = public)
        function obj = ProbMeas2D_CPWADens(vertices, triangles, ...
                dens_vertices)
            % Constructor function
            % Inputs:
            %   vertices: two-column matrix containing the vertices of all
            %   triangles
            %   triangles: three-column matrix indicating the triangulation
            %   dens_triangles: vector containing the density at each
            %   vertex (thus specifying the entire density function via
            %   interpolation)

            % call the superclass constructor to initialize the support
            obj@ProbMeas2D_ConvexPolytope(vertices);

            obj.Dens = struct;

            if any(dens_vertices < 0)
                error('the density is not non-negative');
            end

            % check if there are duplicate vertices
            assert(size(unique(vertices, 'rows'), 1) ...
                == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not
            % appear in any triangle)
            assert(length(unique(triangles(:))) == size(vertices, 1), ...
                'there are redundant vertices');

            % check that the triangles are non-empty (have non-zero area)
            % and satisfy the condition
            [empty, is_simpcover] ...
                = check_triangulation_simplicial_cover( ...
                vertices, triangles);

            assert(all(~empty), 'there are empty triangles');
            assert(is_simpcover, ['the condition of simplicial cover ' ...
                'is violated']);

            obj.Supp.TriangleVertices = vertices;
            obj.Supp.Triangles = triangles;

            tri_num = size(triangles, 1);

            % evaluate the integral of the unnormalized density within each
            % triangular region
            tri_integral = zeros(tri_num, 1);
            for tri_id = 1:tri_num
                tri_vert_id = vertices(triangles(tri_id, :), :);
                trans_shift = tri_vert_id(1, :)';
                trans_mat = tri_vert_id(2:3, :)' - trans_shift;
                trans_density = dens_vertices(triangles(tri_id, :));

                tri_integral(tri_id) = (sum(trans_density) / 6) ...
                    * abs(det(trans_mat));
            end

            % calculate the normalizing constant and normalize the density
            norm_const = sum(tri_integral);
            dens_vertices = dens_vertices / norm_const;
            obj.Dens = struct('VertexDensities', dens_vertices, ...
                'NormConst', norm_const);

            obj.Dens.InterpMat = zeros(tri_num, 3);
            obj.Dens.EdgeIneq = struct;
            
            ineq_w_cell = cell(tri_num, 1);
            ineq_b_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert_id = triangles(tri_id, :)';

                tri_vert_mat = vertices(tri_vert_id, :);

                % calculate the affine interpolation function for each
                % triangle (ignoring the fact that interpolation only makes
                % sense within the triangle)
                obj.Dens.InterpMat(tri_id, :) = ...
                    [tri_vert_mat, ones(3, 1)] ...
                    \ dens_vertices(tri_vert_id);

                % sort the three vertices counterclockwise
                ccorder = convhull(tri_vert_mat(:, 1), ...
                    tri_vert_mat(:, 2), 'Simplify', true);

                % there are now 4 rows, the first and the last are
                % identical
                tri_vert_cc = tri_vert_mat(ccorder, :);

                tri_vert1 = tri_vert_cc(1:3, :);
                tri_vert2 = tri_vert_cc(2:4, :);

                ineq_w_cell{tri_id} = ...
                    [tri_vert2(:, 2) - tri_vert1(:, 2), ...
                     tri_vert1(:, 1) - tri_vert2(:, 1)];
                ineq_b_cell{tri_id} = ...
                    tri_vert2(:, 2) .* tri_vert1(:, 1) ...
                    - tri_vert1(:, 2) .* tri_vert2(:, 1);

                % correct numerical inaccuracies and make sure that the
                % vertices themselves are contained in the triangular
                % region
                ineq_b_cell{tri_id} = max(ineq_b_cell{tri_id}, ...
                    max(ineq_w_cell{tri_id} * tri_vert1', [], 2));
            end
            
            obj.Dens.EdgeIneq.w = vertcat(ineq_w_cell{:});
            obj.Dens.EdgeIneq.b = vertcat(ineq_b_cell{:});

            % density function with the support
            obj.Dens.Func = @obj.computeDensity;

            % prepare information needed for rejection sampling
            obj.RejSamp = struct;
            obj.RejSamp.Range = struct;
            obj.RejSamp.Range.x = [min(vertices(:, 1)), ...
                max(vertices(:, 1))];
            obj.RejSamp.Range.y = [min(vertices(:, 2)), ...
                max(vertices(:, 2))];

            obj.RejSamp.Multiplier = max(dens_vertices) ...
                * (obj.RejSamp.Range.x(2) - obj.RejSamp.Range.x(1)) ...
                * (obj.RejSamp.Range.y(2) - obj.RejSamp.Range.y(1));
            obj.RejSamp.AcceptProbFunc = @(raw_samps)( ...
                obj.Dens.Func(raw_samps) ...
                    / (obj.RejSamp.Multiplier ...
                    / (diff(obj.RejSamp.Range.x) ...
                    * diff(obj.RejSamp.Range.y))));
        end
    
        function inside_mat = checkIfInsideTriangularRegion(obj, pts)
            % Check if each of the input points is inside each of the
            % triangular regions
            % Input: 
            %   pts: two-column matrix containing the inputs
            % Output:
            %   inside_mat: logical matrix where each row represents an
            %   input point and each column represents a triangle region

            region_num = size(obj.Supp.Triangles, 1);
            input_num = size(pts, 1);

            region_ineq_mat = obj.Dens.EdgeIneq.w * pts' ...
                - obj.Dens.EdgeIneq.b <= 1e-14;
            inside_mat = reshape(all(reshape(region_ineq_mat(:), ...
                3, region_num * input_num), 1)', region_num, input_num)';
        end
    end

    methods(Access = protected)
        function inside = doCheckIfInsideSupport(obj, pts)
            % Check if the input points are inside the support of the
            % probability measure
            % Input:
            %   pts: two-column matrix containing the input points
            % Output: 
            %   inside: boolean vector indicating whether each point is
            %   inside the support

            tri_num = size(obj.Supp.Triangles, 1);
            input_num = size(pts, 1);

            tri_ineq_mat = obj.Dens.EdgeIneq.w * pts' ...
                - obj.Dens.EdgeIneq.b <= 1e-14;
            tri_inside_mat = reshape(all(reshape(tri_ineq_mat(:), 3, ...
                tri_num * input_num), 1)', tri_num, input_num);
            inside = any(tri_inside_mat, 1)';
        end

        function dens = computeDensity(obj, x)
            % Compute the probability density function
            % Input: 
            %   x: two-column matrix containing the input points
            % Output: 
            %   dens: vector containing the computed densities

            input_num = size(x, 1);

            % evaluate the density by interpolation in each triangle
            % (ignoring whether the point is in the triangle)
            dens_interp = obj.Dens.InterpMat * [x'; ones(1, input_num)];

            tri_inside_mat = obj.checkIfInsideTriangularRegion(x);

            % normalize the matrix such that each column sums to 1; if a
            % point belongs to more than one triangles, then the densities
            % calculated with respect to each of such triangles should
            % anyways be the same
            tri_inside_mat_sum = sum(tri_inside_mat, 2);
            tri_inside_mat_sum(tri_inside_mat_sum == 0) = 1;
            tri_inside_mat = tri_inside_mat ./ tri_inside_mat_sum;

            dens = sum(dens_interp' .* tri_inside_mat, 2);
        end

        function raw_samps = sampleFromProposalDistribution(obj, ...
                samp_num, rand_stream)
            % Generate random samples from a proposal distribution for
            % rejection sampling. For this probability measure, the
            % proposal distribution is uniform on a square.
            % Inputs:
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling

            width_x = diff(obj.RejSamp.Range.x);
            width_y = diff(obj.RejSamp.Range.y);

            raw_samps = rand(rand_stream, samp_num, 2) ...
                .* [width_x, width_y] ...
                + [obj.RejSamp.Range.x(1), obj.RejSamp.Range.y(1)];
        end

        function prepareOptimalTransport(obj, atoms, probs, angle_num)
            % Prepare quantities for the computation of semi-discrete
            % optimal transport.
            % Inputs:
            %   atoms: two-column matrix containing the locations of the
            %   atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in
            %   the discrete measure
            %   angle_num: number of angles used in the integration

            prepareOptimalTransport@ProbMeas2D_ConvexPolytope(obj, ...
                atoms, probs, angle_num);

            angle_list = obj.OT.Prep.angle_list;

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);
            tri_num = size(obj.Supp.Triangles, 1);
            obj.OT.Prep.radial = cell(atom_num, 1);

            % coefficients in the left-hand side of the inequalities
            % characterizing the triangles with repect to the distance away
            % from the atom
            ineq_coef = obj.Dens.EdgeIneq.w * [cos(angle_list)'; ...
                sin(angle_list)'];

            for atom_id = 1:atom_num
                % differences between the left-hand side of the
                % inequalities characterizing the triangles and the
                % right-hand side evaluated at the atom
                ineq_diff = obj.Dens.EdgeIneq.w ...
                    * obj.OT.DiscMeas.Atoms(atom_id, :)' ...
                    - obj.Dens.EdgeIneq.b;
                
                % intersections of lines going through the atom and the
                % boundaries of triangles
                ineq_isec = -ineq_diff ./ ineq_coef;
                
                % lower boundaries of the intervals
                itvl_lb = zeros(tri_num * 3, angle_num);

                % upper boundaries of the intervals
                itvl_ub = inf(tri_num * 3, angle_num);

                ineq_diff_nneg = ineq_diff >= 0;
                ineq_coef_nneg = ineq_coef >= 0;
                mat_neg_nneg = (~ineq_diff_nneg) & ineq_coef_nneg;
                mat_nneg_neg = ineq_diff_nneg & (~ineq_coef_nneg);
                mat_nneg_nneg = ineq_diff_nneg & ineq_coef_nneg;

                itvl_lb(mat_nneg_neg) = ineq_isec(mat_nneg_neg);
                itvl_ub(mat_neg_nneg) = ineq_isec(mat_neg_nneg);
                itvl_ub(mat_nneg_nneg) = 0;

                % taking the intersection of the three intervals
                tri_itvl_lb = reshape(max(reshape(itvl_lb(:), 3, ...
                    tri_num * angle_num), [], 1), tri_num, angle_num);
                tri_itvl_ub = reshape(min(reshape(itvl_ub(:), 3, ...
                    tri_num * angle_num), [], 1), tri_num, angle_num);

                % set the lower bound of the interval to be at most equal
                % to the upper bound
                tri_itvl_lb = min(tri_itvl_lb, tri_itvl_ub);

                aff_intercept = obj.Dens.InterpMat ...
                    * [obj.OT.DiscMeas.Atoms(atom_id, :)'; 1];
                aff_coefficient = obj.Dens.InterpMat ...
                    * [cos(angle_list'); sin(angle_list'); ...
                    zeros(1, angle_num)];

                obj.OT.Prep.radial{atom_id} = struct( ...
                    'interval_lb', tri_itvl_lb, ...
                    'interval_ub', tri_itvl_ub, ...
                    'intercept', aff_intercept, ...
                    'coefficient', aff_coefficient);
            end
        end

        function [val1, val2] = computeInnerIntegral(obj, ...
                atom_id, upper_limits)
            % Evaluates inner integrals along a list of angles in the polar
            % coordinate system. The integrands are r -> r * p(r) and 
            % r -> r^2 * p(r) where p(r) denotes the probability density
            % function.
            % Inputs:
            %   atom_id: the index of the atom (which is treated as the
            %   origin)
            %   upper_limits: the upper integral limits specified in a
            %   vector (the lower integral limits is always 0)
            % Outputs:
            %   val1: the integral of r -> r * p(r)
            %   val2: the integral of r -> r^2 * p(r)

            radial = obj.OT.Prep.radial{atom_id};
            
            % intersect the intervals (representing line segments in every
            % direction) with (0, upper_limits), and making sure that empty
            % intervals are represented with equal lower and upper limits
            itg_ub = min(radial.interval_ub, upper_limits');
            itg_lb = min(radial.interval_lb, itg_ub);

            upper_limits_p2 = itg_ub .^ 2 - itg_lb .^ 2;
            upper_limits_p3 = itg_ub .^ 3 - itg_lb .^ 3;
            upper_limits_p4 = itg_ub .^ 4 - itg_lb .^ 4;

            val1 = sum(radial.coefficient / 3 .* upper_limits_p3 ...
                + radial.intercept / 2 .* upper_limits_p2, 1)';
            val2 = sum(radial.coefficient / 4 .* upper_limits_p4 ...
                + radial.intercept / 3 .* upper_limits_p3, 1)';
        end
     
        function postProcessOptimalTransport(obj)
            % Post-processing after computing semi-discrete optimal
            % transport including computing the contour of each cell in the
            % Laguerre diagram and preparing for conditional rejection
            % sampling

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);
            angle_num = length(obj.OT.Prep.angle_list);

            obj.OT.Laguerre = struct;
            
            % cell array containing information about rejection sampling
            % from the conditional distributions (conditional on each atom)
            obj.OT.CondRejSamp = cell(atom_num, 1);

            rho = zeros(angle_num, atom_num);

            hyp_a = (obj.OT.Weights - obj.OT.Weights') / 2;

            for atom_id = 1:atom_num
                [rho(:, atom_id), hyp_full] ...
                    = hyperbolic_Dirichlet_tessellation_polar( ...
                    obj.OT.Prep.angle_list, ...
                    hyp_a(:, atom_id), ...
                    obj.OT.Prep.hyp.c(:, atom_id), ...
                    obj.OT.Prep.hyp.intercept(:, atom_id), ...
                    obj.OT.Prep.lin{atom_id}.c, ...
                    obj.OT.Prep.lin{atom_id}.intercept, true);

                radial = obj.OT.Prep.radial{atom_id};
                obj.OT.CondRejSamp{atom_id} = struct;
                itvl_ub_max = max(radial.interval_ub, [], 1)';

                % take the distance to be the smaller of the distance
                % calculated from the Laguerre diagram and the intersection
                % with the triangles
                rho_i = min(rho(:, atom_id), itvl_ub_max);

                % add 10% of the standard deviation to make sure that the
                % circle covers the whole region
                std_rho = std(rho_i);
                radius = max(rho_i) + std_rho * 0.1;

                obj.OT.CondRejSamp{atom_id}.Radius = radius;
                obj.OT.CondRejSamp{atom_id}.NeighborAtoms ...
                    = [find(any(hyp_full < radius, 1))'; atom_id];
                
                % conservative estimate of the maximum density in the cell
                % in each direction
                tri_density_lb = radial.interval_lb .* radial.coefficient;
                tri_density_ub = radial.interval_ub .* radial.coefficient;
                interval_nonempty = radial.interval_lb ...
                    < radial.interval_ub;
                tri_max_density = (radial.intercept ...
                    + max(tri_density_lb, tri_density_ub)) ...
                    .* interval_nonempty;
                max_density = max(max(tri_max_density)) * 1.01 ...
                    / obj.OT.DiscMeas.Probs(atom_id);

                obj.OT.CondRejSamp{atom_id}.MaxDensity = max_density;
                obj.OT.CondRejSamp{atom_id}.Multiplier = max_density ...
                    * pi * radius ^ 2;
                obj.OT.CondRejSamp{atom_id}.AcceptProbFunc = ...
                    @(raw_samps)(obj.computeCondAcceptProb(atom_id, ...
                    raw_samps));
            end

            obj.OT.Laguerre.rho = rho;
            obj.OT.Laguerre.x = obj.OT.DiscMeas.Atoms(:, 1)' ...
                + rho .* cos(obj.OT.Prep.angle_list);
            obj.OT.Laguerre.y = obj.OT.DiscMeas.Atoms(:, 2)' ...
                + rho .* sin(obj.OT.Prep.angle_list);
        end

        function raw_samps = sampleFromCondProposalDistribution(obj, ...
                atom_id, samp_num, rand_stream)
            % Generate random samples from a proposal distribution for
            % conditional rejection sampling. For this procedure, the
            % proposal distribution is uniform on a circle.
            % Inputs:
            %   atom_id: the index of the corresponding atom in the
            %   discrete measure
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling

            % randomly generate points uniformly in the circle
            dists_from_center = sqrt(rand(rand_stream, samp_num, 1)) ...
                * obj.OT.CondRejSamp{atom_id}.Radius;
            angles = rand(rand_stream, samp_num, 1) * 2 * pi;

            % raw samples before the rejection step
            raw_samps = obj.OT.DiscMeas.Atoms(atom_id, :) ...
                + dists_from_center .* [cos(angles), sin(angles)];
        end

        function accept_probs = computeCondAcceptProb(obj, ...
                atom_id, raw_samps)
            % Compute the acceptance probability in conditional rejection
            % sampling from a cell
            % Inputs: 
            %   atom_id: the index of the atom (corresponding to the cell
            %   that is being sampled)
            %   raw_samps: two-column matrix containing the raw samples
            %   before the rejection step 
            % Output:
            %   accept_probs: vector containing the computed acceptance
            %   probabilities

            raw_density = obj.Dens.Func(raw_samps);

            condrejsamp = obj.OT.CondRejSamp{atom_id};
            
            radius = condrejsamp.Radius;
            multiplier = condrejsamp.Multiplier;
            neighbors = condrejsamp.NeighborAtoms;

            neighbor_atoms = obj.OT.DiscMeas.Atoms(neighbors, :);

            weighted_dist = sqrt((raw_samps(:, 1) ...
                - neighbor_atoms(:, 1)') .^ 2 ...
                + (raw_samps(:, 2) - neighbor_atoms(:, 2)') .^ 2) ...
                - obj.OT.Weights(neighbors)';

            accept_probs = raw_density .* all(weighted_dist(:, ...
                neighbors == atom_id) <= weighted_dist, 2) ...
                / (obj.OT.DiscMeas.Probs(atom_id) * multiplier ...
                / (pi * radius ^ 2));
        end
    
        function integrals = doSetSimplicialTestFuncs(obj, ...
                vertices, triangles_cell, check_simpcover)
            % Initialize the test functions with respect to a simplicial
            % cover and compute the integrals of the test functions with
            % respect to the probability measure
            % Inputs: 
            %   vertices: two-column matrix containing the vertices used in
            %   the triangulation; the first rows must be identical to the
            %   vertices of the support
            %   triangles_cell: cell array where the number of cells is
            %   equal to the number of triangular regions; each cell
            %   contains a three-column matrix representing the 
            %   triangulation within a triangular region
            %   check_simpcover: check whether the given triangulation
            %   forms a simplicial cover (default is true)
            % Output:
            %   integrals: vector containing integrals of the test
            %   functions with respect to the probability measure; the
            %   number of test functions is equal to the number of vertices
            %   in the triangulation

            if ~exist('check_simpcover', 'var') || isempty(check_simpcover)
                check_simpcover = true;
            end

            region_num = size(obj.Supp.Triangles, 1);
            assert(length(triangles_cell) == region_num, ...
                'triangulation mis-specified');

            % check if the first vertices are identical to the vertices in
            % the triangular regions
            assert(max(max(abs(vertices(1:size( ...
                obj.Supp.TriangleVertices, 1), :) ...
                - obj.Supp.TriangleVertices))) < 1e-14, ...
                ['the vertices do not begin with ' ...
                'the vertices of the triangular regions']);

            % check if there are duplicate vertices
            assert(size(unique(vertices, 'rows'), 1) ...
                == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not
            % appear in any triangle)
            triangles_agg = vertcat(triangles_cell{:});
            assert(length(unique(triangles_agg(:))) ...
                == size(vertices, 1), 'there are redundant vertices');

            region_inside_mat = obj.checkIfInsideTriangularRegion( ...
                vertices);

            for region_id = 1:region_num
                % check if all vertices are contained in the triangular
                % region
                region_vert_indices = unique(triangles_cell{region_id}(:));
                assert(all(region_inside_mat(region_vert_indices, ...
                    region_id)), ['there are vertices used in the ' ...
                    'triangulation of this region that is outside of ' ...
                    'the region']);
            end

            if check_simpcover
                % check if the triangulation forms a simplicial cover
                [empty, is_simpcover] ...
                    = check_triangulation_simplicial_cover(vertices, ...
                    triangles_agg);
    
                if any(empty)
                    error('some triangles are empty');
                end
    
                if ~is_simpcover
                    error(['the triangulation does not form ' ...
                        'a simplicial cover']);
                end
            end
            
            obj.SimplicialTestFuncs = struct;
            obj.SimplicialTestFuncs.Vertices = vertices;
            obj.SimplicialTestFuncs.Triangles = triangles_agg;

            vert_num = size(vertices, 1);
            tri_num = size(triangles_agg, 1);

            dens_vertices = obj.Dens.Func(vertices);

            % evaluate the integral of three interpolation functions with
            % respect to the probability measure within each triangular
            % region via an affine transformation which places the vertices
            % on (0, 0), (1, 0), and (0, 1)
            tri_integrals = zeros(tri_num, 3);

            % compute the inverse transformation matrix which transform a
            % coordinate into weights in the affine combination
            invtrans_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert = vertices(triangles_agg(tri_id, :), :);
                tri_dens = dens_vertices(triangles_agg(tri_id, :));
                trans_mat = tri_vert(2:3, :)' - tri_vert(1, :)';

                tri_integrals(tri_id, :) = (sum(tri_dens) + tri_dens) ...
                    / 24 * abs(det(trans_mat));

                % compute the inverse transformation matrix
                invtrans_cell{tri_id} = [tri_vert'; ones(1, 3)] ...
                    \ eye(3);
            end

            integrals = accumarray(triangles_agg(:), tri_integrals(:), ...
                [vert_num, 1]);
            obj.SimplicialTestFuncs.Integrals = integrals;

            if abs(sum(integrals) - 1) > 1e-10
                error('integrals are incorrectly computed');
            end


            obj.SimplicialTestFuncs.InvTransMat ...
                = vertcat(invtrans_cell{:});
        end
    end
end