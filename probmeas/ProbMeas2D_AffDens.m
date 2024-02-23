classdef ProbMeas2D_AffDens < ProbMeas2D_ConvexPolytope ...
        & HasTractableQuadraticIntegrals
    % Class for probability measures with affine density function supported
    % on a convex polytope

    methods(Access = public)
        function obj = ProbMeas2D_AffDens(vertices, dens_weight, ...
                dens_intercept)
            % Constructor function
            % Inputs:
            %   vertices: two-column matrix containing the vertices of the
            %   convex polytope which is the support of the probability
            %   measure
            %   dens_weight: weight vector in the affine density function
            %   dens_intercept: constant intercept in the affine density
            %   function

            % call the superclass constructor to initialize the support
            obj@ProbMeas2D_ConvexPolytope(vertices);

            if size(obj.Supp.Vertices, 1) < size(vertices, 1)
                error(['the support is not convex or there are ' ...
                    'redundant vertices']);
            end

            obj.Dens = struct;

            % evaluate the density at the vertices
            dens_vert = sum(dens_weight' .* vertices, 2) + dens_intercept;

            if any(dens_vert < 0)
                error('the density is not non-negative');
            end

            % divide the convex polytope into triangles via Delaunay
            % triangulation
            DT = delaunay(vertices);
            tri_num = size(DT, 1);

            % evaluate the integral of the unnormalized density within each
            % triangular region via an affine transformation which places
            % the vertices on (0, 0), (1, 0), and (0, 1)
            tri_integral = zeros(tri_num, 1);

            % evaluate the integrals of x1, x2, x1^2, x2^2, and x1*x2 with 
            % respect to the mixture of bivariate Gaussian measure within 
            % each triangular region
            tri_lin_integral_x1 = zeros(tri_num, 1);
            tri_lin_integral_x2 = zeros(tri_num, 1);
            tri_quad_integral_x1sq = zeros(tri_num, 1);
            tri_quad_integral_x2sq = zeros(tri_num, 1);
            tri_quad_integral_x1x2 = zeros(tri_num, 1);

            for tri_id = 1:tri_num
                tri_verts = vertices(DT(tri_id, :), :);
                dens_vertices = tri_verts * dens_weight + dens_intercept;
                
                int_vals = ...
                        ProbMeas2D_AffDens.affineIntegralTriangle( ...
                        dens_vertices, tri_verts', eye(6), false(6, 1));

                tri_integral(tri_id) = int_vals(6);
                tri_lin_integral_x1(tri_id) = int_vals(4);
                tri_lin_integral_x2(tri_id) = int_vals(5);
                tri_quad_integral_x1sq(tri_id) = int_vals(1);
                tri_quad_integral_x2sq(tri_id) = int_vals(2);
                tri_quad_integral_x1x2(tri_id) = int_vals(3);
            end

            % calculate the normalizing constant and normalize the density
            norm_const = sum(tri_integral);
            obj.Dens = struct('Weight', dens_weight / norm_const, ...
                'Intercept', dens_intercept / norm_const, ...
                'NormConst', norm_const);

            % calculate the expected value of x1, x2, x1^2, x2^2, x1*x2
            obj.FirstMomentVec = zeros(2, 1);
            obj.FirstMomentVec(1) = sum(tri_lin_integral_x1) / norm_const;
            obj.FirstMomentVec(2) = sum(tri_lin_integral_x2) / norm_const;
            obj.SecondMomentMat = zeros(2, 2);
            obj.SecondMomentMat(1, 1) = ...
                sum(tri_quad_integral_x1sq) / norm_const;
            obj.SecondMomentMat(2, 2) = ...
                sum(tri_quad_integral_x2sq) / norm_const;
            obj.SecondMomentMat(1, 2) = ...
                sum(tri_quad_integral_x1x2) / norm_const;
            obj.SecondMomentMat(2, 1) = obj.SecondMomentMat(1, 2);

            % prepare information needed for rejection sampling
            obj.RejSamp = struct;
            obj.RejSamp.Range = struct;
            obj.RejSamp.Range.x = [min(vertices(:, 1)), ...
                max(vertices(:, 1))];
            obj.RejSamp.Range.y = [min(vertices(:, 2)), ...
                max(vertices(:, 2))];

            obj.RejSamp.Multiplier = (max(sum(obj.Dens.Weight' ...
                .* vertices, 2)) + obj.Dens.Intercept) ...
                * (obj.RejSamp.Range.x(2) - obj.RejSamp.Range.x(1)) ...
                * (obj.RejSamp.Range.y(2) - obj.RejSamp.Range.y(1));
        end

        function dens = densityFunction(obj, pts)
            % Probability density function
            % Input:
            %   pts: two-column matrix containing the input points
            % Output:
            %   dens: the computed density values

            dens = all(pts * obj.Supp.Hyperplane.w' ...
                <= obj.Supp.Hyperplane.b', 2) ...
                .* (pts * obj.Dens.Weight + obj.Dens.Intercept);
        end

        function pts_san = sanitizePoints(obj, pts, tolerance)
            % Fix small numerical inaccuracies in the points such that
            % points that are slightly outside the support are moved into
            % the support
            % Inputs: 
            %   pts: two-column matrix containing the input points
            %   tolerance: numerical tolerance values for constraint
            %   violations (default is 1e-6)
            % Output:
            %   pts_san: two-column matrix containing the sanitized points

            if ~exist('tolerance', 'var') || isempty(tolerance)
                tolerance = 1e-6;
            end

            % partition the support into triangles
            supp_verts = obj.Supp.Vertices;
            supp_tri = delaunay(supp_verts);

            tri_num = size(supp_tri, 1);

            % calculate inequalities representing half-spaces bounding each
            % triangle
            ineq_w_cell = cell(tri_num, 1);
            ineq_b_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_verts = supp_tri(tri_id, :)';
                tri_vert_mat = supp_verts(tri_verts, :);

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
            
            ineq_w = vertcat(ineq_w_cell{:});
            ineq_b = vertcat(ineq_b_cell{:});

            input_num = size(pts, 1);
            constr_vio = ineq_w * pts' - ineq_b;

            vio_mat = reshape(max(reshape(constr_vio(:), ...
                3, tri_num * input_num), [], 1)', ...
                tri_num, input_num)';

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
                tri = supp_tri(min_id, :);
                tri_verts = supp_verts(tri, :);

                % represent the point as a convex combination of the
                % vertices of the triangle
                coefs = [tri_verts'; ones(1, 3)] \ [pt; 1];

                % make all the coefficients positive
                coefs = max(coefs, 1e-6);

                % normalize the coefficients to sum to 1
                coefs = coefs / sum(coefs);

                % update the point
                pts_san(input_id, :) = tri_verts' * coefs;
            end
        end
    end

    methods(Access = protected)
        function inside = doCheckIfInsideSupport(obj, pts)
            % Check if the input points are inside the support of the
            % probability measure. 
            % Input:
            %   pts: two-column matrix containing the input points
            % Output: 
            %   inside: boolean vector indicating whether each point is
            %   inside the support

            inside = all(pts * obj.Supp.Hyperplane.w' ...
                - obj.Supp.Hyperplane.b' <= ...
                ProbMeas2D_ConvexPolytope.INSIDE_TOLERANCE, 2);
        end

        function raw_samps = sampleFromProposalDistribution(obj, ...
                samp_num, rand_stream)
            % Generate random samples from a proposal distribution for
            % rejection sampling. For this probability measure, the
            % proposal distribution is uniform on a square.
            % Inputs:
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling
            % Output:
            %   raw_samps: two-column matrix containing the generated raw
            %   samples (before the rejection step)

            width_x = diff(obj.RejSamp.Range.x);
            width_y = diff(obj.RejSamp.Range.y);

            raw_samps = rand(rand_stream, samp_num, 2) ...
                .* [width_x, width_y] ...
                + [obj.RejSamp.Range.x(1), obj.RejSamp.Range.y(1)];
        end

        function accept_probs = rejSampAcceptProbFunction(obj, raw_samps)
            % Computes the acceptance probability when using rejection 
            % sampling to generate independent samples from the probability
            % measure
            % Input:
            %   raw_samps: two-column matrix containing the raw samples
            %   generated from the proposal distribution before rejection
            % Output:
            %   accept_probs: vector containing the computed acceptance
            %   probability

            accept_probs = obj.densityFunction(raw_samps) ...
                    / (obj.RejSamp.Multiplier ...
                    / (diff(obj.RejSamp.Range.x) ...
                    * diff(obj.RejSamp.Range.y)));
        end

        function prepareOptimalTransport(obj, atoms, probs, angle_list)
            % Prepare quantities for the computation of semi-discrete
            % optimal transport.
            % Inputs:
            %   atoms: two-column matrix containing the locations of the
            %   atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in
            %   the discrete measure
            %   angle_list: vector containing angles used in the 
            %   integration

            prepareOptimalTransport@ProbMeas2D_ConvexPolytope(obj, ...
                atoms, probs, angle_list);

            angle_list = obj.OT.Prep.angle_list;

            obj.OT.Prep.radial = struct;

            % intercept in the density function evaluated along a ray
            obj.OT.Prep.radial.intercept = atoms * obj.Dens.Weight ...
                + obj.Dens.Intercept;

            % coefficient in the density function evaluated along a ray
            obj.OT.Prep.radial.coefficient = ...
                [cos(angle_list), sin(angle_list)] * obj.Dens.Weight;
        end

        function [val1, val2] = computeInnerIntegral(obj, ...
                atom_id, angle_indices, upper_limits)
            % Evaluates inner integrals along a list of angles in the polar
            % coordinate system. The integrands are r -> r * p(r) and 
            % r -> r^2 * p(r) where p(r) denotes the probability density
            % function.
            % Inputs:
            %   atom_id: the index of the atom (which is treated as the
            %   origin)
            %   angle_indices: vector indicating the indices of the
            %   pre-specified angles where the inner integral is evaluated
            %   upper_limits: the upper integral limits specified in a
            %   vector (the lower integral limits is always 0)
            % Outputs:
            %   val1: the integral of r -> r * p(r)
            %   val2: the integral of r -> r^2 * p(r)

            coef = obj.OT.Prep.radial.coefficient(angle_indices);

            upper_limits_p2 = upper_limits .^ 2;
            upper_limits_p3 = upper_limits .^ 3;
            upper_limits_p4 = upper_limits .^ 4;

            val1 = coef / 3 .* upper_limits_p3 ...
                + obj.OT.Prep.radial.intercept(atom_id) / 2 ...
                * upper_limits_p2;
            val2 = coef / 4 .* upper_limits_p4 ...
                + obj.OT.Prep.radial.intercept(atom_id) / 3 ...
                * upper_limits_p3;
        end

        function postProcessOptimalTransport(obj, pp_angle_indices)
            % Post-processing after computing semi-discrete optimal
            % transport including computing the contour of each cell in the
            % Laguerre diagram and preparing for conditional rejection
            % sampling
            % Input:
            %   pp_angle_indices: indices of the angles used in
            %   post-processing (default is all indices)

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);

            % cell array containing information about rejection sampling
            % from the conditional distributions (conditional on each atom)
            obj.OT.CondRejSamp = cell(atom_num, 1);

            hyp_a = (obj.OT.Weights - obj.OT.Weights') / 2;

            angle_samp = obj.OT.Prep.angle_list(pp_angle_indices);

            for atom_id = 1:atom_num
                [rho_i, hyp_full] ...
                    = hyperbolic_Dirichlet_tessellation_polar( ...
                    angle_samp, hyp_a(:, atom_id), ...
                    obj.OT.Prep.hyp.c(:, atom_id), ...
                    obj.OT.Prep.hyp.intercept(:, atom_id), ...
                    obj.OT.Prep.lin{atom_id}.c, ...
                    obj.OT.Prep.lin{atom_id}.intercept, true);
                
                % prepare information about 
                obj.OT.CondRejSamp{atom_id} = struct;

                % add 10% of the standard deviation to make sure that the
                % circle covers the whole region
                std_rho = std(rho_i);
                radius = max(rho_i) + std_rho * 0.1;

                obj.OT.CondRejSamp{atom_id}.Radius = radius;
                obj.OT.CondRejSamp{atom_id}.NeighborAtoms ...
                    = [find(any(hyp_full < radius, 1))'; atom_id];
                
                % conservative estimate of the maximum density in the cell
                % in each direction
                max_density = max(max(obj.OT.Prep.radial.coefficient( ...
                    pp_angle_indices), 0) .* rho_i ...
                    + obj.OT.Prep.radial.intercept(atom_id)) * 1.01 ...
                    / obj.OT.DiscMeas.Probs(atom_id);

                obj.OT.CondRejSamp{atom_id}.MaxDensity = max_density;
                obj.OT.CondRejSamp{atom_id}.Multiplier = max_density ...
                    * pi * radius ^ 2;
            end
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

        function accept_probs = condRejSampAcceptProbFunction(obj, ...
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

            raw_density = obj.densityFunction(raw_samps);
            
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
                vertices, triangles, check_simpcover)
            % Initialize the test functions with respect to a simplicial
            % cover and compute the integrals of the test functions with
            % respect to the probability measure
            % Inputs: 
            %   vertices: two-column matrix containing the vertices used in
            %   the triangulation; the first rows must be identical to the
            %   vertices of the support
            %   triangles: three-column matrix containing the triangulation
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

            % check if the first vertices are identical to the support
            assert(max(max(abs(vertices(1:size(obj.Supp.Vertices, ...
                1), :) - obj.Supp.Vertices))) < 1e-14, ...
                ['the vertices do not begin with ' ...
                'the vertices of the support']);

            % check if all vertices are contained in the support
            assert(all(obj.checkIfInsideSupport(vertices)), ...
                'there are vertices outside the support');

            % check if there are duplicate vertices
            assert(size(unique(vertices, 'rows'), 1) ...
                == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not
            % appear in any triangle)
            assert(length(unique(triangles(:))) == size(vertices, 1), ...
                'there are redundant vertices');

            if check_simpcover
                % check if the triangulation forms a simplicial cover
                [empty, is_simpcover] ...
                    = check_triangulation_simplicial_cover(vertices, ...
                    triangles);
    
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
            obj.SimplicialTestFuncs.Triangles = triangles;

            vert_num = size(vertices, 1);
            tri_num = size(triangles, 1);

            dens_vertices = vertices * obj.Dens.Weight ...
                + obj.Dens.Intercept;

            % evaluate the integral of three interpolation functions with
            % respect to the probability measure within each triangular
            % region via an affine transformation which places the vertices
            % on (0, 0), (1, 0), and (0, 1)
            tri_integrals = zeros(tri_num, 3);

            % compute the inverse transformation matrix which transform a
            % coordinate into weights in the affine combination
            invtrans_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert = vertices(triangles(tri_id, :), :);
                tri_dens = dens_vertices(triangles(tri_id, :));
                trans_mat = tri_vert(2:3, :)' - tri_vert(1, :)';

                tri_integrals(tri_id, :) = (sum(tri_dens) + tri_dens) ...
                    / 24 * abs(det(trans_mat));

                % compute the inverse transformation matrix
                invtrans_cell{tri_id} = [tri_vert'; ones(1, 3)] ...
                    \ eye(3);
            end

            integrals = accumarray(triangles(:), tri_integrals(:), ...
                [vert_num, 1]);
            obj.SimplicialTestFuncs.Integrals = integrals;

            if abs(sum(integrals) - 1) > 1e-10
                error('integrals are incorrectly computed');
            end

            obj.SimplicialTestFuncs.InvTransMat ...
                = vertcat(invtrans_cell{:});
        end
    end

    methods(Static, Access = public)
        function vals = affineIntegralTriangle(dens_vertices, ...
                z_mat, polynomial_coefs, after_transform)
            % Integrate bivariate second-order polynomials inside a
            % triangle over a measure with positive affine density
            % Inputs:
            %   dens_vertices: the affine density function evaluated at the
            %   three vertices of the triangle
            %   z_mat: 2-by-3 matrix containing the three vertices of the
            %   triangle as columns
            %   polynomial_coefs: matrix with six rows where each column
            %   represents a polynomial integrand and each row represents a
            %   coefficient in the polynomial: c1 * x1^2 + c2 * x2^2 
            %   + c3 * x1 * x2 + c4 * x1 + c5 * x2 + c6
            %   after_transform: logical vector indicating whether the
            %   coefficients of each polynomial is applicable to the
            %   variables before transformation or after transforming to
            %   the triangle with vertices [0; 0], [1; 0], and [1; 1]
            % Output: 
            %   vals: vector containing the values of the integrals; each
            %   entry corresponds to a polynomial

            % compute the affine transformation that will transform the
            % triangle formed by the columns of z_mat to the triangled with
            % vertices [0; 0], [1; 0], and [1; 1]
            trans_all = [0, 1, 1; 0, 0, 1] / [z_mat; ones(1, 3)];
            trans_mat = trans_all(:, 1:2);
            trans_shift = trans_all(:, 3);
            trans_mat_inv = inv(trans_mat);

            vals = zeros(size(polynomial_coefs, 2), 1);

            for poly_id = 1:size(polynomial_coefs, 2)
                if ~after_transform(poly_id)
                    % if the coefficients apply to the variables before
                    % transformation, we need to compute the coefficients
                    % after transformation
                    coefs = polynomial_coefs(:, poly_id);
                    coefs_mat = [coefs(1), coefs(3)/2; ...
                        coefs(3)/2, coefs(2)];
                    coefs_vec = coefs(4:5);
                    coefs_const = coefs(6);

                    tcoefs_mat = trans_mat_inv' * coefs_mat ...
                        * trans_mat_inv;
                    tcoefs_vec = trans_mat_inv' * coefs_vec ...
                        - 2 * tcoefs_mat * trans_shift;
                    tcoefs_const = trans_shift' * tcoefs_mat ...
                        * trans_shift ...
                        - coefs_vec' * trans_mat_inv ...
                        * trans_shift + coefs_const; %#ok<MINV>
                    coefs = [tcoefs_mat(1, 1); tcoefs_mat(2, 2); ...
                        2 * tcoefs_mat(1, 2); tcoefs_vec; tcoefs_const];
                else
                    % otherwise, we can directly integrate over the
                    % polynomial with the given coefficients
                    coefs = polynomial_coefs(:, poly_id);
                end

                p1 = dens_vertices(1);
                d21 = dens_vertices(2) - dens_vertices(1);
                d32 = dens_vertices(3) - dens_vertices(2);
                vals(poly_id) = ...
                    ((coefs(1) * d21 ...
                    + (coefs(1) * d32 + coefs(3) * d21) / 2 ...
                    + (coefs(3) * d32 + coefs(2) * d21) / 3 ...
                    + (coefs(2) * d32) / 4) / 5 ...
                    + ((coefs(1) * p1 + coefs(4) * d21) ...
                    + (coefs(3) * p1 + coefs(4) * d32 ...
                    + coefs(5) * d21) / 2 ...
                    + (coefs(2) * p1 + coefs(5) * d32) / 3) / 4 ...
                    + ((coefs(4) * p1 + coefs(6) * d21) ...
                    + (coefs(5) * p1 + coefs(6) * d32) / 2) / 3 ...
                    + (coefs(6) * p1) / 2) ...
                    / abs(det(trans_mat));
            end
        end
    end
end

