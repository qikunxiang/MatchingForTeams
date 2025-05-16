classdef ProbMeas2DAffDens < ProbMeas2DConvexPolytopeWithW2OT ...
        & HasTractableQuadraticIntegrals
    % Class for probability measures with affine density function supported on a convex polytope

    methods(Access = public)
        function obj = ProbMeas2DAffDens(vertices, dens_weight, dens_intercept)
            % Constructor function
            % Inputs:
            %   vertices: two-column matrix containing the vertices of the convex polytope which is the support of the probability
            %   measure
            %   dens_weight: weight vector in the affine density function
            %   dens_intercept: constant intercept in the affine density function

            % call the superclass constructor to initialize the support
            obj@ProbMeas2DConvexPolytopeWithW2OT(vertices);

            if size(obj.Supp.Vertices, 1) < size(vertices, 1)
                error('the support is not convex or there are redundant vertices');
            end

            obj.Dens = struct;

            % evaluate the density at the vertices
            dens_vert = sum(dens_weight' .* vertices, 2) + dens_intercept;

            if any(dens_vert < 0)
                error('the density is not non-negative');
            end

            tri_num = size(obj.Supp.Triangles, 1);

            % evaluate the integral of the unnormalized density within each triangular region via an affine transformation which places
            % the vertices on (0, 0), (1, 0), and (0, 1)
            tri_integral = zeros(tri_num, 1);

            % evaluate the integrals of x1, x2, x1^2, x2^2, and x1*x2 with respect to the mixture of bivariate Gaussian measure within 
            % each triangular region
            tri_lin_integral_x1 = zeros(tri_num, 1);
            tri_lin_integral_x2 = zeros(tri_num, 1);
            tri_quad_integral_x1sq = zeros(tri_num, 1);
            tri_quad_integral_x2sq = zeros(tri_num, 1);
            tri_quad_integral_x1x2 = zeros(tri_num, 1);

            for tri_id = 1:tri_num
                tri_verts = obj.Supp.TriangleVertices(obj.Supp.Triangles(tri_id, :), :);
                dens_vertices = tri_verts * dens_weight + dens_intercept;

                int_vals = ProbMeas2DAffDens.affineIntegralTriangle( ...
                        repmat(dens_vertices', 6, 1), ...
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
            obj.Dens = struct('Weight', dens_weight / norm_const, ...
                'Intercept', dens_intercept / norm_const, ...
                'NormConst', norm_const);

            % calculate the expected value of x1, x2, x1^2, x2^2, x1*x2
            obj.FirstMomentVec = zeros(2, 1);
            obj.FirstMomentVec(1) = sum(tri_lin_integral_x1) / norm_const;
            obj.FirstMomentVec(2) = sum(tri_lin_integral_x2) / norm_const;
            obj.SecondMomentMat = zeros(2, 2);
            obj.SecondMomentMat(1, 1) = sum(tri_quad_integral_x1sq) / norm_const;
            obj.SecondMomentMat(2, 2) = sum(tri_quad_integral_x2sq) / norm_const;
            obj.SecondMomentMat(1, 2) = sum(tri_quad_integral_x1x2) / norm_const;
            obj.SecondMomentMat(2, 1) = obj.SecondMomentMat(1, 2);

            % prepare information needed for rejection sampling
            obj.RejSamp = struct;
            obj.RejSamp.Range = struct;
            obj.RejSamp.Range.x = [min(vertices(:, 1)), max(vertices(:, 1))];
            obj.RejSamp.Range.y = [min(vertices(:, 2)), max(vertices(:, 2))];

            obj.RejSamp.Multiplier = (max(sum(obj.Dens.Weight' .* vertices, 2)) + obj.Dens.Intercept) ...
                * (obj.RejSamp.Range.x(2) - obj.RejSamp.Range.x(1)) * (obj.RejSamp.Range.y(2) - obj.RejSamp.Range.y(1));
        end

        function dens = densityFunction(obj, pts)
            % Probability density function
            % Input:
            %   pts: two-column matrix containing the input points
            % Output:
            %   dens: the computed density values

            dens = all(pts * obj.Supp.Hyperplane.w' <= obj.Supp.Hyperplane.b', 2) .* (pts * obj.Dens.Weight + obj.Dens.Intercept);
        end

        function pts_san = sanitizePoints(obj, pts, tolerance)
            % Fix small numerical inaccuracies in the points such that points that are slightly outside the support are moved into
            % the support
            % Inputs: 
            %   pts: two-column matrix containing the input points
            %   tolerance: numerical tolerance values for constraint violations (default is 1e-6)
            % Output:
            %   pts_san: two-column matrix containing the sanitized points

            if ~exist('tolerance', 'var') || isempty(tolerance)
                tolerance = 1e-6;
            end

            % partition the support into triangles
            supp_verts = obj.Supp.Vertices;
            supp_tri = delaunay(supp_verts);

            tri_num = size(supp_tri, 1);

            % calculate inequalities representing half-spaces bounding each triangle
            ineq_w_cell = cell(tri_num, 1);
            ineq_b_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_verts = supp_tri(tri_id, :)';
                tri_vert_mat = supp_verts(tri_verts, :);

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
            
            ineq_w = vertcat(ineq_w_cell{:});
            ineq_b = vertcat(ineq_b_cell{:});

            input_num = size(pts, 1);
            constr_vio = ineq_w * pts' - ineq_b;

            vio_mat = reshape(max(reshape(constr_vio(:), 3, tri_num * input_num), [], 1)', tri_num, input_num)';

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
            % Check if the input points are inside the support of the probability measure. 
            % Input:
            %   pts: two-column matrix containing the input points
            % Output: 
            %   inside: boolean vector indicating whether each point is inside the support

            inside = all(pts * obj.Supp.Hyperplane.w' - obj.Supp.Hyperplane.b' <= ProbMeas2DConvexPolytope.INSIDE_TOLERANCE, 2);
        end

        function raw_samps = sampleFromProposalDistribution(obj, ...
                samp_num, rand_stream)
            % Generate random samples from a proposal distribution for rejection sampling. For this probability measure, the proposal 
            % distribution is uniform on a square.
            % Inputs:
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling
            % Output:
            %   raw_samps: two-column matrix containing the generated raw samples (before the rejection step)

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

        function prepareOptimalTransport(obj, atoms, probs, angle_list)
            % Prepare quantities for the computation of semi-discrete optimal transport.
            % Inputs:
            %   atoms: two-column matrix containing the locations of the atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in the discrete measure
            %   angle_list: vector containing angles used in the integration

            prepareOptimalTransport@ProbMeas2DConvexPolytope(obj, atoms, probs, angle_list);

            angle_list = obj.OT.Prep.angle_list;

            obj.OT.Prep.radial = struct;

            % intercept in the density function evaluated along a ray
            obj.OT.Prep.radial.intercept = atoms * obj.Dens.Weight + obj.Dens.Intercept;

            % coefficient in the density function evaluated along a ray
            obj.OT.Prep.radial.coefficient = [cos(angle_list), sin(angle_list)] * obj.Dens.Weight;
        end

        function [val1, val2] = computeInnerIntegral(obj, ...
                atom_id, angle_indices, upper_limits)
            % Evaluates inner integrals along a list of angles in the polar coordinate system. The integrands are r -> r * p(r) and 
            % r -> r^2 * p(r) where p(r) denotes the probability density function.
            % Inputs:
            %   atom_id: the index of the atom (which is treated as the origin)
            %   angle_indices: vector indicating the indices of the pre-specified angles where the inner integral is evaluated
            %   upper_limits: the upper integral limits specified in a vector (the lower integral limits is always 0)
            % Outputs:
            %   val1: the integral of r -> r * p(r)
            %   val2: the integral of r -> r^2 * p(r)

            coef = obj.OT.Prep.radial.coefficient(angle_indices);

            upper_limits_p2 = upper_limits .^ 2;
            upper_limits_p3 = upper_limits .^ 3;
            upper_limits_p4 = upper_limits .^ 4;

            val1 = coef / 3 .* upper_limits_p3 + obj.OT.Prep.radial.intercept(atom_id) / 2 * upper_limits_p2;
            val2 = coef / 4 .* upper_limits_p4 + obj.OT.Prep.radial.intercept(atom_id) / 3 * upper_limits_p3;
        end

        function postProcessOptimalTransport(obj, pp_angle_indices)
            % Post-processing after computing semi-discrete optimal transport including computing the contour of each cell in the
            % Laguerre diagram and preparing for conditional rejection sampling
            % Input:
            %   pp_angle_indices: indices of the angles used in post-processing (default is all indices)

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);

            % cell array containing information about rejection sampling from the conditional distributions (conditional on each atom)
            obj.OT.CondRejSamp = cell(atom_num, 1);

            hyp_a = (obj.OT.Weights - obj.OT.Weights') / 2;

            angle_samp = obj.OT.Prep.angle_list(pp_angle_indices);

            for atom_id = 1:atom_num
                [rho_i, hyp_full] = ProbMeas2DConvexPolytope ...
                    .hyperbolic_Dirichlet_tessellation_polar( ...
                    angle_samp, ...
                    hyp_a(:, atom_id), ...
                    obj.OT.Prep.hyp.c(:, atom_id), ...
                    obj.OT.Prep.hyp.intercept(:, atom_id), ...
                    obj.OT.Prep.lin{atom_id}.c, ...
                    obj.OT.Prep.lin{atom_id}.intercept, ...
                    true);
                
                % prepare information about 
                obj.OT.CondRejSamp{atom_id} = struct;

                % add 10% of the standard deviation to make sure that the circle covers the whole region
                std_rho = std(rho_i);
                radius = max(rho_i) + std_rho * 0.1;

                obj.OT.CondRejSamp{atom_id}.Radius = radius;
                obj.OT.CondRejSamp{atom_id}.NeighborAtoms = [find(any(hyp_full < radius, 1))'; atom_id];
                
                % conservative estimate of the maximum density in the cell in each direction
                max_density = max(max(obj.OT.Prep.radial.coefficient(pp_angle_indices), 0) .* rho_i ...
                    + obj.OT.Prep.radial.intercept(atom_id)) * 1.01 / obj.OT.DiscMeas.Probs(atom_id);

                obj.OT.CondRejSamp{atom_id}.MaxDensity = max_density;
                obj.OT.CondRejSamp{atom_id}.Multiplier = max_density * pi * radius ^ 2;
            end
        end

        function raw_samps = sampleFromCondProposalDistribution(obj, ...
                atom_id, samp_num, rand_stream)
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
                atom_id, raw_samps)
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

            weighted_dist = sqrt((raw_samps(:, 1) - neighbor_atoms(:, 1)') .^ 2 ...
                + (raw_samps(:, 2) - neighbor_atoms(:, 2)') .^ 2) - obj.OT.Weights(neighbors)';

            accept_probs = raw_density .* all(weighted_dist(:, neighbors == atom_id) <= weighted_dist, 2) ...
                / (obj.OT.DiscMeas.Probs(atom_id) * multiplier / (pi * radius ^ 2));
        end

        function [int1, int2] = computeW2InnerIntegral(obj, ...
                vertices, triangles, centers, ~)
            % Compute integrals of constant 1 functions as well as the squared distance functions (from given centers) within triangles
            % with respect to the probability measure
            % Inputs:
            %   vertices: two-column matrix containing the vertices in the triangulation
            %   triangles: three-column matrix encoding the triangulation
            %   centers: two-column matrix containing centers with respect to which the second integral is evaluated
            %   mesh_indices: vector containing the indices of the triangles in the triangulation of the support; not useful here since 
            %   the density is globally affine
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

            % since the density is affine on the entire support, it can be directly evaluated using the x- and y-coordinates
            dens_vertices = tri_vertices_x * obj.Dens.Weight(1) + tri_vertices_y * obj.Dens.Weight(2) + obj.Dens.Intercept;

            int_vals = ProbMeas2DAffDens.affineIntegralTriangle( ...
                repmat(dens_vertices, 2, 1), ...
                repmat(tri_vertices_x, 2, 1), ...
                repmat(tri_vertices_y, 2, 1), ...
                [repmat([0, 0, 0, 0, 0, 1], triangle_num, 1); ...
                [ones(triangle_num, 2), zeros(triangle_num, 1), ...
                -2 * centers, sum(centers.^ 2, 2)]], ...
                false(2 * triangle_num, 1));

            int1 = int_vals(1:triangle_num);
            int2 = int_vals(triangle_num + 1:end);
        end

        function postProcessW2OptimalTransport(obj)
            % Post-processing after computing semi-discrete Wasserstein-2 optimal transport including computing the contour of each
            % cell in the Laguerre diagram and preparing for conditional rejection sampling 

            atom_num = size(obj.W2OT.DiscMeas.Atoms, 1);

            % shift the weights to make them positive so that they can be interpreted as squared radii of circles in a power diagram
            squared_radii = obj.W2OT.Weights - min(obj.W2OT.Weights) + 1;

            % invoke the mex function to compute a triangular mesh formed by the intersections of power cells and 
            [output_vertices, output_triangulation, power_cell_indices, ~] = mesh_intersect_power_diagram( ...
                obj.Supp.TriangleVertices, ...
                obj.Supp.Triangles, ...
                obj.W2OT.DiscMeas.Atoms, ...
                squared_radii, ...
                obj.W2OT.BoundingBox.BottomLeft, ...
                obj.W2OT.BoundingBox.TopRight);

            obj.W2OT.CondRejSamp = cell(atom_num, 1);

            for atom_id = 1:atom_num
                power_cell_vertex_indices = output_triangulation(power_cell_indices == atom_id, :);
                power_cell_vertices = output_vertices(power_cell_vertex_indices, :);

                % store the bounding box of the power cell as well as the polyshape representation of the power cell
                info = struct;
                info.BoxBottomLeft = min(power_cell_vertices, [], 1)';
                info.BoxWidthHeight = max(power_cell_vertices, [], 1)' - info.BoxBottomLeft;
                ch = convhull(power_cell_vertices, 'Simplify', true);
                info.PolyShape = polyshape(power_cell_vertices(ch, :), 'Simplify', false, 'KeepCollinearPoints', false);

                % maximum density in the power cell
                max_density = max(power_cell_vertices * obj.Dens.Weight + obj.Dens.Intercept) * 1.01 ...
                    / obj.W2OT.DiscMeas.Probs(atom_id);

                info.MaxDensity = max_density;
                info.Multiplier = max_density * info.BoxWidthHeight(1) * info.BoxWidthHeight(2);

                obj.W2OT.CondRejSamp{atom_id} = info;
            end
        end

        function samps = condRejSampFromW2Coupling(obj, ...
                atom_id, samp_num, rand_stream)
            % Perform rejection sampling of the conditional distribution given an atom in the discrete measure; i.e., perform rejection
            % sampling within a power cell.
            % Inputs:
            %   atom_id: atom index
            %   samp_num: number of samples required
            %   rand_stream: RandStream object used for sampling
            % Output:
            %   samps: two-column matrix containing the generated samples

            condrejsamp = obj.W2OT.CondRejSamp{atom_id};

            % the number of samples generated in each iteration is equal to 150% of samp_num times the expected number of iterations 
            % for each sample, capped at 1e4
            batch_samp_num = min(ceil(samp_num * condrejsamp.Multiplier * 1.5), 1e4);

            samp_counter = 0;
            samps = zeros(samp_num, 2);

            while samp_counter < samp_num
                % raw samples before the rejection step
                raw_samps = condrejsamp.BoxWidthHeight' .* rand(rand_stream, batch_samp_num, 2) + condrejsamp.BoxBottomLeft';

                % compute the acceptance probabilities
                raw_density = (raw_samps * obj.Dens.Weight + obj.Dens.Intercept) .* isinterior(condrejsamp.PolyShape, raw_samps);

                accept_probs = raw_density / (obj.W2OT.DiscMeas.Probs(atom_id) * condrejsamp.Multiplier ...
                    / (condrejsamp.BoxWidthHeight(1) * condrejsamp.BoxWidthHeight(2)));

                if any(accept_probs > 1)
                    warning('acceptance probability larger than 1 (%.4f) has been detected', max(accept_probs));
                end

                % the list of samples accepted
                accepted = rand(rand_stream, batch_samp_num, 1) <= accept_probs;

                acc_samps = raw_samps(accepted, :);
                acc_num = size(acc_samps, 1);

                % the number of new samples to be added to the list
                new_samp_num = min(acc_num, samp_num - samp_counter);

                samps(samp_counter + (1:new_samp_num), :) = acc_samps(1:new_samp_num, :);

                % update the sample counter
                samp_counter = samp_counter + new_samp_num;
            end
        end

        function integrals = doSetSimplicialTestFuncs(obj, ...
                vertices, triangles, check_simpcover)
            % Initialize the test functions with respect to a simplicial cover and compute the integrals of the test functions with
            % respect to the probability measure
            % Inputs: 
            %   vertices: two-column matrix containing the vertices used in the triangulation; the first rows must be identical to the
            %   vertices of the support
            %   triangles: three-column matrix containing the triangulation
            %   check_simpcover: check whether the given triangulation forms a simplicial cover (default is true)
            % Output:
            %   integrals: vector containing integrals of the test functions with respect to the probability measure; the number of 
            %   test functions is equal to the number of vertices in the triangulation

            if ~exist('check_simpcover', 'var') || isempty(check_simpcover)
                check_simpcover = true;
            end

            % check if the first vertices are identical to the support
            assert(max(max(abs(vertices(1:size(obj.Supp.Vertices, 1), :) - obj.Supp.Vertices))) < 1e-14, ...
                'the vertices do not begin with the vertices of the support');

            % check if all vertices are contained in the support
            assert(all(obj.checkIfInsideSupport(vertices)), 'there are vertices outside the support');

            % check if there are duplicate vertices
            assert(size(unique(vertices, 'rows'), 1) == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not appear in any triangle)
            assert(length(unique(triangles(:))) == size(vertices, 1), 'there are redundant vertices');

            if check_simpcover
                % check if the triangulation forms a simplicial cover
                [empty, is_simpcover] = check_triangulation_simplicial_cover(vertices, triangles);
    
                if any(empty)
                    error('some triangles are empty');
                end
    
                if ~is_simpcover
                    error('the triangulation does not form a simplicial cover');
                end
            end
            
            obj.SimplicialTestFuncs = struct;
            obj.SimplicialTestFuncs.Vertices = vertices;
            obj.SimplicialTestFuncs.Triangles = triangles;

            vert_num = size(vertices, 1);
            tri_num = size(triangles, 1);

            dens_vertices = vertices * obj.Dens.Weight + obj.Dens.Intercept;

            % evaluate the integral of three interpolation functions with respect to the probability measure within each triangular
            % region via an affine transformation which places the vertices on (0, 0), (1, 0), and (0, 1)
            tri_integrals = zeros(tri_num, 3);

            % compute the inverse transformation matrix which transform a coordinate into weights in the affine combination
            invtrans_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert = vertices(triangles(tri_id, :), :);
                tri_dens = dens_vertices(triangles(tri_id, :));
                trans_mat = tri_vert(2:3, :)' - tri_vert(1, :)';

                tri_integrals(tri_id, :) = (sum(tri_dens) + tri_dens) / 24 * abs(det(trans_mat));

                % compute the inverse transformation matrix
                invtrans_cell{tri_id} = [tri_vert'; ones(1, 3)] \ eye(3);
            end

            integrals = accumarray(triangles(:), tri_integrals(:), [vert_num, 1]);
            obj.SimplicialTestFuncs.Integrals = integrals;

            if abs(sum(integrals) - 1) > 1e-10
                error('integrals are incorrectly computed');
            end

            obj.SimplicialTestFuncs.InvTransMat = vertcat(invtrans_cell{:});
        end
    end

    methods(Static, Access = public)

        function vals = affineIntegralTriangle(dens_vertices, ...
                vertices_x, vertices_y, polynomial_coefs, after_transform)
            % Integrate bivariate second-order polynomials inside a triangle over a measure with positive affine density
            % Inputs:
            %   dens_vertices: three-column matrix containing the values of the affine density function evaluated at the three vertices
            %   of the triangle
            %   vertices_x: three-column matrix containing the x-coordinates of the three vertices
            %   vertices_y: three-column matrix containing the y-coordinates of the three vertices
            %   polynomial_coefs: six-column matrix containing the coefficients in the integrand polynomial: f(x, y) =  a1 * x^2 
            %   + a2 * y^2 + a3 * x * y + a4 * x + a5 * y + a6
            %   after_transform: logical vector indicating whether the coefficients of each integrand polynomial is applied to the 
            %   variables before transformation or after transforming to the triangle with vertices [0; 0], [1; 0], and [1; 1]
            % Output: 
            %   vals: vector containing the values of the integrals

            % express the affine transformation: [u;v] -> S([u;v]) = [S11, S12; S21, S22] * [u;v] + [S13; S23] that satisfies 
            % S([0;0]) = [x1;y1], S([1;0]) = [x2;y2], and S([1;1]) = [x3;y3]
            S11 = vertices_x(:, 2) - vertices_x(:, 1);
            S12 = vertices_x(:, 3) - vertices_x(:, 2);
            S13 = vertices_x(:, 1);
            S21 = vertices_y(:, 2) - vertices_y(:, 1);
            S22 = vertices_y(:, 3) - vertices_y(:, 2);
            S23 = vertices_y(:, 1);

            % absolute value of the determinant of the Jacobian
            abs_det_term = abs(S11 .* S22 - S12 .* S21);

            % compute c1, c2, c3, c4, c5, c6 such that f(S(x, y)) = c1 * x^2 + c2 * y^2 + c3 * x * y + c4 * x + c5 * y + c6
            a1 = polynomial_coefs(:, 1);
            a2 = polynomial_coefs(:, 2);
            a3 = polynomial_coefs(:, 3);
            a4 = polynomial_coefs(:, 4);
            a5 = polynomial_coefs(:, 5);
            a6 = polynomial_coefs(:, 6);
            
            c1 = a1 .* S11.^2 + a3 .* S11 .* S21 + a2 .* S21.^2;
            c2 = a1 .* S12.^2 + a3 .* S12 .* S22 + a2 .* S22.^2;
            c3 = 2 * a1 .* S11 .* S12 + a3 .* (S12 .* S21 + S11 .* S22) + 2 * a2 .* S21 .* S22;
            c4 = a4 .* S11 + a5 .* S21 + S13 .* (2 * a1 .* S11 + a3 .* S21) + S23 .* (2 * a2 .* S21 + a3 .* S11);
            c5 = a4 .* S12 + a5 .* S22 + S13 .* (2 * a1 .* S12 + a3 .* S22) + S23 .* (2 * a2 .* S22 + a3 .* S12);
            c6 = a1 .* S13.^2 + a3 .* S13 .* S23 + a2 .* S23.^2 + a4 .* S13 + a5 .* S23 + a6;

            % if after_transform is false, the coefficients do not need to be updated
            c1(after_transform) = a1(after_transform);
            c2(after_transform) = a2(after_transform);
            c3(after_transform) = a3(after_transform);
            c4(after_transform) = a4(after_transform);
            c5(after_transform) = a5(after_transform);
            c6(after_transform) = a6(after_transform);

            % evaluate the integral
            p1 = dens_vertices(:, 1);
            pdiff21 = dens_vertices(:, 2) - dens_vertices(:, 1);
            pdiff32 = dens_vertices(:, 3) - dens_vertices(:, 2);
            vals = ((c1 .* pdiff21 + (c1 .* pdiff32 + c3 .* pdiff21) / 2 + (c3 .* pdiff32 + c2 .* pdiff21) / 3 ...
                + (c2 .* pdiff32) / 4) / 5 + ((c1 .* p1 + c4 .* pdiff21) + (c3 .* p1 + c4 .* pdiff32 + c5 .* pdiff21) / 2 ...
                + (c2 .* p1 + c5 .* pdiff32) / 3) / 4 + ((c4 .* p1 + c6 .* pdiff21) + (c5 .* p1 + c6 .* pdiff32) / 2) / 3 ...
                + (c6 .* p1) / 2) .* abs_det_term;
        end
    end
end