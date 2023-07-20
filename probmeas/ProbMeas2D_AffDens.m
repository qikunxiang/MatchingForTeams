classdef ProbMeas2D_AffDens < ProbMeas2D_ConvexPolytope
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

            % evaluate the integral of the unnormalized density within each
            % triangular region via an affine transformation which places
            % the vertices on (0, 0), (1, 0), and (0, 1)
            tri_integral = zeros(size(DT, 1), 1);
            for tri_id = 1:size(DT, 1)
                tri_vert = vertices(DT(tri_id, :), :);
                trans_shift = tri_vert(1, :)';
                trans_mat = tri_vert(2:3, :)' - trans_shift;
                trans_weight = dens_weight' * trans_mat;
                trans_intercept = dens_weight' * trans_shift ...
                    + dens_intercept;

                tri_integral(tri_id) = (sum(trans_weight) / 6 ...
                    + trans_intercept / 2) * abs(det(trans_mat));
            end

            % calculate the normalizing constant and normalize the density
            norm_const = sum(tri_integral);
            obj.Dens = struct('Weight', dens_weight / norm_const, ...
                'Intercept', dens_intercept / norm_const, ...
                'NormConst', norm_const);

            % density function with the support
            obj.Dens.Func = @(z)(all(z * obj.Supp.Hyperplane.w' ...
                <= obj.Supp.Hyperplane.b', 2) ...
                .* (z * obj.Dens.Weight + obj.Dens.Intercept));

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
            obj.RejSamp.AcceptProbFunc = @(raw_samps)( ...
                obj.Dens.Func(raw_samps) ...
                    / (obj.RejSamp.Multiplier ...
                    / (diff(obj.RejSamp.Range.x) ...
                    * diff(obj.RejSamp.Range.y))));
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
                - obj.Supp.Hyperplane.b' <= 1e-14, 2);
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

            obj.OT.Prep.radial = struct;

            % intercept in the density function evaluated along a ray
            obj.OT.Prep.radial.intercept = atoms * obj.Dens.Weight ...
                + obj.Dens.Intercept;

            % coefficient in the density function evaluated along a ray
            obj.OT.Prep.radial.coefficient = ...
                [cos(angle_list), sin(angle_list)] * obj.Dens.Weight;
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

            upper_limits_p2 = upper_limits .^ 2;
            upper_limits_p3 = upper_limits .^ 3;
            upper_limits_p4 = upper_limits .^ 4;

            val1 = obj.OT.Prep.radial.coefficient / 3 ...
                .* upper_limits_p3 ...
                + obj.OT.Prep.radial.intercept(atom_id) / 2 ...
                * upper_limits_p2;
            val2 = obj.OT.Prep.radial.coefficient / 4 ...
                .* upper_limits_p4 ...
                + obj.OT.Prep.radial.intercept(atom_id) / 3 ...
                * upper_limits_p3;
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
                
                % prepare information about 
                obj.OT.CondRejSamp{atom_id} = struct;

                % add 10% of the standard deviation to make sure that the
                % circle covers the whole region
                rho_i = rho(:, atom_id);
                std_rho = std(rho_i);
                radius = max(rho_i) + std_rho * 0.1;

                obj.OT.CondRejSamp{atom_id}.Radius = radius;
                obj.OT.CondRejSamp{atom_id}.NeighborAtoms ...
                    = [find(any(hyp_full < radius, 1))'; atom_id];
                
                % conservative estimate of the maximum density in the cell
                % in each direction
                max_density = max(max(obj.OT.Prep.radial.coefficient, ...
                    0) .* rho_i + obj.OT.Prep.radial.intercept( ...
                    atom_id)) * 1.01 / obj.OT.DiscMeas.Probs(atom_id);

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
end

