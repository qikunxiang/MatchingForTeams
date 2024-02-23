classdef (Abstract) ProbMeas2D_ConvexPolytope < handle
    % Abstract class for probability measures supported within a convex
    % polytope

    properties(Constant)
        % numerical tolerance for deciding whether a point is inside a
        % polytope
        INSIDE_TOLERANCE = 1e-12;
    end

    properties(SetAccess = protected, GetAccess = public)
        % struct storing information about the support of the probability
        % measure
        Supp;

        % struct storing information about the density function of the
        % probability measure
        Dens;

        % struct storing information about rejection sampling from this
        % probability measure
        RejSamp;

        % struct storing information about the optimal transference plan
        % and the way to sample from the conditional measure
        OT = struct;

        % struct storing information about test functions for this
        % probability measure with respect to a simplicial cover
        SimplicialTestFuncs = [];
    end

    methods(Access = public)
        function obj = ProbMeas2D_ConvexPolytope(vertices)
            % Constructor function which sets the Supp property
            % Input:
            %   vertices: two-column matrix containing the vertices of the
            %   convex polytope which is the support of the probability
            %   measure

            % compute the convex hull of the vertices
            ccorder = convhull(vertices(:, 1), vertices(:, 2), ...
                'Simplify', true);

            % rearrange the vertices in counterclockwise order
            vertices_cc = vertices(ccorder, :);
            vertices_cc1 = vertices_cc(1:end - 1, :);
            vertices_cc2 = vertices_cc(2:end, :);

            obj.Supp = struct;
            
            % only store the vertices in the convex hull, in
            % counterclockwise order
            obj.Supp.Vertices = vertices_cc1;

            % store the maximum norm of points in the support of the
            % measure
            obj.Supp.MaxNorm = sqrt(max(sum(obj.Supp.Vertices .^ 2, 2)));

            obj.Supp.Hyperplane = struct;
            obj.Supp.Hyperplane.w = ...
                [vertices_cc2(:, 2) - vertices_cc1(:, 2), ...
                 vertices_cc1(:, 1) - vertices_cc2(:, 1)];
            obj.Supp.Hyperplane.b = ...
                vertices_cc2(:, 2) .* vertices_cc1(:, 1) ...
                - vertices_cc1(:, 2) .* vertices_cc2(:, 1);

            % correct numerical inaccuracies and make sure that the
            % vertices themselves are contained in all of the half-spaces
            obj.Supp.Hyperplane.b = max(obj.Supp.Hyperplane.b, ...
                max(obj.Supp.Hyperplane.w * vertices_cc1', [], 2));
        end

        function inside = checkIfInsideSupport(obj, pts, batch_size)
            % Check if the input points are inside the support of the
            % probability measure. Evaluation is done in batches if
            % necessary.
            % Inputs:
            %   pts: two-column matrix containing the input points
            %   batch_size: maximum number of inputs to be evaluated
            %   together via a vectorized routine (default is 1e4)
            % Output: 
            %   inside: boolean vector indicating whether each point is
            %   inside the support

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            inside_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                [inside_cell{batch_id}] ...
                    = obj.doCheckIfInsideSupport( ...
                    pts(((batch_id - 1) * batch_size + 1):min(batch_id ...
                    * batch_size, input_num), :));
            end

            inside = vertcat(inside_cell{:});
        end

        function samps = randSample(obj, samp_num, rand_stream)
            % Rejection sampling the probability measure.
            % Inputs:
            %   samp_num: number of samples to genenerate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream) 
            % Output:
            %   samps: two-column matrix containing the generated samples

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            % the number of samples generated in each iteration is equal to
            % 150% of samp_num times the expected number of iterations for
            % each sample, capped at 1e5
            batch_samp_num = min(ceil(samp_num ...
                * obj.RejSamp.Multiplier * 1.5), 1e5);

            samp_counter = 0;
            samps = zeros(samp_num, 2);

            while samp_counter < samp_num
                % raw samples before rejection
                raw_samps = obj.sampleFromProposalDistribution( ...
                    batch_samp_num, rand_stream);

                % compute the acceptance probabilities
                accept_probs = obj.rejSampAcceptProbFunction(raw_samps);

                % the list of samples accepted
                accepted = rand(rand_stream, batch_samp_num, 1) ...
                    <= accept_probs;

                acc_samps = raw_samps(accepted, :);
                acc_num = size(acc_samps, 1);

                % the number of new samples to be added to the list
                new_samp_num = min(acc_num, samp_num - samp_counter);

                samps(samp_counter + (1:new_samp_num), :) ...
                    = acc_samps(1:new_samp_num, :);

                % update the sample counter
                samp_counter = samp_counter + new_samp_num;
            end
        end

        function [OT_weights, OT_cost, OT_output] ...
                = computeOptimalTransport(obj, atoms, probs, angle_num, ...
                init_weights, pp_angle_indices, optim_options)
            % Compute the semi-discrete optimal transport
            % Inputs:
            %   atoms: two-column matrix containing the locations of the
            %   atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in
            %   the discrete measure
            %   angle_num: number of angles used in the integration
            %   (default is 1e4)
            %   init_weights: the initial weights in the optimization
            %   (default is the all-zero vector)
            %   pp_angle_indices: indices of the angles used in
            %   post-processing (default is all indices)
            %   optim_options: struct containing options used in the
            %   optimization
            % Outputs:
            %   OT_weights: the optimized weights in the Laguerre diagram
            %   OT_cost: the optimal transport cost
            %   OT_output: struct containing additional information with
            %   fields
            %       exitflag: the exit flag of the fminunc function
            %       optim_output: struct output by the fminunc function
            %       final_grad: the gradient at termination

            if ~exist('angle_num', 'var') || isempty(angle_num)
                angle_num = 1e4;
            end

            if ~exist('pp_angle_indices', 'var') ...
                    || isempty(pp_angle_indices)
                pp_angle_indices = (1:angle_num)';
            end

            atom_num = size(atoms, 1);

            if ~exist('init_weights', 'var') || isempty(init_weights)
                init_weights = zeros(atom_num - 1, 1);
            end

            % the multiplier increases the angle number when numerical
            % issues arise
            angle_num_multiplier = 1;

            while true
                angle_num_rt = angle_num * angle_num_multiplier;
                angle_list_rt = linspace(0, 2 * pi, angle_num_rt)';
                pp_angle_indices_rt = (pp_angle_indices - 1) ...
                    * angle_num_multiplier + 1;

                % pre-process the quantities
                obj.prepareOptimalTransport(atoms, probs, angle_list_rt);
                
                options = optimoptions('fminunc', 'Display', 'none', ...
                    'Algorithm', 'quasi-newton', ...
                    'SpecifyObjectiveGradient', true, ...
                    'MaxFunctionEvaluations', 1e5, ...
                    'MaxIterations', 1e5, ...
                    'FunctionTolerance', 0, ...
                    'StepTolerance', 0, ...
                    'OptimalityTolerance', 1e-5);
    
                % set the additional optimization options
                if exist('optim_options', 'var') && ~isempty(optim_options)
                    optim_fields = fieldnames(optim_options);
                    optim_values = struct2cell(optim_options);
    
                    for f_id = 1:length(optim_fields)
                        options.(optim_fields{f_id}) = optim_values{f_id};
                    end
                end
    
                OT_obj_func = @(w)(obj.evaluateOTMinObjective(w));
    
                OT_output = struct;
                [OT_weights, objective_min, OT_output.exitflag, ...
                    OT_output.optim_output, OT_output.final_grad] ...
                    = fminunc(OT_obj_func, init_weights, options);
                OT_cost = -objective_min;
    
                if OT_output.exitflag == 1
                    obj.OT.Coupled = true;
                    obj.OT.Weights = [0; OT_weights];
                    obj.OT.Cost = OT_cost;

                    % call post-processing function
                    obj.postProcessOptimalTransport(pp_angle_indices_rt);

                    break;
                end

                warning(['optimization solver stopped when the ' ...
                    'first-order condition fails']);

                % double the number of angles for numerical integration and
                % try again
                angle_num_multiplier = angle_num_multiplier * 2;
            end
        end

        function info = saveOptimalTransportInfo(obj)
            % Save optimal transport related information in a struct that
            % can be used later
            % Output:
            %   info: struct that copies the fields from obj.OT

            info = struct;
            info.atoms = obj.OT.DiscMeas.Atoms;
            info.probs = obj.OT.DiscMeas.Probs;
            info.angle_list = obj.OT.Prep.angle_list;
            info.weights = obj.OT.Weights;
            info.cost = obj.OT.Cost;
            info.condrejsamp = obj.OT.CondRejSamp;
        end

        function loadOptimalTransportInfo(obj, info)
            % Load optimal transport related information from a previously
            % saved struct
            % Input:
            %   info: struct that will be copied to the fields of obj.OT

            obj.prepareOptimalTransport(info.atoms, info.probs, ...
                info.angle_list);

            obj.OT.Weights = info.weights;
            obj.OT.Cost = info.cost;
            obj.OT.CondRejSamp = info.condrejsamp;
            obj.OT.Coupled = true;
        end

        function testOTObjective(obj, atoms, probs, angle_num, ...
                OT_weights, grid_pts_x, grid_pts_y)
            % Test the objective and gradient in the computation of
            % semi-discrete optimal transport by comparing with
            % two-dimensional integration over a grid
            % Inputs: 
            %   atoms: two-column matrix containing the locations of the
            %   atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in
            %   the discrete measure
            %   angle_num: number of angles used in the integration
            %   (default is 1e4)
            %   OT_weights: weight vector in the objective (including the
            %   first component)
            %   grid_pts_x: grid points on the x-axis for the integration
            %   grid_pts_y: grid points on the y-axis for the integration

            if ~exist('angle_num', 'var') || isempty(angle_num)
                angle_num = 1e4;
            end

            atom_num = size(atoms, 1);

            angle_list = linspace(0, 2 * pi, angle_num)';

            % pre-process the quantities
            obj.prepareOptimalTransport(atoms, probs, angle_list);

            % the objective value and gradient computed by the integration
            % in the polar coordinates
            [val, grad] = obj.evaluateOTMinObjective(OT_weights(2:end));

            [grid_x, grid_y] = meshgrid(grid_pts_x, grid_pts_y);

            % approximate the probability measure with a discrete measure
            % supported on the grid
            grid_density = obj.densityFunction([grid_x(:), grid_y(:)]);
            grid_density = grid_density / sum(grid_density);

            weighted_dist = sqrt((grid_x(:) - atoms(:, 1)') .^ 2 ...
                + (grid_y(:) - atoms(:, 2)') .^ 2) - OT_weights';
            [min_weighted_dist, min_id] = min(weighted_dist, [], 2);

            val_grid = -(OT_weights' * probs) - min_weighted_dist' ...
                * grid_density;
            grad_grid = accumarray(min_id, grid_density, [atom_num, 1]) ...
                - probs;
            grad_grid = grad_grid(2:end);

            fprintf('objective via integration: %.6f\n', val);
            fprintf('objective via grid: %.6f\n', val_grid);
            fprintf('difference: %e\n', abs(val - val_grid));

            fprintf('gradients:\n');
            display([grad, grad_grid]);

            fprintf('max diff. in gradients: %e\n', ...
                max(abs(grad - grad_grid)));
        end

        function testOTGradient(obj, atoms, probs, angle_num, ...
                OT_weights, nudge)
            % Check the gradient of the objective function in the
            % computation of semi-discrete optimal transport
            % Inputs: 
            %   atoms: two-column matrix containing the locations of the
            %   atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in
            %   the discrete measure
            %   angle_num: number of angles used in the integration
            %   (default is 1e4)
            %   OT_weights: weight vector in the objective (including the
            %   first component)
            %   nudge: the size of the disturbance in finite differencing
            %   (default is 1e-6)


            if ~exist('angle_num', 'var') || isempty(angle_num)
                angle_num = 1e4;
            end

            if ~exist('nudge', 'var') || isempty(nudge)
                nudge = 1e-6;
            end

            atom_num = size(atoms, 1);

            angle_list = linspace(0, 2 * pi, angle_num);

            % pre-process the quantities
            obj.prepareOptimalTransport(atoms, probs, angle_list);

            % the objective value and gradient computed by the integration
            % in the polar coordinates
            [val, grad] = obj.evaluateOTMinObjective(OT_weights(2:end));

            grad_fd = zeros(atom_num - 1, 1);

            for comp_id = 2:atom_num
                OT_weights_nudged = OT_weights;
                OT_weights_nudged(comp_id) = OT_weights_nudged(comp_id) ...
                    + nudge;

                val_nudged = obj.evaluateOTMinObjective( ...
                    OT_weights_nudged(2:end));

                grad_fd(comp_id - 1) = (val_nudged - val) / nudge;
            end

            fprintf('gradients:\n');
            display([grad, grad_fd]);

            fprintf('max diff. in gradients: %e\n', ...
                max(abs(grad - grad_fd)));
        end

        function [x, y, rho] = computeLaguerreDiagram(obj, ...
                anchors, weights, angle_num)
            % Compute the contour of each cell in the Laguerre diagram with
            % given weights
            % Inputs: 
            %   anchors: two-column matrix containing the locations of the
            %   anchor points
            %   weights: weights in the Laguerre diagram
            %   angle_num: number of angles used in the integration
            %   (default is 1e4)
            % Outputs:
            %   x: matrix containing the x-coordinates of the points 
            %   describing the Laguerre cells; each column represents a 
            %   cell
            %   y: matrix containing the y-coordinates of the points 
            %   describing the Laguerre cells; each column represents a 
            %   cell
            %   rho: matrix containing the distances from the anchor point
            %   to the points describing Laguerre cells; each column
            %   represents a cell

            if ~exist('angle_num', 'var') || isempty(angle_num)
                angle_num = 1e4;
            end

            anchor_num = size(anchors, 1);

            % pre-process the quantities
            angle_list = linspace(0, 2 * pi, angle_num)';

            assert(all(obj.checkIfInsideSupport(anchors)), ...
                'some atoms are outside the support of the measure');


            % prepare information related to the hyperbolas
            hyp = ProbMeas2D_ConvexPolytope.hyperbolaPreprocess(anchors);

            % prepare information related to the linear boundaries
            lin = ProbMeas2D_ConvexPolytope.linePreprocess( ...
                anchors, obj.Supp.Hyperplane.w, obj.Supp.Hyperplane.b);

            rho = zeros(angle_num, anchor_num);

            hyp_a = (weights - weights') / 2;

            for anchor_id = 1:anchor_num
                rho(:, anchor_id) ...
                    = hyperbolic_Dirichlet_tessellation_polar( ...
                    angle_list, hyp_a(:, anchor_id), ...
                    hyp.c(:, anchor_id), hyp.intercept(:, anchor_id), ...
                    lin{anchor_id}.c, lin{anchor_id}.intercept);
            end

            x = anchors(:, 1)' + rho .* cos(angle_list);
            y = anchors(:, 2)' + rho .* sin(angle_list);
        end

        function samp_cell = conditionalRandSample(obj, ...
                samp_num_list, rand_stream)
            % Randomly generate samples from the conditional distributions
            % given the coupled atoms in the discrete measure
            % Inputs: 
            %   samp_num_list: vector containing number of samples coupled
            %   with each of the atoms in the discrete measure
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream) 
            % Output:
            %   samp_cell: cell array where each cell contains a two-column
            %   array containing the generated samples

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~isfield(obj.OT, 'Coupled') || ~obj.OT.Coupled
                error('must set the coupled discrete measure first');
            end

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);
            
            assert(length(samp_num_list) == atom_num, ...
                'number of atoms mismatch');

            samp_cell = cell(atom_num, 1);

            for atom_id = 1:atom_num
                % the number of samples generated in each iteration is
                % equal to 150% of samp_num_list(atom_id) times the
                % expected number of iterations for each sample, capped at
                % 1e4
                batch_samp_num = min(ceil(samp_num_list(atom_id) ...
                    * obj.OT.CondRejSamp{atom_id}.Multiplier * 1.5), 1e4);

                samp_counter = 0;
                samp_cell{atom_id} = zeros(samp_num_list(atom_id), 2);

                while samp_counter < samp_num_list(atom_id)
                    % raw samples before the rejection step
                    raw_samps = obj.sampleFromCondProposalDistribution( ...
                        atom_id, batch_samp_num, rand_stream);

                    % compute the acceptance probabilities
                    accept_probs = ...
                        obj.condRejSampAcceptProbFunction(atom_id, ...
                        raw_samps);

                    if any(accept_probs > 1)
                        warning(['acceptance probability ' ...
                            'larger than 1 (%.4f) has been detected'], ...
                            max(accept_probs));
                    end

                    % the list of samples accepted
                    accepted = rand(rand_stream, batch_samp_num, 1) ...
                        <= accept_probs;

                    acc_samps = raw_samps(accepted, :);
                    acc_num = size(acc_samps, 1);

                    % the number of new samples to be added to the list
                    new_samp_num = min(acc_num, samp_num_list(atom_id) ...
                        - samp_counter);

                    samp_cell{atom_id}(samp_counter ...
                        + (1:new_samp_num), :) ...
                        = acc_samps(1:new_samp_num, :);

                    % update the sample counter
                    samp_counter = samp_counter + new_samp_num;
                end
            end
        end

        function [samps, disc_samps] = randSampleCoupled(obj, samp_num, ...
                rand_stream)
            % Randomly generate samples from the coupling of the
            % probability measure and the discrete measure
            % Inputs: 
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream) 
            % Output:
            %   samps: two-column matrix containing the generated samples
            %   disc_samps: vector containing the corresponding sampled
            %   atom indices in the discrete measure

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~isfield(obj.OT, 'Coupled') || ~obj.OT.Coupled
                error('must set the coupled discrete measure first');
            end

            atom_num = length(obj.OT.DiscMeas.Probs);

            % generate samples from the discrete measure
            disc_samps = randsample(rand_stream, atom_num, samp_num, ...
                true, obj.OT.DiscMeas.Probs);

            % count the number of samples for each cell
            samp_num_list = accumarray(disc_samps, ones(samp_num, 1), ...
                [atom_num, 1]);

            samp_cell = obj.conditionalRandSample(samp_num_list, ...
                rand_stream);
            
            samps = vertcat(samp_cell{:});

            % order the samples according the discrete samples
            [~, ordering] = sort(disc_samps, 'ascend');
            samps(ordering, :) = samps;
        end

        function [vals, inside] = evaluateSimplicialTestFuncs(obj, pts, ...
                batch_size)
            % Evaluate the test functions at given locations; the input
            % locations must be inside the support of the measure
            % Input:
            %   pts: two-column matrix containing the input points
            %   batch_size: maximum number of inputs to be evaluated
            %   together via a vectorized routine (default is 1e4)
            % Output:
            %   vals: sparse matrix containing the computed function values
            %   where each row corresponds to an input and each column
            %   corresponds to a test function
            %   inside: sparse boolean matrix indicating whether an input
            %   point is inside each of the triangles; each row corresponds
            %   to an input and each column corresponds to a triangle

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            assert(all(obj.checkIfInsideSupport(pts)), ...
                'some points are not inside the support');

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);
            inside_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                [vals_cell{batch_id}, inside_cell{batch_id}] ...
                    = obj.doEvaluateSimplicialTestFuncs( ...
                    pts(((batch_id - 1) * batch_size + 1):min(batch_id ...
                    * batch_size, input_num), :));
            end

            vals = vertcat(vals_cell{:});
            inside = vertcat(inside_cell{:});
        end

        function vals = evaluateWeightedSumOfSimplicialTestFuncs(obj, ...
                pts, coefficients, batch_size)
            % Evaluate a weighted sum of simplicial test functions at given
            % locations; the input locations must be inside the support of
            % the measure 
            % Input:
            %   pts: two-column matrix containing the input points
            %   coefficients: vector containing the coefficients of the
            %   test functions in the weighted sum
            %   batch_size: maximum number of inputs to be evaluated
            %   together via a vectorized routine (default is 1e4)
            % Output:
            %   vals: sparse matrix containing the computed function values
            %   where each row corresponds to an input and each column
            %   corresponds to a test function

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            assert(all(obj.checkIfInsideSupport(pts)), ...
                'some points are not inside the support');

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                [raw_vals, ~] ...
                    = obj.doEvaluateSimplicialTestFuncs( ...
                    pts(((batch_id - 1) * batch_size + 1):min(batch_id ...
                    * batch_size, input_num), :));
                vals_cell{batch_id} = raw_vals * coefficients;
            end

            vals = vertcat(vals_cell{:});
        end

        function testSimplicialTestFuncs(obj, grid_pts_x, grid_pts_y)
            % Test the computation of the integrals of test functions with
            % respect to the probability measure as well as the evaluation
            % of test functions

            [grid_x, grid_y] = meshgrid(grid_pts_x, grid_pts_y);
            grid_pts = [grid_x(:), grid_y(:)];

            % approximate the probability measure with a discrete measure
            % supported on the grid
            grid_density = obj.densityFunction(grid_pts);
            grid_density = grid_density / sum(grid_density);

            % only retain points that are inside the support
            inside_support = grid_density > 0;
            grid_pts = grid_pts(inside_support, :);
            grid_density = grid_density(inside_support);

            testfunc_vals = obj.evaluateSimplicialTestFuncs(grid_pts);

            integrals_grid = sum(testfunc_vals .* grid_density, 1)';

            fprintf('integrals:\n');
            display(full([obj.SimplicialTestFuncs.Integrals, ...
                integrals_grid]));

            fprintf('max diff. in integrals: %e\n', ...
                max(abs(obj.SimplicialTestFuncs.Integrals ...
                - integrals_grid)));
        end
    
        function integrals = setSimplicialTestFuncs(obj, varargin)
            % Initialize the test functions with respect to a simplicial 
            % cover and compute the integrals of the test functions with 
            % respect to the probability measure

            integrals = obj.doSetSimplicialTestFuncs(varargin{:});

            % compute the mesh size
            vertices = obj.SimplicialTestFuncs.Vertices;
            triangles = obj.SimplicialTestFuncs.Triangles;

            v1 = vertices(triangles(:, 1), :);
            v2 = vertices(triangles(:, 2), :);
            v3 = vertices(triangles(:, 3), :);

            edge1sq = sum((v1 - v2) .^ 2, 2);
            edge2sq = sum((v2 - v3) .^ 2, 2);
            edge3sq = sum((v3 - v1) .^ 2, 2);
            
            obj.SimplicialTestFuncs.MeshSize = ...
                sqrt(max(max(max(edge1sq, edge2sq), edge3sq)));
        end
    end

    methods(Static, Access = protected)
        function hyp = hyperbolaPreprocess(anchors)
            % Computes quantities for the polar representations of
            % hyperbolas where the foci are among the given anchor points.
            % Input:
            %   anchors: two-column matrix containing the locations of the
            %   anchor points
            % Output:
            %   hyp: struct with fields intercept and c which stores
            %   information about the polar representations of hyperbolas

            % column i contains x(i) - x
            diff_x = anchors(:, 1)' - anchors(:, 1);
            diff_y = anchors(:, 2)' - anchors(:, 2);

            hyp = struct;

            % rotation angle (in radian) in order to transform each
            % hyperbola into its canonical orientation 
            hyp.intercept = atan(diff_y ./ diff_x);
            flip_list = anchors(:, 1)' < anchors(:, 1);
            hyp.intercept(flip_list) = hyp.intercept(flip_list) + pi;

            % distance between two foci (same as the pairwise distance
            % matrix for the anchors points) divided by 2
            hyp.c = sqrt(diff_x .^ 2 + diff_y .^ 2) / 2;
        end

        function lin = linePreprocess(anchors, hp_w, hp_b)
            % Computes polar representations of the characterizing
            % hyperplanes of a convex polytope using each of the anchor
            % points as origin. The representation expresses each
            % characterizing hyperplane (which is a straight line) as a
            % line parallel to the x-axis (i.e., y = c) rotated by an angle
            % (in radian).
            % Inputs: 
            %   anchors: two-column matrix containing the locations of the
            %   anchor points
            %   hp_w: two-column matrix containing the weight vectors of
            %   the characterizing hyperplanes
            %   hp_b: vector containing the intercepts of the
            %   characterizing hyperplanes
            % Output:
            %   lin: cell array containing structs with fields intercept
            %   and c describing the polar representations

            % shifting the hyperplanes such that an anchor point is the new
            % origin; this only changes the intercepts of the hyperplanes
            hp_b_shifted = hp_b - hp_w * anchors';

            % some anchor points may lie on the wrong side of a hyperplane 
            % due to small numerical inaccuracies; set these atoms to be on
            % the hyperplane
            closeto0_list = hp_b_shifted < 0 & hp_b_shifted > -1e-13;
            hp_b_shifted(closeto0_list) = 0;

            % iterate through the anchor points
            anchor_num = size(anchors, 1);

            % cell array where each cell corresponds to an anchor point; 
            % each cell contains a struct with fields c and intercept which
            % are vectors representing the characterizing hyperplanes in 
            % polar coordinates using the anchor point as the origin 
            lin = cell(anchor_num, 1);

            for anchor_id = 1:anchor_num
                % adjust the signs of the weight vector and the intercept
                % such that the intercept is always positive
                hyp_b_sign_list = (hp_b_shifted(:, anchor_id) >= 0) ...
                    * 2 - 1;
                hyp_b_adjusted = hp_b_shifted(:, anchor_id) ...
                    .* hyp_b_sign_list;
                hyp_w_adjusted = hp_w .* hyp_b_sign_list;

                lin{anchor_id} = struct;
                lin{anchor_id}.c = hyp_b_adjusted ...
                    ./ sqrt(sum(hyp_w_adjusted .^ 2, 2));
                lin{anchor_id}.intercept = -atan(hyp_w_adjusted(:, 1) ...
                    ./ hyp_w_adjusted(:, 2)) ...
                    + pi * (hyp_w_adjusted(:, 2) < 0);
            end
        end
    end

    methods(Access = protected)
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

            assert(all(obj.checkIfInsideSupport(atoms)), ...
                'some atoms are outside the support of the measure')

            obj.OT.DiscMeas = struct('Atoms', atoms, 'Probs', probs);
            obj.OT.Prep = struct('angle_list', angle_list);

            % prepare information related to the hyperbolas
            obj.OT.Prep.hyp = ...
                ProbMeas2D_ConvexPolytope.hyperbolaPreprocess(atoms);

            % prepare information related to the linear boundaries
            obj.OT.Prep.lin = ProbMeas2D_ConvexPolytope.linePreprocess( ...
                atoms, obj.Supp.Hyperplane.w, obj.Supp.Hyperplane.b);
        end

        function [val, grad] = evaluateOTMinObjective(obj, OT_weights, ...
                batch_size)
            % The objective function of the minimization problem when
            % solving semi-discrete optimal transport (OT) problem.
            % Inputs:
            %   OT_weights: the input of the objective function
            %   corresponding to the dual weights in the Laguerre diagram
            %   (the first weight is set to 0 for identification purpose,
            %   and thus the input here corresponds to the second weight
            %   onwards)
            %   batch_size: the maximum number of angles to be handled in
            %   the vectorized procedure (default is 1e4)
            % Outputs:
            %   val: the computed objective value
            %   grad: the computed gradient of the objective function

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            % add back the first weight
            OT_weights = [0; OT_weights];
            OTweights_half = OT_weights / 2;

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);
            hyp_a = OTweights_half - OTweights_half';

            angle_list = obj.OT.Prep.angle_list;
            angle_num = length(angle_list);
            batch_num = ceil(angle_num / batch_size);

            int1 = zeros(atom_num, 1);
            int2 = zeros(atom_num, 1);

            for atom_id = 1:atom_num
                for batch_id = 1:batch_num
                    batch_indices = ((batch_id - 1) * batch_size ...
                        + 1):min(batch_id * batch_size, angle_num);
                    angle_batch = angle_list(batch_indices);

                    % use the polar coordinate representation to compute 
                    % the upper integration limits for the inner integral
                    upper_limits = ...
                        hyperbolic_Dirichlet_tessellation_polar( ...
                        angle_batch, hyp_a(:, atom_id), ...
                        obj.OT.Prep.hyp.c(:, atom_id), ...
                        obj.OT.Prep.hyp.intercept(:, atom_id), ...
                        obj.OT.Prep.lin{atom_id}.c, ...
                        obj.OT.Prep.lin{atom_id}.intercept);

                    % evaluate the outer integral via the Trapezoidal rule
                    [inner1, inner2] = obj.computeInnerIntegral( ...
                        atom_id, batch_indices, upper_limits);

                    int_vals = trapz(angle_batch, [inner1, inner2]);
                    int1(atom_id) = int1(atom_id) + int_vals(1);
                    int2(atom_id) = int2(atom_id) + int_vals(2);
                end
            end

            val = -(OT_weights' * obj.OT.DiscMeas.Probs) ...
                - sum(int2) + OT_weights' * int1;
            grad = int1 - obj.OT.DiscMeas.Probs;

            % the gradient wrt the first weight is discarded
            grad = grad(2:end);
        end

        function [vals, inside] = doEvaluateSimplicialTestFuncs(obj, pts)
            % Function that actually evaluates the test functions at the
            % given locations 
            % Input: 
            %   pts: two-column matrix containing the input points
            % Output:
            %   vals: sparse matrix containing the computed function values
            %   where each row corresponds to an input and each column
            %   corresponds to a test function
            %   inside: sparse boolean matrix indicating whether an input
            %   point is inside each of the triangles; each row corresponds
            %   to an input and each column corresponds to a triangle

            input_num = size(pts, 1);
            vert_num = size(obj.SimplicialTestFuncs.Vertices, 1);
            tri_num = size(obj.SimplicialTestFuncs.Triangles, 1);

            % apply the inverse transformation to recover the coefficients
            % in the affine combination (with respect to the three vertices
            % of each triangle)
            coef_mat = obj.SimplicialTestFuncs.InvTransMat ...
                * [pts'; ones(1, input_num)];

            % allow a small tolerance to avoid numerical precision issues
            inside = reshape(all(reshape(sparse(coef_mat(:) ...
                > -ProbMeas2D_ConvexPolytope.INSIDE_TOLERANCE), ...
                3, tri_num * input_num), 1), tri_num, input_num)';

            % normalize the matrix such that each row sums up to 1; this is
            % to deal with the case where a point belongs to more than one
            % triangles; in such a case, we compute the test function with
            % respect to each triangle that the point belongs to and then
            % take the average
            inside_norm = inside ./ sum(inside, 2);

            % compute the weights in the convex combination
            conv_coef_mat = sparse(coef_mat ...
                .* repelem(inside_norm', 3, 1));

            % compute the input indices and the vertex indices
            % corresponding to each non-zero entry in the convex
            % combination coefficient matrix
            input_id_mat = repmat((1:input_num), tri_num * 3, 1);
            triangles = obj.SimplicialTestFuncs.Triangles';
            vert_id_mat = repmat(triangles(:), 1, input_num);
            [conv_coef_r, conv_coef_c, conv_coef_vals] ...
                = find(conv_coef_mat);
            conv_coef_id = sub2ind([tri_num * 3, input_num], ...
                conv_coef_r, conv_coef_c);

            % sum up the coefficients and return a sparse matrix
            vals = accumarray([input_id_mat(conv_coef_id), ...
                vert_id_mat(conv_coef_id)], conv_coef_vals, ...
                [input_num, vert_num], [], [], true);
        end
    end

    methods(Abstract, Access = public)
        % Probability density function
        dens = densityFunction(obj, pts);

        % Fix small numerical inaccuracies in the points such that points
        % that are slightly outside the support are moved into the support
        pts_san = sanitizePoints(obj, pts, tolerance);
    end

    methods(Abstract, Access = protected)
        % Check if the input points are inside the support of the
        % probability measure
        inside = doCheckIfInsideSupport(obj, pts);

        % Generate random samples from a proposal distribution for
        % rejection sampling
        raw_samps = sampleFromProposalDistribution(obj, samp_num, ...
            rand_stream);

        % Computes the acceptance probability when using rejection sampling
        % to generate independent samples from the probability measure
        accept_probs = rejSampAcceptProbFunction(obj, raw_samps);

        % Evaluates inner integrals along a list of angles in the polar
        % coordinate system. The integrands are r -> r * p(r) and
        % r -> r^2 * p(r) where p(r) denotes the probability density
        % function.
        [val1, val2] = computeInnerIntegral(obj, atom_id, ...
            angle_indices, upper_limits)

        % Generate random samples from a proposal distribution for
        % conditional rejection sampling. For this procedure, the proposal
        % distribution is uniform on a circle.
        raw_samps = sampleFromCondProposalDistribution(obj, atom_id, ...
            samp_num, rand_stream);

        % Compute the acceptance probability in conditional rejection
        % sampling from a cell
        accept_probs = condRejSampAcceptProbFunction(obj, atom_id, ...
            raw_samps);

        % Function that initializes the test functions with respect to a 
        % simplicial cover and computes the integrals of the test functions 
        % with respect to the probability measure
        integrals = doSetSimplicialTestFuncs(obj, vargin);

        % Post-processing after computing semi-discrete optimal transport 
        % including computing the contour of each cell in the Laguerre 
        % diagram and preparing for conditional rejection sampling
        postProcessOptimalTransport(obj, pp_angle_indices);
    end
end

