classdef (Abstract) ProbMeas2D_ConvexPolytope < handle
    % Abstract class for probability measures supported within a convex
    % polytope

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

        % whether the measure has been coupled with a discrete measure
        Coupled = false;

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
                accept_probs = obj.RejSamp.AcceptProbFunc(raw_samps);

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
                optim_options)
            % Compute the semi-discrete optimal transport
            % Inputs:
            %   atoms: two-column matrix containing the locations of the
            %   atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in
            %   the discrete measure
            %   angle_num: number of angles used in the integration
            %   (default is 1e4)
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

            atom_num = size(atoms, 1);

            % pre-process the quantities
            obj.prepareOptimalTransport(atoms, probs, angle_num);
            
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

            % call fminunc to optimize the weights
            
            % maximum number of trials
            rand_stream = RandStream('combRecursive', 'Seed', 5555);
            max_trial = 10;
            trial_counter = 0;
            OT_weights = zeros(atom_num - 1, 1);
            perterbation = zeros(atom_num - 1, 1);

            while trial_counter < max_trial
                OT_output = struct;
                [OT_weights, objective_min, OT_output.exitflag, ...
                    OT_output.optim_output, OT_output.final_grad] ...
                    = fminunc(OT_obj_func, OT_weights + perterbation, ...
                    options);
                OT_cost = -objective_min;
                trial_counter = trial_counter + 1;
    
                if OT_output.exitflag ~= 1
                    warning(['optimization solver stopped when the ' ...
                        'first-order condition fails']);

                    % retry with the current weights plus a small random
                    % perturbation
                    perterbation = randn(rand_stream, atom_num - 1, 1) ...
                        * 1e-4;
                else
                    break;
                end
            end

            obj.OT.Weights = [0; OT_weights];
            obj.OT.Cost = OT_cost;
            obj.Coupled = true;

            % call post-processing function
            obj.postProcessOptimalTransport();
        end

        function setOTWeights(obj, atoms, probs, OT_weights, OT_cost, ...
                angle_num)
            % Set the pre-computed optimal transport weights and cost
            % Inputs:
            %   atoms: two-column matrix containing the locations of the
            %   atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in
            %   the discrete measure
            %   OT_weights: the optimized weights in the Laguerre diagram
            %   OT_cost: the optimal transport cost
            %   angle_num: number of angles used in the Laguerre diagram
            %   (default is 1e4)

            if ~exist('angle_num', 'var') || isempty(angle_num)
                angle_num = 1e4;
            end

            % pre-process the quantities
            obj.prepareOptimalTransport(atoms, probs, angle_num);

            obj.OT.Weights = OT_weights;
            obj.OT.Cost = OT_cost;
            obj.Coupled = true;

            % call post-processing function to compute additional
            % quantities
            obj.postProcessOptimalTransport();
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

            % pre-process the quantities
            obj.prepareOptimalTransport(atoms, probs, angle_num);

            % the objective value and gradient computed by the integration
            % in the polar coordinates
            [val, grad] = obj.evaluateOTMinObjective(OT_weights(2:end));

            [grid_x, grid_y] = meshgrid(grid_pts_x, grid_pts_y);

            % approximate the probability measure with a discrete measure
            % supported on the grid
            grid_density = obj.Dens.Func([grid_x(:), grid_y(:)]);
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

            % pre-process the quantities
            obj.prepareOptimalTransport(atoms, probs, angle_num);

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
                anchors, angle_num, weights)
            % Compute the contour of each cell in the Laguerre diagram with
            % given weights
            % Inputs: 
            %   atoms: two-column matrix containing the locations of the
            %   anchor points
            %   angle_num: number of angles used in the integration
            %   (default is 1e4)
            %   weights: weights in the Laguerre diagram

            if ~exist('angle_num', 'var') || isempty(angle_num)
                angle_num = 1e4;
            end

            anchor_num = size(anchors, 1);
            probs = ones(anchor_num, 1) / anchor_num;

            % pre-process the quantities
            obj.prepareOptimalTransport(anchors, probs, angle_num);

            [x, y, rho] = obj.doComputeLaguerreDiagram(weights);
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

            if ~obj.Coupled
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
                        obj.OT.CondRejSamp{atom_id}.AcceptProbFunc( ...
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

            if ~obj.Coupled
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
            grid_density = obj.Dens.Func(grid_pts);
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

    methods(Access = protected)
        function prepareOptimalTransport(obj, atoms, probs, angle_num)
            % Prepare quantities for the computation of semi-discrete
            % optimal transport.
            % Inputs:
            %   atoms: two-column matrix containing the locations of the
            %   atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in
            %   the discrete measure
            %   angle_num: number of angles used in the integration

            angle_list = linspace(0, 2 * pi, angle_num)';

            assert(all(obj.checkIfInsideSupport(atoms)), ...
                'some atoms are outside the support of the measure')

            obj.OT.DiscMeas = struct('Atoms', atoms, 'Probs', probs);
            obj.OT.Prep = struct('angle_list', angle_list);

            % prepare information related to the hyperbolas
            obj.hyperbolaPreprocess();

            % prepare information related to the linear boundaries
            obj.linePreprocess();
        end
        
        function hyperbolaPreprocess(obj)
            % Computes quantities for the polar representations of
            % hyperbolas where the foci are among the given anchor points.

            % column i contains x(i) - x
            diff_x = obj.OT.DiscMeas.Atoms(:, 1)' ...
                - obj.OT.DiscMeas.Atoms(:, 1);
            diff_y = obj.OT.DiscMeas.Atoms(:, 2)' ...
                - obj.OT.DiscMeas.Atoms(:, 2);

            obj.OT.Prep.hyp = struct;

            % rotation angle (in radian) in order to transform each
            % hyperbola into its canonical orientation 
            obj.OT.Prep.hyp.intercept = atan(diff_y ./ diff_x);
            flip_list = obj.OT.DiscMeas.Atoms(:, 1)' ...
                < obj.OT.DiscMeas.Atoms(:, 1);
            obj.OT.Prep.hyp.intercept(flip_list) ...
                = obj.OT.Prep.hyp.intercept(flip_list) + pi;

            % distance between two foci (same as the pairwise distance
            % matrix for the anchors points) divided by 2
            obj.OT.Prep.hyp.c = sqrt(diff_x .^ 2 + diff_y .^ 2) / 2;
        end

        function linePreprocess(obj)
            % Computes polar representations of the characterizing
            % hyperplanes of a convex polytope using each of the anchor
            % points as origin. The representation expresses each
            % characterizing hyperplane (which is a straight line) as a
            % line parallel to the x-axis (i.e., y = c) rotated by an angle
            % (in radian).

            % shifting the hyperplanes such that an anchor point is the new
            % origin; this only changes the intercepts of the hyperplanes
            hyp_b_shifted = obj.Supp.Hyperplane.b ...
                - obj.Supp.Hyperplane.w * obj.OT.DiscMeas.Atoms';

            % some atoms may lie on the wrong side of a hyperplane due to
            % small numerical inaccuracies; set these atoms to be on the
            % hyperplane
            closeto0_list = hyp_b_shifted < 0 & hyp_b_shifted > -1e-13;
            hyp_b_shifted(closeto0_list) = 0;

            % iterate through the anchor points
            atom_num = size(obj.OT.DiscMeas.Atoms, 1);

            % cell array where each cell corresponds to an atom in the
            % discrete measure; each cell contains a struct with fields c
            % and intercept which are vectors representing the
            % characterizing hyperplanes in polar coordinates using the
            % atom as the origin 
            obj.OT.Prep.lin = cell(atom_num, 1);

            for atom_id = 1:atom_num
                % adjust the signs of the weight vector and the intercept
                % such that the intercept is always positive
                hyp_b_sign_list = (hyp_b_shifted(:, atom_id) >= 0) * 2 - 1;
                hyp_b_adjusted = hyp_b_shifted(:, atom_id) ...
                    .* hyp_b_sign_list;
                hyp_w_adjusted = obj.Supp.Hyperplane.w .* hyp_b_sign_list;

                obj.OT.Prep.lin{atom_id} = struct;
                obj.OT.Prep.lin{atom_id}.c = hyp_b_adjusted ...
                    ./ sqrt(sum(hyp_w_adjusted .^ 2, 2));
                obj.OT.Prep.lin{atom_id}.intercept ...
                    = -atan(hyp_w_adjusted(:, 1) ...
                    ./ hyp_w_adjusted(:, 2)) ...
                    + pi * (hyp_w_adjusted(:, 2) < 0);
            end
        end

        function [val, grad] = evaluateOTMinObjective(obj, OT_weights)
            % The objective function of the minimization problem when
            % solving semi-discrete optimal transport (OT) problem.
            % Inputs:
            %   OT_weights: the input of the objective function
            %   corresponding to the dual weights in the Laguerre diagram
            %   (the first weight is set to 0 for identification purpose,
            %   and thus the input here corresponds to the second weight
            %   onwards)
            % Outputs:
            %   val: the computed objective value
            %   grad: the computed gradient of the objective function

            % add back the first weight
            OT_weights = [0; OT_weights];
            OTweights_half = OT_weights / 2;

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);
            hyp_a = OTweights_half - OTweights_half';

            int1 = zeros(atom_num, 1);
            int2 = zeros(atom_num, 1);

            for atom_id = 1:atom_num
                % use the polar coordinate representation to compute the
                % upper integration limits for the inner integral
                upper_limits = hyperbolic_Dirichlet_tessellation_polar( ...
                    obj.OT.Prep.angle_list, hyp_a(:, atom_id), ...
                    obj.OT.Prep.hyp.c(:, atom_id), ...
                    obj.OT.Prep.hyp.intercept(:, atom_id), ...
                    obj.OT.Prep.lin{atom_id}.c, ...
                    obj.OT.Prep.lin{atom_id}.intercept);

                % evaluate the outer integral via the Trapezoidal rule
                [inner1, inner2] = obj.computeInnerIntegral(atom_id, ...
                    upper_limits);

                int_vals = trapz(obj.OT.Prep.angle_list, [inner1, inner2]);
                int1(atom_id) = int_vals(1);
                int2(atom_id) = int_vals(2);
            end

            val = -(OT_weights' * obj.OT.DiscMeas.Probs) ...
                - sum(int2) + OT_weights' * int1;
            grad = int1 - obj.OT.DiscMeas.Probs;

            % the gradient wrt the first weight is discarded
            grad = grad(2:end);
        end

        function [x, y, rho] = doComputeLaguerreDiagram(obj, weights)
            % Function that actually computes the Laguerre diagram given
            % weights
            % Input:
            %   weights: weights in the Laguerre diagram
            %   full_info: boolean indicating whether to request full
            %   distance information (default is false)
            % Outputs:
            %   x: matrix in which each column contains the x-coordinates
            %   of the contour line of a cell in the Laguerre diagram
            %   y: matrix in which each column contains the y-coordinates
            %   of the contour line of a cell in the Laguerre diagram
            %   rho: matrix in which each column contains the distance from
            %   the anchor point in each cell in the Laguerre diagram

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);
            angle_num = length(obj.OT.Prep.angle_list);

            rho = zeros(angle_num, atom_num);

            hyp_a = (weights - weights') / 2;

            for atom_id = 1:atom_num
                rho(:, atom_id) ...
                    = hyperbolic_Dirichlet_tessellation_polar( ...
                    obj.OT.Prep.angle_list, ...
                    hyp_a(:, atom_id), ...
                    obj.OT.Prep.hyp.c(:, atom_id), ...
                    obj.OT.Prep.hyp.intercept(:, atom_id), ...
                    obj.OT.Prep.lin{atom_id}.c, ...
                    obj.OT.Prep.lin{atom_id}.intercept);
            end

            x = obj.OT.DiscMeas.Atoms(:, 1)' ...
                + rho .* cos(obj.OT.Prep.angle_list);
            y = obj.OT.DiscMeas.Atoms(:, 2)' ...
                + rho .* sin(obj.OT.Prep.angle_list);
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
            inside = reshape(all(reshape(sparse(coef_mat(:) > -1e-12), ...
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

    methods(Abstract, Access = protected)
        % Check if the input points are inside the support of the
        % probability measure
        inside = doCheckIfInsideSupport(obj, pts);

        % Generate random samples from a proposal distribution for
        % rejection sampling
        raw_samps = sampleFromProposalDistribution(obj, samp_num, ...
            rand_stream);

        % Evaluates inner integrals along a list of angles in the polar
        % coordinate system. The integrands are r -> r * p(r) and
        % r -> r^2 * p(r) where p(r) denotes the probability density
        % function.
        [val1, val2] = computeInnerIntegral(obj, atom_id, upper_limits)

        % Generate random samples from a proposal distribution for
        % conditional rejection sampling. For this procedure, the proposal
        % distribution is uniform on a circle.
        raw_samps = sampleFromCondProposalDistribution(obj, atom_id, ...
            samp_num, rand_stream);

        % Compute the acceptance probability in conditional rejection
        % sampling from a cell
        accept_probs = computeCondAcceptProb(obj, atom_id, raw_samps);

        % Function that initializes the test functions with respect to a 
        % simplicial cover and computes the integrals of the test functions 
        % with respect to the probability measure
        integrals = doSetSimplicialTestFuncs(obj, vargin);

        % Post-processing after computing semi-discrete optimal transport 
        % including computing the contour of each cell in the Laguerre 
        % diagram and preparing for conditional rejection sampling
        postProcessOptimalTransport(obj);
    end
end

