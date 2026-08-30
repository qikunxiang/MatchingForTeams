classdef (Abstract) ProbMeas1DInterval < handle
    % Abstract class for probability measures supported on a one-dimensional interval
    
    properties(Constant)
        % numerical tolerance for deciding if a point is inside an interval
        INSIDE_TOLERANCE = 1e-14;
    end

    properties(SetAccess = protected, GetAccess = public)
        % struct storing information about the support of the probability measure
        Supp;

        % struct storing information about the density function of the probability measure
        Dens;

        % struct storing information about the optimal transport coupling with a discrete measure
        OT = struct;

        % struct storing information about test functions for this probability measure with respect to a simplicial cover
        SimplicialTestFuncs = [];
    end
    
    methods
        function obj = ProbMeas1DInterval(int_lb, int_ub)
            % Constructor which sets the supporting interval
            % Inputs:
            %   int_lb: the lower boundary of the interval
            %   int_ub: the upper boundary of the interval

            assert(int_ub > int_lb, 'the interval is empty');

            obj.Supp = struct;
            obj.Supp.LowerBound = int_lb;
            obj.Supp.UpperBound = int_ub;
        end
        
        function inside = checkIfInsideSupport(obj, pts)
            % Check if the input points are inside the support of the probability measure
            % Inputs:
            %   pts: vector containing the input points
            % Output: 
            %   inside: boolean vector indicating whether each point is inside the support

            inside = pts >= obj.Supp.LowerBound - ProbMeas1DInterval.INSIDE_TOLERANCE & pts <= obj.Supp.UpperBound ...
                + ProbMeas1DInterval.INSIDE_TOLERANCE;
        end

        function setCoupledDiscreteMeasure(obj, atoms, probs)
            % Couple this probability measure with a discrete measure via optimal transport
            % Inputs:
            %   atoms: vector containing the atoms in the discrete measure
            %   probs: vector containing the corresponding probabilities in the discrete measure

            obj.OT.DiscMeas = struct;

            assert(all(obj.checkIfInsideSupport(atoms)), 'some atoms are outside the support');

            obj.OT.DiscMeas.Atoms = atoms;
            probs = probs / sum(probs);
            obj.OT.DiscMeas.Probs = probs;
            
            % sort the atoms into ascending order
            [~, sorted_indices] = sort(atoms, 'ascend');
            sorted_probs = probs(sorted_indices);

            % record the ordering of the atoms
            [~, atoms_order] = sort(sorted_indices, 'ascend');

            % compute the cumulative probabilities (of atoms from left to right)
            cum_probs = cumsum(sorted_probs);

            obj.OT.DiscMeas.CumProbsSorted = cum_probs(atoms_order);

            obj.OT.Coupled = true;
        end

        function samps = randSample(obj, samp_num, rand_stream, varargin)
            % Rejection sampling the probability measure.
            % Inputs:
            %   samp_num: number of samples to genenerate
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream) 
            % Output:
            %   samps: vector containing the generated samples

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            % first generate samples from Uniform[0, 1]
            u = rand(rand_stream, samp_num, 1);

            % then apply the inverse CDF transform
            samps = obj.evaluateInverseCDF(u, varargin{:});
        end

        function samp_cell = conditionalRandSample(obj, atom_list, samp_num_list, rand_stream, varargin)
            % Randomly generate samples from the conditional distributions given the coupled atoms in the discrete measure
            % Inputs: 
            %   atom_list: vector containing the atom indices to be sampling from
            %   samp_num_list: vector containing number of samples coupled with each of the atoms in atom_list
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream) 
            % Output:
            %   samp_cell: cell array where each cell contains a vector containing the generated samples

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~isfield(obj.OT, 'Coupled') || ~obj.OT.Coupled
                error('must set the coupled discrete measure first');
            end

            cell_num = length(atom_list);
            
            assert(length(samp_num_list) == cell_num, 'number of atoms mismatch');

            atom_cum_probs = obj.OT.DiscMeas.CumProbsSorted;
            atom_probs = obj.OT.DiscMeas.Probs;
            samp_cell = cell(cell_num, 1);

            for cell_id = 1:cell_num
                chosen_atom_id = atom_list(cell_id);
                u = atom_cum_probs(chosen_atom_id) - atom_probs(chosen_atom_id) * rand(rand_stream, samp_num_list(cell_id), 1);
                samp_cell{cell_id} = obj.evaluateInverseCDF(u, varargin{:});
            end
        end

        function [samps, disc_samps] = randSampleCoupled(obj, samp_num, rand_stream, varargin)
            % Randomly generate samples from the coupling of the probability measure and the discrete measure
            % Inputs: 
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream) 
            % Output:
            %   samps: vector containing the generated samples
            %   disc_samps: vector containing the corresponding sampled atom indices in the discrete measure

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~isfield(obj.OT, 'Coupled') || ~obj.OT.Coupled
                error('must set the coupled discrete measure first');
            end

            atom_num = length(obj.OT.DiscMeas.Probs);

            % generate samples from the discrete measure
            disc_samps = randsample(rand_stream, atom_num, samp_num, true, obj.OT.DiscMeas.Probs);

            % count the number of samples for each cell
            samp_num_list = accumarray(disc_samps, ones(samp_num, 1), [atom_num, 1]);

            samp_cell = obj.conditionalRandSample((1:atom_num)', samp_num_list, rand_stream, varargin{:});
            
            samps = vertcat(samp_cell{:});

            % order the samples according the discrete samples
            [~, ordering] = sort(disc_samps, 'ascend');
            samps(ordering, :) = samps;
        end

        function vals = evaluateSimplicialTestFuncs(obj, pts, batch_size)
            % Evaluate the test functions at given locations; the input locations must be inside the support of the measure
            % Input:
            %   pts: vector containing the input points
            %   batch_size: maximum number of inputs to be evaluated together via a vectorized routine (default is 1e4)
            % Output:
            %   vals: sparse matrix containing the computed function values where each row corresponds to an input and each column
            %   corresponds to a test function

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            assert(all(obj.checkIfInsideSupport(pts)), 'some points are not inside the support');

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                vals_cell{batch_id} = obj.doEvaluateSimplicialTestFuncs(pts(((batch_id - 1) * batch_size + 1):min(batch_id ...
                    * batch_size, input_num), :));
            end

            vals = vertcat(vals_cell{:});
        end

        function vals = evaluateWeightedSumOfSimplicialTestFuncs(obj, pts, coefficients, batch_size)
            % Evaluate a weighted sum of simplicial test functions at given locations; the input locations must be inside the support 
            % of the measure 
            % Input:
            %   pts: vector containing the input points
            %   coefficients: vector containing the coefficients of the
            %   test functions in the weighted sum
            %   batch_size: maximum number of inputs to be evaluated together via a vectorized routine (default is 1e4)
            % Output:
            %   vals: sparse matrix containing the computed function values where each row corresponds to an input and each column
            %   corresponds to a test function

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            assert(all(obj.checkIfInsideSupport(pts)), 'some points are not inside the support');

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                raw_vals = obj.doEvaluateSimplicialTestFuncs(pts(((batch_id - 1) * batch_size + 1):min(batch_id * batch_size, ...
                    input_num), :));
                vals_cell{batch_id} = raw_vals * coefficients;
            end

            vals = vertcat(vals_cell{:});
        end

        function testSimplicialTestFuncs(obj, grid_pts)
            % Test the computation of the integrals of test functions with respect to the probability measure as well as the evaluation
            % of test functions
            % Input:
            %   grid_pts: vector containing the points used for the test

            % approximate the probability measure with a discrete measure supported on the grid
            grid_density = obj.densityFunction(grid_pts);
            grid_density = grid_density / sum(grid_density);

            % only retain points that are inside the support
            inside_support = grid_density > 0;
            grid_pts = grid_pts(inside_support, :);
            grid_density = grid_density(inside_support);

            testfunc_vals = obj.evaluateSimplicialTestFuncs(grid_pts);

            integrals_grid = sum(testfunc_vals .* grid_density, 1)';

            fprintf('integrals:\n');
            display(full([obj.SimplicialTestFuncs.Integrals, integrals_grid]));

            fprintf('max diff. in integrals: %e\n', ...
                max(abs(obj.SimplicialTestFuncs.Integrals - integrals_grid)));
        end
    
        function integrals = setSimplicialTestFuncs(obj, knots)
            % Initialize the test functions with respect to a simplicial cover and compute the integrals of the test functions with 
            % respect to the probability measure
            % Input: 
            %   knots: vector containing the knots in the simplicial test functions; the first and last knots must agree with
            %   obj.Supp.LowerBound and obj.Supp.UpperBound

            assert(abs(knots(1) - obj.Supp.LowerBound) < ProbMeas1DInterval.INSIDE_TOLERANCE...
                && abs(knots(end) - obj.Supp.UpperBound) < ProbMeas1DInterval.INSIDE_TOLERANCE, ...
                'the first and the last knots are different from the end points of the support of the measure');

            knot_diff = diff(knots);

            assert(all(knot_diff > 0), 'knots are not in ascending order');

            integrals = obj.integrateSimplicialTestFuncs(knots);

            obj.SimplicialTestFuncs = struct;
            obj.SimplicialTestFuncs.Knots = knots;
            obj.SimplicialTestFuncs.Integrals = integrals;
            obj.SimplicialTestFuncs.KnotDiff = diff(obj.SimplicialTestFuncs.Knots);
            obj.SimplicialTestFuncs.MeshSize = max(obj.SimplicialTestFuncs.KnotDiff);
        end
    end

    methods(Abstract, Access = public)
        % Compute the probability density function
        dens = densityFunction(obj, x);

        % Evaluate the inverse cumulative density function
        x = evaluateInverseCDF(obj, p, varargin);

        % Compute the optimal transport cost to the coupled discrete measure
        cost = computeOTCost(obj);
    end

    methods(Access = protected)
        function vals = doEvaluateSimplicialTestFuncs(obj, pts)
            % Function that actually evaluates the test functions at the given locations 
            % Input: 
            %   pts: vector containing the input points
            % Output:
            %   vals: sparse matrix containing the computed function values where each row corresponds to an input and each column
            %   corresponds to a test function

            testfuncs = obj.SimplicialTestFuncs;
            knots = testfuncs.Knots;
            knot_diff = testfuncs.KnotDiff;
            knot_num = length(knots);
            input_num = length(pts);

            % compute the sub-intervals the inputs are in
            int_indices = min(sum(pts - knots' >= -ProbMeas1DInterval.INSIDE_TOLERANCE, 2), knot_num - 1);
            
            % compute the test function corresponding to the left end point
            % of the sub-interval
            int_lengths = knot_diff(int_indices);
            testfunc_right = (pts - knots(int_indices)) ./ int_lengths;

            % assemble the sparse matrix
            vals = sparse(repmat((1:input_num)', 2, 1), [int_indices; int_indices + 1], [1 - testfunc_right; testfunc_right], ...
                input_num, knot_num);
        end
    end

    methods(Abstract, Access = protected)
        % Compute the integrals of the test functions with respect to the probability measure 
        integrals = integrateSimplicialTestFuncs(obj, knots);
    end
end

