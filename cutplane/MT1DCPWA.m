classdef (Abstract) MT1DCPWA < LSIPMinCuttingPlaneAlgo
    % Abstract class for matching for teams problems with one-dimensional 
    % marginals, a quality space which is a polytope, and cost functions
    % that have the form c_i(x_i, z) = l_i(x_i - s' * z), where l_i is a
    % continuous piece-wise affine (CPWA) function

    properties(Constant)
        % numerical tolerance for deciding whether a point is inside a
        % polytope
        INSIDE_TOLERANCE = 1e-12;
    end

    properties(GetAccess = public, SetAccess = protected)
        % cell array containing the marginals
        Marginals;

        % cell array containing the cost functions
        CostFuncs;

        % struct containing information about the quality space
        Quality;
    end

    methods(Access = public)
        function obj = MT1DCPWA(marginals, costfuncs, quality, varargin)
            % Constructor method
            % Inputs:
            %   marginals: cell array containing objects of 
            %   ProbMeas1D_Interval
            %   costfuncs: cell array containing structs with fields
            %   weights, knots, and values
            %   quality: struct with fields:
            %       dim: the dimensionality of the quality space
            %       aux_num: number of auxiliary variables
            %       ineq_A: sparse matrix containing inequality constraints
            %       (A * x <= rhs)
            %       ineq_rhs: vector containing the right-hand side of
            %       inequality constraints
            %       eq_A: sparse matrix containing equality constraints 
            %       (A * x = rhs)
            %       eq_rhs: vector containing the right-hand side of
            %       equality constraints
            %       lb: vector containing lower bounds of the variables
            %       (including auxiliary variables)
            %       ub: vector containing upper bounds of the variables
            %       (including auxiliary variables)

            obj@LSIPMinCuttingPlaneAlgo(varargin{:});

            % set the default options for reducing constraints
            if ~isfield(obj.Options, 'reduce') ...
                    || isempty(obj.Options.reduce)
                obj.Options.reduce = struct;
            end

            if ~isfield(obj.Options.reduce, 'thres') ...
                    || isempty(obj.Options.reduce.thres)
                obj.Options.reduce.thres = inf;
            end

            if ~isfield(obj.Options.reduce, 'freq') ...
                    || isempty(obj.Options.reduce.freq)
                obj.Options.reduce.freq = 20;
            end

            if ~isfield(obj.Options.reduce, 'max_iter') ...
                    || isempty(obj.Options.reduce.max_iter)
                obj.Options.reduce.max_iter = inf;
            end

            % set the marginals
            obj.Marginals = marginals;

            % set the quality space
            obj.Quality = struct;
            obj.Quality.Dim = quality.dim;
            obj.Quality.AuxVarNum = quality.aux_num;

            % build the representation of the quality space
            if ~isfield(quality, 'eq_A') || isempty(quality.eq_A) ...
                    || ~isfield(quality, 'eq_rhs') ...
                    || isempty(quality.eq_rhs)
                quality.eq_A = sparse(0, quality.dim + quality.aux_num);
                quality.eq_rhs = zeros(0, 1);
            end

            if ~isfield(quality, 'ineq_A') || isempty(quality.ineq_A) ...
                    || ~isfield(quality, 'ineq_rhs') ...
                    || isempty(quality.ineq_rhs)
                quality.ineq_A = sparse(0, quality.dim + quality.aux_num);
                quality.ineq_rhs = zeros(0, 1);
            end

            if ~isfield(quality, 'lb') || isempty(quality.lb)
                quality.lb = -inf(quality.dim + quality.aux_num, 1);
            end

            if ~isfield(quality, 'ub') || isempty(quality.ub)
                quality.ub = inf(quality.dim + quality.aux_num, 1);
            end

            ineq_num = length(quality.ineq_rhs);
            eq_num = length(quality.eq_rhs);

            obj.Quality.GurobiModel = struct;
            obj.Quality.GurobiModel.A = sparse([quality.eq_A; ...
                quality.ineq_A]);
            obj.Quality.GurobiModel.rhs = [quality.eq_rhs; ...
                quality.ineq_rhs];
            obj.Quality.GurobiModel.sense = [repmat('=', eq_num, 1); ...
                repmat('<', ineq_num, 1)];
            obj.Quality.GurobiModel.lb = quality.lb;
            obj.Quality.GurobiModel.ub = quality.ub;

            % check if the quality space is bounded
            check_options = struct;
            check_options.OutputFlag = 0;
            for var_id = 1:obj.Quality.Dim
                for dir_sign = [-1, 1]
                    check_model = obj.Quality.GurobiModel;
                    check_model.modelsense = 'min';
                    check_model.obj = zeros(obj.Quality.Dim ...
                        + obj.Quality.AuxVarNum, 1);
                    check_model.obj(var_id) = dir_sign;

                    check_result = gurobi(check_model, check_options);

                    if ~strcmp(check_result.status, 'OPTIMAL')
                        error('the quality space is empty or unbounded');
                    end
                end
            end

            % set the cost functions
            marg_num = length(marginals);
            obj.CostFuncs = cell(marg_num, 1);

            obj.Storage.MarginalLowerBounds = zeros(marg_num, 1);
            obj.Storage.MarginalUpperBounds = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                cf = costfuncs{marg_id};
                obj.CostFuncs{marg_id} = struct;
                obj.CostFuncs{marg_id}.Weights = cf.weights;
                obj.CostFuncs{marg_id}.Knots = cf.knots;
                obj.CostFuncs{marg_id}.KnotDiff = diff(cf.knots);

                % check the inputs 
                assert(length(cf.weights) == obj.Quality.Dim, ...
                    'weights in the cost function are mis-specified');
                assert(all(obj.CostFuncs{marg_id}.KnotDiff > 0), ...
                    ['knots in the cost function are not ' ...
                    'in ascending order']);

                obj.CostFuncs{marg_id}.Values = cf.values;
                obj.CostFuncs{marg_id}.ValueDiff = diff(cf.values);

                % Lipschitz constant wrt change in the first argument (the
                % input from the support of the marginal)
                obj.CostFuncs{marg_id}.LipschitzConstMarg = max(abs( ...
                    obj.CostFuncs{marg_id}.ValueDiff ...
                    ./ obj.CostFuncs{marg_id}.KnotDiff));

                % Lipschitz constant wrt change in the second argument (the
                % input from the quality space)
                obj.CostFuncs{marg_id}.LipschitzConstQuality = ...
                    obj.CostFuncs{marg_id}.LipschitzConstMarg ...
                    * sqrt(cf.weights' * cf.weights);

                % check the domain of the CPWA function
                check_model = obj.Quality.GurobiModel;
                check_model.obj = [cf.weights; ...
                    zeros(obj.Quality.AuxVarNum, 1)];
                check_model.modelsense = 'min';

                check_result = gurobi(check_model, check_options);

                if ~strcmp(check_result.status, 'OPTIMAL')
                    error('unexpected error');
                end

                assert(obj.CostFuncs{marg_id}.Knots(end) >= ...
                    obj.Marginals{marg_id}.Supp.UpperBound ...
                    - check_result.objval, ...
                    'domain of the CPWA function is too small');

                check_model.modelsense = 'max';

                check_result = gurobi(check_model, check_options);

                if ~strcmp(check_result.status, 'OPTIMAL')
                    error('unexpected error');
                end

                assert(obj.CostFuncs{marg_id}.Knots(1) <= ...
                    obj.Marginals{marg_id}.Supp.LowerBound ...
                    - check_result.objval, ...
                    'domain of the CPWA function is too small');

                obj.Storage.MarginalLowerBounds(marg_id) = ...
                    obj.Marginals{marg_id}.Supp.LowerBound;
                obj.Storage.MarginalUpperBounds(marg_id) = ...
                    obj.Marginals{marg_id}.Supp.UpperBound;
            end

            % set up the gurobi model for computing the minimum cost over
            % the quality space
            min_cost_model = struct;
            lb_cell = cell(marg_num, 1);
            ub_cell = cell(marg_num, 1);
            vtype_cell = cell(marg_num, 1);
            obj_cell = cell(marg_num, 1);
            obj_con = 0;

            % compute the indices corresponding to the variables
            quality_indices = (1:(obj.Quality.Dim ...
                + obj.Quality.AuxVarNum))';
            
            % variables between 0 and 1 that are used to interpolate within
            % the intervals
            zeta_indices = cell(marg_num, 1);

            % binary-valued variables for indicating the interval
            iota_indices = cell(marg_num, 1);

            var_counter = quality_indices(end);
            
            for marg_id = 1:marg_num
                costfunc_knot_num = length(obj.CostFuncs{marg_id}.Knots);
                costfunc_vals = obj.CostFuncs{marg_id}.Values;

                zeta_indices{marg_id} = var_counter ...
                    + (1:costfunc_knot_num - 1)';
                var_counter = var_counter + costfunc_knot_num - 1;
                iota_indices{marg_id} = var_counter ...
                    + (1:costfunc_knot_num - 2)';
                var_counter = var_counter + costfunc_knot_num - 2;

                lb_cell{marg_id} = [zeros(costfunc_knot_num - 1, 1); ...
                    -inf(costfunc_knot_num - 2, 1)];
                ub_cell{marg_id} = [1; inf(costfunc_knot_num - 2, 1); ...
                    inf(costfunc_knot_num - 2, 1)];
                vtype_cell{marg_id} = [repmat('C', ...
                    costfunc_knot_num - 1, 1); ...
                    repmat('B', costfunc_knot_num - 2, 1)];

                obj_con = obj_con + costfunc_vals(1);
                obj_cell{marg_id} = [diff(costfunc_vals); ...
                    zeros(costfunc_knot_num - 2, 1)];
            end

            quality_model = obj.Quality.GurobiModel;
            min_cost_model.objcon = obj_con;
            min_cost_model.obj = [zeros(obj.Quality.Dim ...
                + obj.Quality.AuxVarNum, 1); ...
                vertcat(obj_cell{:})];
            min_cost_model.lb = [quality_model.lb; vertcat(lb_cell{:})];
            min_cost_model.ub = [quality_model.ub; vertcat(ub_cell{:})];
            min_cost_model.vtype = [repmat('C', obj.Quality.Dim ...
                + obj.Quality.AuxVarNum, 1); ...
                vertcat(vtype_cell{:})];

            rhs_input_indices = zeros(marg_num, 1);
            A_cell = cell(marg_num, 1);
            rhs_cell = cell(marg_num, 1);
            sense_cell = cell(marg_num, 1);
            constr_counter = 0;

            for marg_id = 1:marg_num
                costfunc_knots = obj.CostFuncs{marg_id}.Knots;
                costfunc_knot_num = length(costfunc_knots);

                % constraint: iota_{i,t} - zeta_{i,t} <= 0 and
                % constraint: iota_{i,t} - zeta_{i, t+1} >= 0
                A_ordering_r = ...
                    repelem((1:(costfunc_knot_num - 2) * 2)', 2, 1);
                A_ordering_c_mat = [iota_indices{marg_id}'; ...
                    zeta_indices{marg_id}(1:end - 1)'; ...
                    iota_indices{marg_id}'; ...
                    zeta_indices{marg_id}(2:end)'];
                A_ordering_c = A_ordering_c_mat(:);
                A_ordering_v = repmat([1; -1; 1; -1], ...
                    costfunc_knot_num - 2, 1);
                A_ordering = sparse(A_ordering_r, A_ordering_c, ...
                    A_ordering_v, 2 * (costfunc_knot_num - 2), ...
                    var_counter);
                rhs_ordering = zeros(2 * (costfunc_knot_num - 2), 1);
                sense_ordering = repmat(['<'; '>'], ...
                    costfunc_knot_num - 2, 1);

                % constraint linking zeta and the inputs
                A_link = sparse(ones(costfunc_knot_num - 1 ...
                    + obj.Quality.Dim, 1), ...
                    [zeta_indices{marg_id}; ...
                    quality_indices(1:obj.Quality.Dim)], ...
                    [diff(costfunc_knots); ...
                    obj.CostFuncs{marg_id}.Weights], ...
                    1, var_counter);
                rhs_link = -costfunc_knots(1);
                sense_link = '=';

                rhs_input_indices(marg_id) = constr_counter ...
                    + 2 * costfunc_knot_num - 3;
                constr_counter = rhs_input_indices(marg_id);

                A_cell{marg_id} = [A_ordering; A_link];
                rhs_cell{marg_id} = [rhs_ordering; rhs_link];
                sense_cell{marg_id} = [sense_ordering; sense_link];
            end

            min_cost_model.A = [vertcat(A_cell{:}); ...
                [quality_model.A, sparse(size(quality_model.A, 1), ...
                var_counter - obj.Quality.Dim - obj.Quality.AuxVarNum)]];
            min_cost_model.rhs = [vertcat(rhs_cell{:}); quality_model.rhs];
            min_cost_model.sense = [vertcat(sense_cell{:}); ...
                quality_model.sense];
            min_cost_model.rhs_input_indices = rhs_input_indices;
            min_cost_model.quality_indices = (1:obj.Quality.Dim)';

            obj.Quality.MinCostGurobiModel = min_cost_model;
        end

        function inside = checkIfInsideQualitySpace(obj, pts, batch_size)
            % Check if the points are inside the quality space. Compute in 
            % batches if necessary to avoid excessive memory use.
            % Inputs:
            %   pts: matrix containing the input points, where each row
            %   corresponds to an input
            %   batch_size: maximum number of inputs to be evaluated
            %   together via a vectorized routine (default is 10300)
            % Output:
            %   inside: boolean vector indicating whether each of the
            %   points is inside the quality space

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 10300;
            end

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            inside_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                [inside_cell{batch_id}] ...
                    = obj.doCheckIfInsideQualitySpace( ...
                    pts(((batch_id - 1) * batch_size + 1):min(batch_id ...
                    * batch_size, input_num), :));
            end

            inside = vertcat(inside_cell{:});
        end

        function vals = evaluateCostFunc(obj, marg_id, marg_inputs, ...
                quality_inputs, batch_size)
            % Evaluate the cost function between a marginal and the quality
            % space. Computation is done in batches if necessary. Warning:
            % this function does not check if the inputs are inside their
            % corresponding domains. 
            % Inputs: 
            %   marg_id: index of the selected marginal
            %   marg_inputs: vector containing inputs in the support of the
            %   marginal
            %   quality_inputs: matrix where each row is an input in the
            %   quality space
            %   batch_size: number of inputs to be handled together
            %   (default is 1e4)
            % Output:
            %   vals: vector containing the computed cost values

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            input_num = length(marg_inputs);
            assert(size(quality_inputs, 1) == input_num, ...
                'input mis-specified');

            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                cf = obj.CostFuncs{marg_id};
                batch_list = ((batch_id - 1) * batch_size ...
                    + 1):min(batch_id * batch_size, input_num);
                batch_CPWA_inputs = marg_inputs(batch_list) ...
                    - quality_inputs(batch_list, :) * cf.Weights;
                vals_cell{batch_id} = sum(min(max((batch_CPWA_inputs ...
                    - cf.Knots(1:end - 1)') ./ cf.KnotDiff', 0), 1) ...
                    .* cf.ValueDiff', 2) + cf.Values(1);
            end

            vals = vertcat(vals_cell{:});
        end

        function vals = evaluateSumOfCostFuncs(obj, marg_input_mat, ...
                quality_inputs, batch_size)
            % Evaluate the sum of cost functions. Computation is done in 
            % batches if necessary. Warning: this function does not check 
            % if the inputs are inside their corresponding domains. 
            % Inputs: 
            %   marg_input_mat: matrix where each row is an input and each
            %   column corresponds to a marginal
            %   quality_inputs: matrix where each row is an input in the
            %   quality space
            %   batch_size: number of inputs to be handled together
            %   (default is 1e4)
            % Output:
            %   vals: vector containing the computed cost values

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            marg_num = length(obj.Marginals);
            input_num = size(marg_input_mat, 1);
            assert(size(quality_inputs, 1) == input_num, ...
                'input mis-specified');
            assert(size(marg_input_mat, 2) == marg_num, ...
                'input mis-specified');

            vals = zeros(input_num, 1);

            for marg_id = 1:marg_num
                vals = vals + obj.evaluateCostFunc(marg_id, ...
                    marg_input_mat(:, marg_id), quality_inputs, ...
                    batch_size);
            end
        end
        
        function setSimplicialTestFuncs(obj, args_cell)
            % Set the simplicial test functions for all marginals at the
            % same time
            % Input:
            %   args_cell: cell array where each cell is a cell array
            %   containing all inputs to the method setSimplicialTestFuncs
            %   of each marginal

            for marg_id = 1:length(obj.Marginals)
                obj.Marginals{marg_id}.setSimplicialTestFuncs( ...
                    args_cell{marg_id}{:});
            end

            % after setting the simplicial functions, initialize the
            % quantities for the cutting-plane algorithm
            obj.initializeSimplicialTestFuncs();
        end

        function setLSIPSolutions(obj, primal_sol, dual_sol, ...
                LSIP_UB, LSIP_LB)
            % Set the primal and dual solution of the LSIP problem in the
            % runtime environment. This is to allow certain quantities to
            % be computed using stored versions of the primal and dual
            % solutions without executing the cutting-plane algorithm
            % again. 
            % Inputs: 
            %   primal_sol: struct containing information about the primal
            %   solution of the LSIP problem
            %   dual_sol: struct containing information about the dual
            %   solution of the LSIP problem
            %   LSIP_UB: the upper bound for the optimal value of the LSIP
            %   problem
            %   LSIP_LB: the lower bound for the optimal value of the LSIP
            %   problem

            if isempty(obj.Runtime)
                obj.Runtime = struct;
            end

            obj.Runtime.PrimalSolution = primal_sol;
            obj.Runtime.DualSolution = dual_sol;
            obj.Runtime.LSIP_UB = LSIP_UB;
            obj.Runtime.LSIP_LB = LSIP_LB;
        end

        function [UB_disc_list, samps] ...
                = getMTUpperBoundDiscreteWRepetition(obj, samp_num, ...
                rep_num, rand_stream, batch_size)
            % Compute the upper bound for the matching for teams problem
            % based on the discrete measure on the quality space with 
            % repetition of Monte Carlo integration.
            % Inputs:
            %   samp_num: number of samples for Monte Carlo integration
            %   rep_num: number of repetitions
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: number of samples to be handled simultaneously
            %   during the computation of the sum of the cost functions 
            %   (default is 1e4)
            % Output:
            %   UB_disc_list: vector containing the computed upper bounds 
            %   for the matching for teams problem based on the discrete 
            %   measure on the quality space
            %   samps: one particular set of Monte Carlo samples used in 
            %   the approximation of the bounds

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') ...
                    || isempty(batch_size)
                batch_size = 1e4;
            end

            UB_disc_list = zeros(rep_num, 1);

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, ...
                    ['--- Monte Carlo sampling (discrete measure) ' ...
                    'starts ---\n']);
            end

            for rep_id = 1:rep_num
                [UB_disc_list(rep_id), samps] ...
                    = obj.getMTUpperBoundDiscrete(samp_num, ...
                    rand_stream, batch_size);

                % display output
                if obj.Options.display
                    fprintf(['%s: ' ...
                        'Monte Carlo sampling (discrete measure) ' ...
                        'repetition %3d done\n'], class(obj), rep_id);
                end

                % write log
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, ['%s: ' ...
                        'Monte Carlo sampling (discrete measure) ' ...
                        'repetition %3d done\n'], class(obj), rep_id);
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, ...
                    ['--- Monte Carlo sampling (discrete measure) ' ...
                    'ends ---\n\n']);
                fclose(log_file);
            end
        end
    
        function [UB_disc_list, UB_cont_list, samps] ...
                = getMTUpperBoundsWRepetition(obj, samp_num, rep_num, ...
                rand_stream, batch_size)
            % Compute two upper bounds for the matching for teams problem 
            % with repetition of Monte Carlo integration.
            % Inputs:
            %   samp_num: number of samples for Monte Carlo integration
            %   rep_num: number of repetitions
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: number of samples to be handled simultaneously
            %   during the computation of the sum of the cost functions 
            %   (default is 1e4)
            % Output:
            %   UB_disc_list: vector containing the computed upper bounds 
            %   for the matching for teams problem based on the discrete 
            %   measure on the quality space
            %   UB_cont_list: vector containing the computed upper bounds 
            %   for the matching for teams problem based on the continuous 
            %   measure on the quality space
            %   samps: one particular set of Monte Carlo samples used in 
            %   the approximation of the bounds

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') ...
                    || isempty(batch_size)
                batch_size = 1e4;
            end

            UB_disc_list = zeros(rep_num, 1);
            UB_cont_list = zeros(rep_num, 1);

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, ...
                    '--- Monte Carlo sampling starts ---\n');
            end

            for rep_id = 1:rep_num
                [UB_disc_list(rep_id), UB_cont_list(rep_id), samps] ...
                    = obj.getMTUpperBounds(samp_num, rand_stream, ...
                    batch_size);

                % display output
                if obj.Options.display
                    fprintf(['%s: ' ...
                        'Monte Carlo sampling repetition %3d done\n'], ...
                        class(obj), rep_id);
                end

                % write log
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, ['%s: ' ...
                        'Monte Carlo sampling repetition %3d done\n'], ...
                        class(obj), rep_id);
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, ...
                    '--- Monte Carlo sampling ends ---\n\n');
                fclose(log_file);
            end
        end
    end

    methods(Abstract, Access = public)
        % Update the simplicial test functions after an execution of the 
        % cutting-plane algorithm. Besides setting the new simplicial test 
        % functions, a new coupling of discretized marginals is generated. 
        % This coupling can be used to generate initial constraints for the 
        % new LSIP problem with the updated test functions.
        [coup_indices, coup_probs, opt_quality] ...
                = updateSimplicialTestFuncs(obj, args_cell);

        % Initialize some quantities related to the simplicial test 
        % function of the marginals
        initializeSimplicialTestFuncs(obj);

        % Evaluate the parametric function on the support of a marginal
        % parametrized by simplicial test functions that is optimized
        % using the cutting plane algorithm. This corresponds to the
        % primal solution of the LSIP problem. The computation is done
        % is batches if necessary.
        vals = evaluateOptParametricFunc(obj, marg_id, pts, batch_size);
        
        % Evaluate the transfer functions resulted from optimized
        % parametric functions. These transfer functions is part of
        % approximate matching equilibria. The computation is done in
        % batches if necessary.
        vals_mat = evaluateOptTransferFuncs(obj, pts, ref_pt, batch_size);
        
        % Retrieve the computed lower bound for the matching for teams
        % problem
        LB = getMTLowerBound(obj);
        
        % Retrieve the optimized discrete quality distribution. This is
        % a part of an approximate matching equilibrium.
        distri = getOptDiscreteQualityDistr(obj);
        
        % Generate independent random sample from the coupling between the 
        % discrete measure on the quality space and a continuous marginal
        coup = randSampleFromOptDiscreteCoupling(obj, marg_id, ...
            samp_num, rand_stream);
        
        % Generate independent random sample from the continuous measure on 
        % the quality space
        ql_samps = randSampleFromOptContinuousQualityDistr(obj, ...
            samp_num, rand_stream);
        
        % Generate independent random sample from the coupling between the 
        % continuous measure on the quality space and a continuous marginal
        coup = randSampleFromOptContinuousCoupling(obj, marg_id, ...
            samp_num, rand_stream);

        % Compute the upper bound for the matching for teams problem based 
        % on the discrete measure on the quality space. The bounds are 
        % approximated by Monte Carlo integration.
        [UB_disc, UB_cont, samps] = getMTUpperBoundDiscrete(obj, ...
            samp_num, rand_stream);
        
        % Compute two upper bounds for the matching for teams problem based 
        % on the discrete measure on the quality space and based on the 
        % continuous measure on the quality space. The bounds are 
        % approximated by Monte Carlo integration.
        [UB_disc, UB_cont, samps] = getMTUpperBounds(obj, samp_num, ...
            rand_stream);

        % Compute the error bound for the objective value of the matching
        % for teams problem based on the Lipschitz constant of the cost
        % functions and the optimal transport distances between measures
        EB = getMTErrorBoundBasedOnOT(obj);

        % Compute the theoretical error bound for the objective value of
        % the matching for teams problem via the Lipschitz constant of the 
        % cost functions and the mesh sizes of the simplicial covers
        EB = getMTTheoreticalErrorBound(obj, LSIP_tolerance);
    end

    methods(Access = protected)
        function inside = doCheckIfInsideQualitySpace(obj, pts)
            % Check if the points are inside the quality space
            % Input:
            %   pts: matrix containing the input points, where each row
            %   corresponds to an input
            % Output:
            %   inside: boolean vector indicating whether each of the
            %   points is inside the quality space

            check_model = obj.Quality.GurobiModel;
            lb_list = all(pts >= check_model.lb(1:obj.Quality.Dim)' ...
                - MT1DCPWA.INSIDE_TOLERANCE, 2);
            ub_list = all(pts <= check_model.ub(1:obj.Quality.Dim)' ...
                + MT1DCPWA.INSIDE_TOLERANCE, 2);

            if obj.Quality.AuxVarNum == 0
                % if the quality space is characterized by inequality and
                % equality constraints only without any auxiliary 
                % variables, we can simply check whether these constraints 
                % are satisfied
                constr_diff =  pts * check_model.A' - check_model.rhs';

                ineq_list = all(constr_diff(:, check_model.sense ...
                    == '<') <= MT1DCPWA.INSIDE_TOLERANCE, 2);
                eq_list = all(abs(constr_diff(:, check_model.sense ...
                    == '=')) <= MT1DCPWA.INSIDE_TOLERANCE, 2);

                inside = lb_list & ub_list & ineq_list & eq_list;
            else
                % if there are auxiliary variables, we solve an LP problem
                % to determine whether the point is inside the quality
                % space
                inside = lb_list & ub_list;
                rem_num = sum(inside);
                rem_pts = pts(inside, :);
                rem_inside = false(rem_num, 1);

                aux_model = struct;
                aux_model.obj = zeros(obj.Quality.AuxVarNum, 1);
                rhs_mat = check_model.rhs ...
                    - check_model.A(:, 1:obj.Quality.Dim) * rem_pts';
                aux_model.A = check_model.A(:, obj.Quality.Dim + 1:end);
                aux_model.sense = check_model.sense;
                aux_model.lb = check_model.lb(obj.Quality.Dim + 1:end);
                aux_model.ub = check_model.ub(obj.Quality.Dim + 1:end);
                
                aux_options = struct;
                aux_options.OutputFlag = 0;

                for input_id = 1:rem_num
                    aux_model.rhs = rhs_mat(:, input_id);
                    aux_result = gurobi(aux_model, aux_options);

                    if strcmp(aux_result.status, 'OPTIMAL')
                        rem_inside(input_id) = true;
                    end
                end

                inside(inside) = rem_inside;
            end
        end

        function initializeBeforeRun(obj)
            % Initialize the algorithm by computing some static quantities

            if ~obj.Storage.SimplicialTestFuncsInitialized
                obj.initializeSimplicialTestFuncs();
            end
        end

        function bisect_cell = compute1DBisection(obj, intv_num)
            % Compute a bisection of knots in a one-dimensional continuous
            % piece-wise affine (CPWA) function used for formulating it
            % into a mixed-integer programming problem
            % Input:
            %   intv_num: number of intervals in the CPWA function
            % Output:
            %   bisect_cell: cell array containing the bisection, the
            %   number of rows is equal to ceil(log2(intv_num)) and the
            %   number of columns is 2

            bit_len = ceil(log2(intv_num));

            if ~isfield(obj.Storage, 'ReflexiveBinarySequences') ...
                || isempty(obj.Storage.ReflexiveBinarySequences)
                obj.Storage.ReflexiveBinarySequences = cell(20, 1);
            end

            if bit_len == 0
                bisect_cell = cell(0, 2);
                return;
            end

            if bit_len > 20
                error('the number of intervals is too large');
            end

            if isempty(obj.Storage.ReflexiveBinarySequences{bit_len})
                seq_mat = zeros(2 ^ bit_len, bit_len);
                
                % fill the first two rows
                seq_mat(1:2, end) = [0; 1];

                for bit_id = 2:bit_len
                    % flip the previous matrix upside down and append to
                    % the end while adding a column of 1's to the left
                    prev_row_num = 2 ^ (bit_id - 1);
                    seq_mat(prev_row_num + (1:prev_row_num), ...
                        end - bit_id + 1:end) ...
                        = [ones(prev_row_num, 1), flipud( ...
                        seq_mat(1:prev_row_num, end - bit_id + 2:end))];
                end

                obj.Storage.ReflexiveBinarySequences{bit_len} = seq_mat;
            end

            seq_mat = obj.Storage.ReflexiveBinarySequences{bit_len}( ...
                1:intv_num, :);

            bisect_cell = cell(bit_len, 2);

            for bit_id = 1:bit_len
                rflx_col = seq_mat(:, bit_id);
                bisect_cell{bit_id, 1} = ...
                    find([rflx_col(1) == 0; ...
                    rflx_col(1:end - 1) == 0 ...
                    & rflx_col(2:end) == 0; ...
                    rflx_col(end) == 0]);
                bisect_cell{bit_id, 2} = ...
                    find([rflx_col(1) == 1; ...
                    rflx_col(1:end - 1) == 1 ...
                    & rflx_col(2:end) == 1; ...
                    rflx_col(end) == 1]);
            end
        end

        function [min_pts, min_vals] = computeMinCostOverQuality(obj, ...
                input_mat)
            % Solve the minimization of the sum of cost functions where the
            % arguments in the support of marginals are given.
            % Inputs:
            %   input_mat: matrix where each row corresponds to an input 
            %   and each column corresponds to a marginal
            % Outputs:
            %   min_pts: the computed minimizers stored in a matrix where
            %   each row contains a minimizer
            %   min_vals: the computed minimal values

            input_num = size(input_mat, 1);
            min_pts = zeros(input_num, obj.Quality.Dim);
            min_vals = zeros(input_num, 1);

            model = obj.Quality.MinCostGurobiModel;
            MIP_options = struct;
            MIP_options.OutputFlag = 0;
            MIP_options.IntFeasTol = 1e-6;
            MIP_options.FeasibilityTol = 1e-8;
            MIP_options.OptimalityTol = 1e-8;
            MIP_options.NodefileStart = 2;
            MIP_options.MIPGap = 1e-4;
            MIP_options.MIPGapAbs = 1e-10;

            for input_id = 1:input_num
                model.rhs(model.rhs_input_indices) = ...
                    obj.Quality.MinCostGurobiModel.rhs( ...
                    model.rhs_input_indices) + input_mat(input_id, :)';
                result = gurobi(model, MIP_options);

                if ~strcmp(result.status, 'OPTIMAL')
                    error('MIP solver failed');
                end

                min_pts(input_id, :) = result.x(model.quality_indices)';
                min_vals(input_id) = result.objval;
            end
        end

        function updateLSIPUB(obj, min_lb, optimizers) %#ok<INUSD> 
            % Update the LSIP upper bound after each call to the global
            % minimization oracle
            % Inputs:
            %   min_lb: the lower bound for the global minimization problem
            %   optimizers: a set of approximate optimizers of the global
            %   minimization problem

            obj.Runtime.LSIP_UB = min(obj.Runtime.LSIP_UB, ...
                obj.Runtime.LSIP_LB - min_lb);
        end
    end
end

