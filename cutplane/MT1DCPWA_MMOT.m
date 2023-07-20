classdef MT1DCPWA_MMOT < MT1DCPWA
    % Class for matching for teams problems with one-dimensional marginals,
    % a quality space which is a polytope, and cost functions that have the
    % form c_i(x_i, z) = l_i(x_i - s' * z), where l_i is a continuous 
    % piece-wise affine (CPWA) function. The problem is solved via the 
    % multi-marginal optimal transport (MMOT) formulation.

    methods(Access = public)
        function obj = MT1DCPWA_MMOT(marginals, costfuncs, quality, ...
                varargin)
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

            obj@MT1DCPWA(marginals, costfuncs, quality, varargin{:});

            % set the default option for the formulation of global
            % minimization problem (options are 'LOG' and 'INC')
            if ~isfield(obj.Options, 'global_formulation') ...
                    || isempty(obj.Options.global_formulation)
                obj.Options.global_formulation = 'LOG';
            end

            % set the default option for sanitizing the entries in the
            % inequality constraints
            if ~isfield(obj.Options, 'sanitation_threshold') ...
                    || isempty(obj.Options.sanitation_threshold)
                obj.Options.sanitation_threshold = 0;
            end

            % find a point in the quality space as a reference point
            model = obj.Quality.GurobiModel;
            model.obj = zeros(obj.Quality.Dim + obj.Quality.AuxVarNum, 1);
            options = struct;
            options.OutputFlag = 0;
            result = gurobi(model, options);
            
            if ~strcmp(result.status, 'OPTIMAL')
                error('unexpected error with the quality model');
            end

            obj.Quality.SamplePoint = result.x(1:obj.Quality.Dim);

            % this flag is used to track if the function
            % obj.initializeSimplicialTestFuncs has been called
            obj.Storage.SimplicialTestFuncsInitialized = false;

            marg_num = length(obj.Marginals);
            marg_testfunc_set = false(marg_num, 1);

            for marg_id = 1:marg_num
                marg_testfunc_set(marg_id) = ~isempty( ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs);
            end

            if all(marg_testfunc_set)
                % initialize the simplicial test functions at the end of 
                % the constructor if they have already been set
                obj.initializeSimplicialTestFuncs();
            end
        end

        function [coup_indices, coup_qualities, coup_costs, coup_probs] ...
                = updateSimplicialTestFuncs(obj, args_cell)
            % Update the simplicial test functions after an execution of
            % the cutting-plane algorithm. Besides setting the new
            % simplicial test functions, a new coupling of discretized
            % marginals is generated via reassembly of the dual solution
            % from the cutting-plane algorithm with the new discretized
            % marginals. This coupling can be used to generate initial
            % constraints for the new LSIP problem with the updated test
            % functions.
            % Input: 
            %   args_cell: cell array where each cell is a cell array
            %   containing all inputs to the methods setSimplicialTestFuncs
            %   of each marginal
            % Outputs:
            %   coup_indices: the indices of knots in the coupled discrete
            %   measure
            %   coup_qualities: 
            %   coup_costs: 
            %   coup_probs: the probabilities of atoms in the coupled
            %   discrete measure

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            marg_num = length(obj.Marginals);

            % retrieve the dual solution resulted from the cutting-plane
            % algorithm
            old_coup_atoms = obj.Runtime.DualSolution.InputAtoms;
            old_coup_probs = obj.Runtime.DualSolution.Probabilities;

            % set the new simplicial test functions
            obj.setSimplicialTestFuncs(args_cell);

            % the optimal couplings between the original discrete marginals
            % and the new discrete marginals where the cost function is
            % given by the Euclidean distance
            new_marg_atoms_cell = cell(marg_num, 1);
            new_marg_probs_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                new_marg_atoms_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Knots;
                new_marg_probs_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;
                
            end

            % perform discrete reassembly to get the new coupling
            [coup_indices, coup_probs] ...
                = discrete_reassembly_1D(old_coup_atoms, ...
                old_coup_probs, new_marg_atoms_cell, new_marg_probs_cell);

            % transform knot indices into actual coordinates
            coup_inputs = zeros(length(coup_probs), marg_num);

            for marg_id = 1:marg_num
                coup_inputs(:, marg_id) = ...
                    new_marg_atoms_cell{marg_id}(coup_indices(:, marg_id));
            end

            % solve for the optimizers in the quality space
            [coup_qualities, coup_costs] = ...
                obj.computeMinCostOverQuality(coup_inputs);
        end

        function initializeSimplicialTestFuncs(obj)
            % Initialize some quantities related to the simplicial test 
            % functions of the marginals

            marg_num = length(obj.Marginals);

            obj.Storage.MargKnotNumList = zeros(marg_num, 1);
            deci_logical_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                testfunc = obj.Marginals{marg_id}.SimplicialTestFuncs;
                testfunc_knot_num = length(testfunc.Knots);
                obj.Storage.MargKnotNumList(marg_id) = testfunc_knot_num;

                % the coefficient corresponding to the first test function
                % will not be included in the decision variable of the LP 
                % problem for identification purposes
                deci_logical_cell{marg_id} = [0; ...
                    ones(testfunc_knot_num - 1, 1)];
            end

            obj.Storage.TotalKnotNum = sum(obj.Storage.MargKnotNumList);

            % compute the offset of the knots in the marginals in the
            % vector containing all knots
            knot_num_cumsum = cumsum(obj.Storage.MargKnotNumList);
            obj.Storage.MargKnotNumOffsets = [0; ...
                knot_num_cumsum(1:end - 1)];

            % store the indices to place the decision variables of the LP 
            % problem in a vector containing the coefficients of the test 
            % functions
            obj.Storage.DeciVarIndicesInTestFuncs = find( ...
                vertcat(deci_logical_cell{:}));

            if strcmp(obj.Options.global_formulation, 'INC')
                % build the gurobi models for the global minimization 
                % problems; the one-dimensional CPWA functions are modeled 
                % by the incremental (INC) formulation and the 
                % two-dimensional CPWA functions are modeled by the convex 
                % combination (CC) formulation
                gl_model = struct;
    
                lb_cell = cell(marg_num, 1);
                ub_cell = cell(marg_num, 1);
                vtype_cell = cell(marg_num, 1);
                obj_cell = cell(marg_num, 1);
                obj_con = 0;
    
                % compute the indices corresponding to the variables
                input_indices = (1:marg_num)';
                quality_indices = marg_num + (1:(obj.Quality.Dim ...
                    + obj.Quality.AuxVarNum))';
                
                % variables between 0 and 1 that are used to interpolate 
                % within the intervals for the cost functions
                zeta_indices = cell(marg_num, 1);
    
                % binary-valued variables for indicating the interval for 
                % the cost functions
                iota_indices = cell(marg_num, 1);
    
                % variables between 0 and 1 that are used to interpolate 
                % within the intervals for the parametric potential 
                % functions
                xi_indices = cell(marg_num, 1);
    
                % binary-valued variables for indicating the interval for 
                % the parametric potential functions
                eta_indices = cell(marg_num, 1);
    
                var_counter = quality_indices(end);
                
                for marg_id = 1:marg_num
                    costfunc_knot_num = ...
                        length(obj.CostFuncs{marg_id}.Knots);
                    costfunc_vals = obj.CostFuncs{marg_id}.Values;
    
                    testfunc = obj.Marginals{marg_id}.SimplicialTestFuncs;
                    testfunc_knot_num = length(testfunc.Knots);
    
                    zeta_indices{marg_id} = var_counter ...
                        + (1:costfunc_knot_num - 1)';
                    var_counter = var_counter + costfunc_knot_num - 1;
                    iota_indices{marg_id} = var_counter ...
                        + (1:costfunc_knot_num - 2)';
                    var_counter = var_counter + costfunc_knot_num - 2;
    
                    xi_indices{marg_id} = var_counter ...
                        + (1:testfunc_knot_num - 1)';
                    var_counter = var_counter + testfunc_knot_num - 1;
                    eta_indices{marg_id} = var_counter ...
                        + (1:testfunc_knot_num - 2)';
                    var_counter = var_counter + testfunc_knot_num - 2;
    
                    lb_cell{marg_id} = ...
                        [zeros(costfunc_knot_num - 1, 1); ...
                        -inf(costfunc_knot_num - 2, 1); ...
                        zeros(testfunc_knot_num - 1, 1); ...
                        -inf(testfunc_knot_num - 2, 1)];
                    ub_cell{marg_id} = ...
                        [1; inf(costfunc_knot_num - 2, 1); ...
                        inf(costfunc_knot_num - 2, 1); ...
                        1; inf(testfunc_knot_num - 2, 1); ...
                        inf(testfunc_knot_num - 2, 1)];
                    vtype_cell{marg_id} = ...
                        [repmat('C', costfunc_knot_num - 1, 1); ...
                        repmat('B', costfunc_knot_num - 2, 1); ...
                        repmat('C', testfunc_knot_num - 1, 1); ...
                        repmat('B', testfunc_knot_num - 2, 1)];
    
                    obj_con = obj_con + costfunc_vals(1);
                    obj_cell{marg_id} = ...
                        [diff(costfunc_vals); ...
                        zeros(costfunc_knot_num - 2, 1); ...
                        zeros(testfunc_knot_num - 1, 1); ...
                        zeros(testfunc_knot_num - 2, 1)];
                end
    
                quality_model = obj.Quality.GurobiModel;
                gl_model.modelsense = 'min';
                gl_model.objcon = obj_con;
                gl_model.obj = [zeros(marg_num + obj.Quality.Dim ...
                    + obj.Quality.AuxVarNum, 1); ...
                    vertcat(obj_cell{:})];
                gl_model.lb = [-inf(marg_num, 1); quality_model.lb; ...
                    vertcat(lb_cell{:})];
                gl_model.ub = [inf(marg_num, 1); quality_model.ub; ...
                    vertcat(ub_cell{:})];
                gl_model.vtype = [repmat('C', marg_num ...
                    + obj.Quality.Dim + obj.Quality.AuxVarNum, 1); ...
                    vertcat(vtype_cell{:})];
    
                A_cell = cell(marg_num, 1);
                rhs_cell = cell(marg_num, 1);
                sense_cell = cell(marg_num, 1);
    
                % construct the matrix for computing the vector of 
                % coefficients that will be placed into the objective 
                % vector of the mixed-integer model from a vector of 
                % coefficients
                obj_coef_mat_cell = cell(marg_num, 1);

                % construct the matrix and the intercept vector for 
                % directly extracting the values of the test functions for 
                % the marginal from the results of the solver
                tfcoef_extract_mat_cell = cell(marg_num, 1);
                tfcoef_extract_intercept_cell = cell(marg_num, 1);
    
                for marg_id = 1:marg_num
                    costfunc_knots = obj.CostFuncs{marg_id}.Knots;
                    costfunc_knot_num = length(costfunc_knots);
    
                    testfunc = obj.Marginals{marg_id}.SimplicialTestFuncs;
                    testfunc_knots = testfunc.Knots;
                    testfunc_knot_num = length(testfunc_knots);
    
                    % constraint: iota_{i,t} - zeta_{i,t} <= 0 and
                    % constraint: iota_{i,t} - zeta_{i, t+1} >= 0
                    A_cf_ord_r = ...
                        repelem((1:(costfunc_knot_num - 2) * 2)', 2, 1);
                    A_cf_ord_c_mat = [iota_indices{marg_id}'; ...
                        zeta_indices{marg_id}(1:end - 1)'; ...
                        iota_indices{marg_id}'; ...
                        zeta_indices{marg_id}(2:end)'];
                    A_cf_ord_c = A_cf_ord_c_mat(:);
                    A_cf_ord_v = repmat([1; -1; 1; -1], ...
                        costfunc_knot_num - 2, 1);
                    A_cf_ord = sparse(A_cf_ord_r, A_cf_ord_c, ...
                        A_cf_ord_v, 2 * (costfunc_knot_num - 2), ...
                        var_counter);
                    rhs_cf_ord = zeros(2 * (costfunc_knot_num - 2), 1);
                    sense_cf_ord = repmat(['<'; '>'], ...
                        costfunc_knot_num - 2, 1);
    
                    % constraint linking zeta and the inputs
                    A_cf_link = sparse(ones(costfunc_knot_num - 1 ...
                        + obj.Quality.Dim + 1, 1), ...
                        [zeta_indices{marg_id}; ...
                        quality_indices(1:obj.Quality.Dim); ...
                        input_indices(marg_id)], ...
                        [diff(costfunc_knots); ...
                        obj.CostFuncs{marg_id}.Weights; ...
                        -1], 1, var_counter);
                    rhs_cf_link = -costfunc_knots(1);
                    sense_cf_link = '=';
    
                    % constraint: eta_{i,t} - xi_{i,t} <= 0 and
                    % constraint: eta_{i,t} - xi_{i, t+1} >= 0
                    A_tf_ord_r = ...
                        repelem((1:(testfunc_knot_num - 2) * 2)', 2, 1);
                    A_tf_ord_c_mat = [eta_indices{marg_id}'; ...
                        xi_indices{marg_id}(1:end - 1)'; ...
                        eta_indices{marg_id}'; ...
                        xi_indices{marg_id}(2:end)'];
                    A_tf_ord_c = A_tf_ord_c_mat(:);
                    A_tf_ord_v = repmat([1; -1; 1; -1], ...
                        testfunc_knot_num - 2, 1);
                    A_tf_ord = sparse(A_tf_ord_r, A_tf_ord_c, ...
                        A_tf_ord_v, 2 * (testfunc_knot_num - 2), ...
                        var_counter);
                    rhs_tf_ord = zeros(2 * (testfunc_knot_num - 2), 1);
                    sense_tf_ord = repmat(['<'; '>'], ...
                        testfunc_knot_num - 2, 1);
    
                    % constraint linking xi and the inputs
                    A_tf_link = sparse( ...
                        ones(testfunc_knot_num - 1 + 1, 1), ...
                        [xi_indices{marg_id}; input_indices(marg_id)], ...
                        [diff(testfunc_knots); -1], ...
                        1, var_counter);
                    rhs_tf_link = -testfunc_knots(1);
                    sense_tf_link = '=';
    
                    A_cell{marg_id} = [A_cf_ord; A_cf_link; ...
                        A_tf_ord; A_tf_link];
                    rhs_cell{marg_id} = [rhs_cf_ord; rhs_cf_link; ...
                        rhs_tf_ord; rhs_tf_link];
                    sense_cell{marg_id} = [sense_cf_ord; sense_cf_link; ...
                        sense_tf_ord; sense_tf_link];
    
                    obj_coef_r = repelem((1:testfunc_knot_num - 1)', 2, 1);
                    obj_coef_c_mat = [1:testfunc_knot_num - 1; ...
                        2:testfunc_knot_num];
                    obj_coef_v = repmat([1; -1], testfunc_knot_num - 1, 1);
                    obj_coef_mat_cell{marg_id} = sparse( ...
                        obj_coef_r, obj_coef_c_mat(:), obj_coef_v, ...
                        testfunc_knot_num - 1, testfunc_knot_num);

                    tfcoef_extract_mat_cell{marg_id} = sparse( ...
                        [(1:testfunc_knot_num - 1)'; ...
                        (2:testfunc_knot_num)'], ...
                        repmat(xi_indices{marg_id}, 2, 1), ...
                        [-ones(testfunc_knot_num - 1, 1); ...
                        ones(testfunc_knot_num - 1, 1)], ...
                        testfunc_knot_num, var_counter);
                    tfcoef_extract_intercept_cell{marg_id} = ...
                        [1; zeros(testfunc_knot_num - 1, 1)];
                end
    
                gl_model.A = [vertcat(A_cell{:}); ...
                    [sparse(size(quality_model.A, 1), marg_num), ...
                    quality_model.A, ...
                    sparse(size(quality_model.A, 1), var_counter ...
                    - marg_num - obj.Quality.Dim ...
                    - obj.Quality.AuxVarNum)]];
                gl_model.rhs = [vertcat(rhs_cell{:}); quality_model.rhs];
                gl_model.sense = [vertcat(sense_cell{:}); ...
                    quality_model.sense];
                gl_model.input_indices = (1:marg_num)';
                gl_model.quality_indices = marg_num + (1:obj.Quality.Dim)';
    
                % store the indices in the objective vector where the
                % coefficients of the test functions will be placed
                gl_model.obj_coef_indices = vertcat(xi_indices{:});
                gl_model.obj_coef_mat = blkdiag(obj_coef_mat_cell{:});

                % matrix and intercept vector for directly extracting
                % the values of the test functions for the marginal
                % from the results of the solver
                gl_model.tfcoef_extract_mat = ...
                    vertcat(tfcoef_extract_mat_cell{:});
                gl_model.tfcoef_extract_intercept = ...
                    vertcat(tfcoef_extract_intercept_cell{:});
    
                obj.Storage.GlobalMinGurobiModel = gl_model;
            elseif strcmp(obj.Options.global_formulation, 'LOG')
                % build the gurobi models for the global minimization 
                % problems; the one-dimensional CPWA functions are modeled 
                % by the logarithmic convex combination (LOG) formulation
                gl_model = struct;
    
                lb_cell = cell(marg_num, 1);
                ub_cell = cell(marg_num, 1);
                vtype_cell = cell(marg_num, 1);
                obj_cell = cell(marg_num, 1);
                obj_con = 0;
    
                % compute the indices corresponding to the variables
                input_indices = (1:marg_num)';
                quality_indices = marg_num + (1:(obj.Quality.Dim ...
                    + obj.Quality.AuxVarNum))';
                
                % variables between 0 and 1 that are used to interpolate 
                % within the intervals for the cost functions
                zeta_indices = cell(marg_num, 1);
    
                % binary-valued variables for indicating the interval for 
                % the cost functions
                iota_indices = cell(marg_num, 1);
    
                % variables between 0 and 1 that are used to interpolate 
                % within the intervals for the parametric potential 
                % functions
                xi_indices = cell(marg_num, 1);
    
                % binary-valued variables for indicating the interval for 
                % the parametric potential functions
                eta_indices = cell(marg_num, 1);
    
                costfunc_bit_len = zeros(marg_num, 1);
                costfunc_bisect_cell = cell(marg_num, 1);
                testfunc_bit_len = zeros(marg_num, 1);
                testfunc_bisect_cell = cell(marg_num, 1);

                var_counter = quality_indices(end);
                
                for marg_id = 1:marg_num
                    costfunc_knot_num = ...
                        length(obj.CostFuncs{marg_id}.Knots);
                    costfunc_vals = obj.CostFuncs{marg_id}.Values;
    
                    testfunc = obj.Marginals{marg_id}.SimplicialTestFuncs;
                    testfunc_knot_num = length(testfunc.Knots);

                    % compute the bisections of the knots in the cost
                    % function
                    costfunc_bit_len(marg_id) = ...
                        ceil(log2(costfunc_knot_num - 1));
                    costfunc_bisect_cell{marg_id} = ...
                        obj.compute1DBisection(costfunc_knot_num - 1);

                    % compute the bisections in the knots of the test
                    % functions
                    testfunc_bit_len(marg_id) = ...
                        ceil(log2(testfunc_knot_num - 1));
                    testfunc_bisect_cell{marg_id} = ...
                        obj.compute1DBisection(testfunc_knot_num - 1);
    
                    zeta_indices{marg_id} = var_counter ...
                        + (1:costfunc_knot_num)';
                    var_counter = var_counter + costfunc_knot_num;
                    iota_indices{marg_id} = var_counter ...
                        + (1:costfunc_bit_len(marg_id))';
                    var_counter = var_counter + costfunc_bit_len(marg_id);
    
                    xi_indices{marg_id} = var_counter ...
                        + (1:testfunc_knot_num)';
                    var_counter = var_counter + testfunc_knot_num;
                    eta_indices{marg_id} = var_counter ...
                        + (1:testfunc_bit_len(marg_id))';
                    var_counter = var_counter + testfunc_bit_len(marg_id);
    
                    lb_cell{marg_id} = ...
                        [zeros(costfunc_knot_num, 1); ...
                        -inf(costfunc_bit_len(marg_id), 1); ...
                        zeros(testfunc_knot_num, 1); ...
                        -inf(testfunc_bit_len(marg_id), 1)];
                    ub_cell{marg_id} = ...
                        [inf(costfunc_knot_num, 1); ...
                        inf(costfunc_bit_len(marg_id), 1); ...
                        inf(testfunc_knot_num, 1); ...
                        inf(testfunc_bit_len(marg_id), 1)];
                    vtype_cell{marg_id} = ...
                        [repmat('C', costfunc_knot_num, 1); ...
                        repmat('B', costfunc_bit_len(marg_id), 1); ...
                        repmat('C', testfunc_knot_num, 1); ...
                        repmat('B', testfunc_bit_len(marg_id), 1)];
    
                    obj_cell{marg_id} = ...
                        [costfunc_vals; ...
                        zeros(costfunc_bit_len(marg_id), 1); ...
                        zeros(testfunc_knot_num, 1); ...
                        zeros(testfunc_bit_len(marg_id), 1)];
                end
    
                quality_model = obj.Quality.GurobiModel;
                gl_model.modelsense = 'min';
                gl_model.objcon = obj_con;
                gl_model.obj = [zeros(marg_num + obj.Quality.Dim ...
                    + obj.Quality.AuxVarNum, 1); ...
                    vertcat(obj_cell{:})];
                gl_model.lb = [-inf(marg_num, 1); quality_model.lb; ...
                    vertcat(lb_cell{:})];
                gl_model.ub = [inf(marg_num, 1); quality_model.ub; ...
                    vertcat(ub_cell{:})];
                gl_model.vtype = [repmat('C', marg_num ...
                    + obj.Quality.Dim + obj.Quality.AuxVarNum, 1); ...
                    vertcat(vtype_cell{:})];
    
                A_cell = cell(marg_num, 1);
                rhs_cell = cell(marg_num, 1);
                sense_cell = cell(marg_num, 1);
    
                % construct the matrix for computing the vector of 
                % coefficients that will be placed into the objective 
                % vector of the mixed-integer model from a vector of 
                % coefficients
                obj_coef_mat_cell = cell(marg_num, 1);

                % construct the matrix and the intercept vector for 
                % directly extracting the values of the test functions for 
                % the marginal from the results of the solver
                tfcoef_extract_mat_cell = cell(marg_num, 1);
                tfcoef_extract_intercept_cell = cell(marg_num, 1);
    
                for marg_id = 1:marg_num
                    costfunc_weights = obj.CostFuncs{marg_id}.Weights;
                    costfunc_knots = obj.CostFuncs{marg_id}.Knots;
                    costfunc_knot_num = length(costfunc_knots);
    
                    testfunc = obj.Marginals{marg_id}.SimplicialTestFuncs;
                    testfunc_knots = testfunc.Knots;
                    testfunc_knot_num = length(testfunc_knots);
    
                    % constraint that the sum of all zeta coefficients must
                    % be equal to 1
                    A_cf_sum = sparse( ...
                        ones(costfunc_knot_num, 1), ...
                        zeta_indices{marg_id}, ...
                        ones(costfunc_knot_num, 1), ...
                        1, var_counter);
                    rhs_cf_sum = 1;
                    sense_cf_sum = '=';

                    % constraints identifying non-zero zeta's with iota's
                    A_cf_id_r_cell = cell(costfunc_bit_len(marg_id), 1);
                    A_cf_id_c_cell = cell(costfunc_bit_len(marg_id), 1);
                    A_cf_id_v_cell = cell(costfunc_bit_len(marg_id), 1);

                    for bit_id = 1:costfunc_bit_len(marg_id)
                        cf_bit0_list = ...
                            costfunc_bisect_cell{marg_id}{bit_id, 1};
                        cf_bit0_num = length(cf_bit0_list);
                        cf_bit1_list = ...
                            costfunc_bisect_cell{marg_id}{bit_id, 2};
                        cf_bit1_num = length(cf_bit1_list);

                        A_cf_id_r_cell{bit_id} = (bit_id - 1) * 2 ...
                            + [ones(cf_bit0_num + 1, 1); ...
                            2 * ones(cf_bit1_num + 1, 1)];
                        A_cf_id_c_cell{bit_id} = ...
                            [zeta_indices{marg_id}(cf_bit0_list); ...
                            iota_indices{marg_id}(bit_id); ...
                            zeta_indices{marg_id}(cf_bit1_list); ...
                            iota_indices{marg_id}(bit_id)];
                        A_cf_id_v_cell{bit_id} = ...
                            [ones(cf_bit0_num, 1); -1; ...
                            ones(cf_bit1_num, 1); 1];
                    end

                    A_cf_id = sparse(vertcat(A_cf_id_r_cell{:}), ...
                        vertcat(A_cf_id_c_cell{:}), ...
                        vertcat(A_cf_id_v_cell{:}), ...
                        2 * costfunc_bit_len(marg_id), var_counter);
                    rhs_cf_id = repmat([0; 1], ...
                        costfunc_bit_len(marg_id), 1);
                    sense_cf_id = repmat('<', ...
                        2 * costfunc_bit_len(marg_id), 1);

                    % constraint linking zeta and the inputs
                    A_cf_link = sparse(ones(costfunc_knot_num ...
                        + obj.Quality.Dim + 1, 1), ...
                        [zeta_indices{marg_id}; quality_indices; ...
                        input_indices(marg_id)], ...
                        [costfunc_knots; costfunc_weights; -1], ...
                        1, var_counter);
                    rhs_cf_link = 0;
                    sense_cf_link = '=';

                    % constraint that the sum of all xi coefficients must
                    % be equal to 1
                    A_tf_sum = sparse( ...
                        ones(testfunc_knot_num, 1), ...
                        xi_indices{marg_id}, ...
                        ones(testfunc_knot_num, 1), ...
                        1, var_counter);
                    rhs_tf_sum = 1;
                    sense_tf_sum = '=';

                    % constraints identifying non-zero xi's with eta's
                    A_tf_id_r_cell = cell(testfunc_bit_len(marg_id), 1);
                    A_tf_id_c_cell = cell(testfunc_bit_len(marg_id), 1);
                    A_tf_id_v_cell = cell(testfunc_bit_len(marg_id), 1);

                    for bit_id = 1:testfunc_bit_len(marg_id)
                        tf_bit0_list = ...
                            testfunc_bisect_cell{marg_id}{bit_id, 1};
                        tf_bit0_num = length(tf_bit0_list);
                        tf_bit1_list = ...
                            testfunc_bisect_cell{marg_id}{bit_id, 2};
                        tf_bit1_num = length(tf_bit1_list);

                        A_tf_id_r_cell{bit_id} = (bit_id - 1) * 2 ...
                            + [ones(tf_bit0_num + 1, 1); ...
                            2 * ones(tf_bit1_num + 1, 1)];
                        A_tf_id_c_cell{bit_id} = ...
                            [xi_indices{marg_id}(tf_bit0_list); ...
                            eta_indices{marg_id}(bit_id); ...
                            xi_indices{marg_id}(tf_bit1_list); ...
                            eta_indices{marg_id}(bit_id)];
                        A_tf_id_v_cell{bit_id} = ...
                            [ones(tf_bit0_num, 1); -1; ...
                            ones(tf_bit1_num, 1); 1];
                    end

                    A_tf_id = sparse(vertcat(A_tf_id_r_cell{:}), ...
                        vertcat(A_tf_id_c_cell{:}), ...
                        vertcat(A_tf_id_v_cell{:}), ...
                        2 * testfunc_bit_len(marg_id), var_counter);
                    rhs_tf_id = repmat([0; 1], ...
                        testfunc_bit_len(marg_id), 1);
                    sense_tf_id = repmat('<', ...
                        2 * testfunc_bit_len(marg_id), 1);

                    % constraint linking xi and the inputs
                    A_tf_link = sparse(ones(testfunc_knot_num + 1, 1), ...
                        [xi_indices{marg_id}; input_indices(marg_id)], ...
                        [testfunc_knots; -1], ...
                        1, var_counter);
                    rhs_tf_link = 0;
                    sense_tf_link = '=';
    
                    A_cell{marg_id} = [A_cf_sum; A_cf_id; A_cf_link; ...
                        A_tf_sum; A_tf_id; A_tf_link];
                    rhs_cell{marg_id} = [rhs_cf_sum; rhs_cf_id; ...
                        rhs_cf_link; ...
                        rhs_tf_sum; rhs_tf_id; rhs_tf_link];
                    sense_cell{marg_id} = [sense_cf_sum; sense_cf_id; ...
                        sense_cf_link; ...
                        sense_tf_sum; sense_tf_id; sense_tf_link];
    
                    obj_coef_mat_cell{marg_id} = -speye(testfunc_knot_num);

                    % matrix and intercept vector for directly extracting
                    % the values of the test functions for the marginal
                    % from the results of the solver
                    tfcoef_extract_mat_cell{marg_id} = sparse( ...
                        (1:testfunc_knot_num)', ...
                        xi_indices{marg_id}, 1, ...
                        testfunc_knot_num, var_counter);
                    tfcoef_extract_intercept_cell{marg_id} = ...
                        zeros(testfunc_knot_num, 1);
                end
    
                gl_model.A = [vertcat(A_cell{:}); ...
                    [sparse(size(quality_model.A, 1), marg_num), ...
                    quality_model.A, ...
                    sparse(size(quality_model.A, 1), var_counter ...
                    - marg_num - obj.Quality.Dim ...
                    - obj.Quality.AuxVarNum)]];
                gl_model.rhs = [vertcat(rhs_cell{:}); quality_model.rhs];
                gl_model.sense = [vertcat(sense_cell{:}); ...
                    quality_model.sense];
                gl_model.input_indices = (1:marg_num)';
                gl_model.quality_indices = marg_num + (1:obj.Quality.Dim)';
    
                % store the indices in the objective vector where the
                % coefficients of the test functions will be placed
                gl_model.obj_coef_indices = vertcat(xi_indices{:});
                gl_model.obj_coef_mat = blkdiag(obj_coef_mat_cell{:});

                % matrix and intercept vector for directly extracting
                % the values of the test functions for the marginal
                % from the results of the solver
                gl_model.tfcoef_extract_mat = ...
                    vertcat(tfcoef_extract_mat_cell{:});
                gl_model.tfcoef_extract_intercept = ...
                    vertcat(tfcoef_extract_intercept_cell{:});
    
                obj.Storage.GlobalMinGurobiModel = gl_model;
            else
                error(['unknown formulation of ' ...
                    'the global minimization problem']);
            end

            obj.Storage.SimplicialTestFuncsInitialized = true;

            % updating the simplicial test functions will invalidate all
            % quantities in the runtime environment, thus all variables in
            % the runtime environment need to be flushed
            obj.Runtime = [];
        end

        function [coup_indices, coup_qualities, coup_costs, coup_probs] ...
                = generateHeuristicCoupling(obj)
            % Heuristically couple the marginals by applying comonotone 
            % coupling of the marginals
            % Outputs:
            %   coup_indices: the atoms in the input space in the coupled 
            %   discrete measure represented as indices of knots in the
            %   simplicial test functions
            %   coup_qualities: the atoms in the quality space in the 
            %   coupled discrete measure
            %   coup_costs: the corresponding values of the cost functions
            %   evaluated at the combinations of inputs and qualities
            %   coup_probs: the probabilities of atoms in the coupled
            %   discrete measure

            marg_num = length(obj.Marginals);

            % retrieve some information from the marginals
            atoms_cell = cell(marg_num, 1);
            probs_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                atoms_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Knots;
                probs_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;
            end

            [coup_indices, coup_probs] = comonotone_coupling( ...
                atoms_cell, probs_cell);

            % transform knot indices into actual coordinates
            coup_inputs = zeros(length(coup_probs), marg_num);

            for marg_id = 1:marg_num
                coup_inputs(:, marg_id) = ...
                    atoms_cell{marg_id}(coup_indices(:, marg_id));
            end

            [coup_qualities, coup_costs] = ...
                obj.computeMinCostOverQuality(coup_inputs);
        end
    
        function vals = evaluateOptParametricFunc(obj, ...
                marg_id, pts, batch_size)
            % Evaluate the parametric function on the support of a marginal
            % parametrized by simplicial test functions that is optimized
            % using the cutting plane algorithm. This corresponds to the
            % primal solution of the LSIP problem. The computation is done
            % is batches if necessary.
            % Inputs: 
            %   marg_id: the index of the marginal
            %   pts: vector containing the input points
            %   batch_size: the maximum number of input points to be
            %   handled at the same time in the vectorized procedure
            %   (default is 1e4)
            % Output:
            %   vals: vector containing the computed function values
            
            if ~isfield(obj.Runtime, 'PrimalSolution') ...
                    || isempty(obj.Runtime.PrimalSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            marg_num = length(obj.Marginals);
            marg = obj.Marginals{marg_id};

            % divide the constant term by marg_num and add that to each of
            % the parametric functions
            vals = marg.evaluateWeightedSumOfSimplicialTestFuncs(pts, ...
                obj.Runtime.PrimalSolution.Coefficients{marg_id}, ...
                batch_size) + obj.Runtime.PrimalSolution.Constant ...
                / marg_num;
        end

        function vals_mat = evaluateOptTransferFuncs(obj, pts, ref_pt, ...
                check_inside, batch_size)
            % Evaluate the transfer functions resulted from optimized
            % parametric functions. These transfer functions is part of
            % approximate matching equilibria. The computation is done
            % is batches if necessary.
            % Inputs:
            %   pts: matrix containing the input points
            %   ref_pt: vector indicating a reference point where the 
            %   transfer function will evaluate to 0 (default is the sample
            %   point from the quality space)
            %   check_inside: boolean indicating whether to check if the
            %   points are inside the quality space (default is true)
            %   batch_size: the maximum number of input points to be
            %   handled at the same time in the vectorized procedure
            %   (default is 1e4)
            % Output:
            %   vals_mat: matrix where each column contains the computed
            %   transfer function corresponding to a marginal at the input
            %   points

            if ~isfield(obj.Runtime, 'PrimalSolution') ...
                    || isempty(obj.Runtime.PrimalSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 10001;
            end

            if ~exist('check_inside', 'var') || isempty(check_inside)
                check_inside = true;
            end

            if ~exist('ref_pt', 'var') || isempty(ref_pt)
                ref_pt = obj.Quality.SamplePoint;
            else
                assert(obj.checkIfInsideQualitySpace(ref_pt'), ...
                    'the reference point is not in the quality space');
            end

            marg_num = length(obj.Marginals);

            % add the reference point
            pts = [ref_pt'; pts];
            pt_num = size(pts, 1);

            if check_inside
                assert(all(obj.checkIfInsideQualitySpace(pts)), ...
                    'some points are outside the quality space');
            end

            batch_num = ceil(pt_num / batch_size);
            vals_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                cur_pts = pts(((batch_id - 1) * batch_size ...
                    + 1):min(batch_id * batch_size, pt_num), :);
                cur_batch_size = size(cur_pts, 1);
                
                % only evaluate the first (marg_num - 1) transfer
                % functions, the last one is determined by the balance
                % condition
                vals_batch = zeros(cur_batch_size, marg_num - 1);

                for marg_id = 1:marg_num - 1
                    marg = obj.Marginals{marg_id};
                    marg_knots = marg.SimplicialTestFuncs.Knots;
                    marg_knot_num = length(marg_knots);
                    coef = ...
                        obj.Runtime.PrimalSolution.Coefficients{marg_id};

                    % knots in the cost function for the input in the
                    % support of the marginal while the quality is fixed
                    shifted_cf_knots = obj.CostFuncs{marg_id}.Knots' ...
                        + cur_pts * obj.CostFuncs{marg_id}.Weights;

                    shifted_cf_knots_vec = shifted_cf_knots(:);
                    shifted_cf_knots_inside = ...
                        marg.checkIfInsideSupport(shifted_cf_knots_vec);
                    parfunc_vals = ...
                        -inf(length(shifted_cf_knots_inside), 1);
                    parfunc_vals(shifted_cf_knots_inside) = ...
                        marg.evaluateWeightedSumOfSimplicialTestFuncs( ...
                        shifted_cf_knots_vec(shifted_cf_knots_inside), ...
                        coef);

                    % values attained at the knots of the cost functions
                    vals1 = obj.CostFuncs{marg_id}.Values' - reshape( ...
                        parfunc_vals, cur_batch_size, ...
                        length(obj.CostFuncs{marg_id}.Knots));

                    % values attained at the knots of the test functions
                    vals2 = reshape(obj.evaluateCostFunc(marg_id, ...
                        repelem(marg_knots, cur_batch_size, 1), ...
                        repmat(cur_pts, marg_knot_num, 1)), ...
                        cur_batch_size, marg_knot_num) - coef';

                    vals_batch(:, marg_id) = min([vals1, vals2], [], 2);
                end

                vals_cell{batch_id} = vals_batch;

                if batch_id == 1
                    % store the values of the transfer functions before
                    % shifting evaluated at the reference point
                    ref_vals = vals_cell{batch_id}(1, :);

                    % shift the transfer functions to make the values
                    % vanish at the reference point, then remove the
                    % reference point
                    vals_cell{batch_id} = vals_cell{batch_id}(2:end, :) ...
                        - ref_vals;
                else
                    % shift the transfer functions to make the values
                    % vanish at the reference point
                    vals_cell{batch_id} = vals_cell{batch_id} - ref_vals;
                end

                % the last transfer function is obtained from the balance
                % condition, the resulting transfer functions are
                % guaranteed to add up to 0
                vals_cell{batch_id} = [vals_cell{batch_id}, ...
                    -sum(vals_cell{batch_id}, 2)];
            end

            vals_mat = vertcat(vals_cell{:});
        end

        function LB = getMTLowerBound(obj)
            % Retrieve the computed lower bound for the matching for teams
            % problem
            % Output:
            %   LB: the computed lower bound

            if ~isfield(obj.Runtime, 'PrimalSolution') ...
                    || isempty(obj.Runtime.PrimalSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            LB = -obj.Runtime.LSIP_UB;
        end

        function distri = getOptDiscreteQualityDistr(obj)
            % Retrieve the optimized discrete quality distribution. This is
            % a part of an approximate matching equilibrium.
            % Output:
            %   distri: struct containing fields Probabilities and Atoms

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            distri = struct;
            distri.Probabilities = obj.Runtime.DualSolution.Probabilities;
            distri.Atoms = obj.Runtime.DualSolution.QualityAtoms;
        end

        function [coup_indices, coup_qualities, coup_costs, coup_probs] ...
                = generateDiscreteCoupling(obj, quality_vertices)
            % Generate a feasible dual solution by solving the discretized
            % version of the problem with fixed atoms in the quality space
            % Input:
            %   quality_vertices: fixed atoms in the quality space used for
            %   formulating a discrete version of the problem
            % Outputs:
            %   coup_indices: the atoms in the input space in the coupled 
            %   discrete measure represented as indices of knots in the
            %   simplicial test functions
            %   coup_qualities: the atoms in the quality space in the 
            %   coupled discrete measure
            %   coup_costs: the corresponding values of the cost functions
            %   evaluated at the combinations of inputs and qualities
            %   coup_probs: the probabilities of atoms in the coupled
            %   discrete measure

            marg_num = length(obj.Marginals);
            quality_vert_num = size(quality_vertices, 1);
            marg_knot_num_list = obj.Storage.MargKnotNumList;
            decivar_num = sum(marg_knot_num_list) ...
                + (quality_vert_num - 1) * marg_num;

            % we build the dual version of the problem (maximization) since 
            % it is easier to code
            objective_cell = cell(marg_num, 1);
            A_eq_cell = cell(marg_num, 1);
            A_ineq_cell = cell(marg_num, 1);
            rhs_ineq_cell = cell(marg_num, 1);

            constr_num_list = zeros(marg_num, 1);
            marg_indices_cell = cell(marg_num, 1);
            marg_points_cell = cell(marg_num, 1);
            quality_indices_cell = cell(marg_num, 1);
            quality_points_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_knots = marg.SimplicialTestFuncs.Knots;
                marg_integral = marg.SimplicialTestFuncs.Integrals;
                marg_knot_num = marg_knot_num_list(marg_id);
                
                % the integrals are placed into the objective vector; note
                % that the first test function for the marginal is removed
                % for identification
                objective_cell{marg_id} = [1; ...
                    marg_integral(2:end); ...
                    zeros(quality_vert_num - 1, 1)];

                % matrix for the equality constraints requiring that the
                % transfer functions must sum up to 0; note that the first
                % test function on the quality space is removed for
                % identification
                A_eq_cell{marg_id} = [sparse(quality_vert_num - 1, ...
                    marg_knot_num), speye(quality_vert_num - 1)];

                % matrix for the inequality constraints will have all
                % possible combinations of the vertices in the test
                % functions for the marginal and the vertices in the test
                % functions on the quality space
                [marg_grid, q_grid] = meshgrid(1:marg_knot_num, ...
                    1:quality_vert_num);
                marg_grid_indices = marg_grid(:);
                q_grid_indices = q_grid(:);
                A_marg_full = sparse( ...
                    1:marg_knot_num * quality_vert_num, ...
                    marg_grid_indices, 1, ...
                    marg_knot_num * quality_vert_num, marg_knot_num);
                A_quality_full = sparse( ...
                    1:marg_knot_num * quality_vert_num, ...
                    q_grid_indices, 1, ...
                    marg_knot_num * quality_vert_num, quality_vert_num);

                % remove the first test function for the marginal and the
                % first test function on the quality space; also prepend a
                % column of 1s
                A_ineq_cell{marg_id} = [sparse(ones(marg_knot_num ...
                    * quality_vert_num, 1)), ...
                    A_marg_full(:, 2:end), ...
                    A_quality_full(:, 2:end)];

                % the right-hand side of the inequality constraints are
                % computed from the corresponding coordinates of the
                % vertices
                marg_grid_pts = marg_knots(marg_grid_indices, :);
                quality_grid_pts = quality_vertices(q_grid_indices, :);
                rhs_ineq_cell{marg_id} = ...
                    obj.evaluateCostFunc(marg_id, marg_grid_pts, ...
                    quality_grid_pts);

                constr_num_list(marg_id) = marg_knot_num ...
                    * quality_vert_num;
                marg_indices_cell{marg_id} = marg_grid_indices;
                marg_points_cell{marg_id} = marg_grid_pts;
                quality_indices_cell{marg_id} = q_grid_indices;
                quality_points_cell{marg_id} = quality_grid_pts;
            end
            
            % build a LP model in gurobi
            model = struct;
            model.modelsense = 'max';
            model.objcon = 0;

            % the coefficients corresponding to the first test function of
            % each marginal is not included in the decision variables for
            % identification purposes
            model.obj = vertcat(objective_cell{:});

            model.lb = -inf(decivar_num, 1);
            model.ub = inf(decivar_num, 1);

            % store the equality constraints as fields of the model
            A_eq = horzcat(A_eq_cell{:});
            rhs_eq = zeros(quality_vert_num - 1, 1);

            % assemble the constraints in the model
            model.A = [A_eq; blkdiag(A_ineq_cell{:})];
            model.rhs = [rhs_eq; vertcat(rhs_ineq_cell{:})];
            model.sense = [repmat('=', length(rhs_eq), 1); ...
                repmat('<', length(model.rhs) - length(rhs_eq), 1)];

            % parameters of the LP solver
            LP_options = struct;
            LP_options.OutputFlag = 1;
            LP_options.FeasibilityTol = 1e-9;

            result = gurobi(model, LP_options);

            if ~strcmp(result.status, 'OPTIMAL')
                error('unexpected error in LP');
            end

            probs = result.pi(quality_vert_num:end);
            constr_counter = 0;

            marg_coup_atoms_cell = cell(marg_num, 1);
            marg_coup_probs_cell = cell(marg_num, 1);

            % unwrap the dual optimal solution
            for marg_id = 1:marg_num
                coup_probs = probs(constr_counter ...
                    + (1:constr_num_list(marg_id)));

                % retain only those atoms with positive probabilities
                pos_list = coup_probs > 0;
                coup_probs = coup_probs(pos_list);
                coup_quality_indices = ...
                    quality_indices_cell{marg_id}(pos_list);

                % compute the coupling between the measure on the quality
                % space and the marginal
                marg_coup_atoms_cell{marg_id} = ...
                    [coup_quality_indices, ...
                    marg_indices_cell{marg_id}(pos_list, :)];
                marg_coup_probs_cell{marg_id} = coup_probs;

                if marg_id == 1
                    % use the first coupling as reference (due to the
                    % constraints in the LP problem, all those couplings
                    % will have identical marginals on the quality space,
                    % up to small numerical differences)
                    quality_meas_probs = accumarray( ...
                        coup_quality_indices, coup_probs, ...
                        [quality_vert_num, 1]);

                    % remove all the atoms in the measure on the quality
                    % space with zero probabilities
                    quality_pos_atom_list = find(quality_meas_probs > 0);
                    quality_meas_probs = ...
                        quality_meas_probs(quality_pos_atom_list);
                    quality_meas_atom_num = length(quality_meas_probs);

                    % construct a mapping to identify the indices of atoms
                    % in the quality space for the rest of the couplings
                    quality_pos_atom_mapping = nan(quality_vert_num, 1);
                    quality_pos_atom_mapping(quality_pos_atom_list) ...
                        = 1:quality_meas_atom_num;
                else
                    % retain only those atoms that are in the reference
                    % (additional atoms may occur due to numerical
                    % inaccuracies, but these additional atoms will have
                    % probabilities very close to zero)
                    marg_coup_keep_list = ismember( ...
                        coup_quality_indices, quality_pos_atom_list);
                    marg_coup_probs_cell{marg_id} = ...
                        marg_coup_probs_cell{marg_id}(marg_coup_keep_list);

                    % use the mapping to re-label the atoms using the
                    % indices in the reference measure
                    marg_coup_atoms_cell{marg_id} = ...
                        marg_coup_atoms_cell{marg_id}( ...
                        marg_coup_keep_list, :);
                end

                marg_coup_atoms_cell{marg_id}(:, 1) = ...
                    quality_pos_atom_mapping( ...
                    marg_coup_atoms_cell{marg_id}(:, 1));

                % update the constraint counter
                constr_counter = constr_counter + constr_num_list(marg_id);
            end

            % perform a discrete reassembly (which has the effect of
            % constructing a binding)
            [coup_indices, coup_probs] = discrete_reassembly( ...
                repmat((1:quality_meas_atom_num)', 1, marg_num), ...
                quality_meas_probs, ...
                marg_coup_atoms_cell, marg_coup_probs_cell, 1:marg_num);

            % fill in the values of the actual knots
            coup_atom_num = length(coup_probs);
            coup_inputs = zeros(coup_atom_num, marg_num);

            for marg_id = 1:marg_id
                coup_inputs(:, marg_id) = obj.Marginals{marg_id ...
                    }.SimplicialTestFuncs.Knots(coup_indices(:, marg_id));
            end

            % compute the optimized quality points and cost function values
            [coup_qualities, coup_costs] = ...
                obj.computeMinCostOverQuality(coup_inputs);
        end
        
        function samps = randSampleFromOptJointDistr(obj, ...
                samp_num, rand_stream)
            % Generate independent random sample from the joint
            % distribution consisting of the coupling of the discretized
            % marginals, the corresponding minimizers in the quality space,
            % the coupled continuous marginals, and the corresponding
            % continuous measure in the quality space
            % Inputs: 
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            % Output:
            %   samps: struct containing the following fields:
            %       DiscreteInputs: cell array containing the coordinates
            %       of samples from the coupling of the discretized
            %       marginals
            %       DiscreteQualities: matrix containing the coordinates of
            %       samples from the discrete measure on the quality space
            %       ContinuousInputs: cell array containing the coordinates
            %       of samples from the continuous marginals
            %       ContinuousQualities: matrix containing the coordinates 
            %       of samples from the continuous measure on the quality 
            %       space

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            marg_num = length(obj.Marginals);
            [samps, atom_indices] ...
                = obj.doRandSampleFromPartialReassembly( ...
                (1:marg_num)', samp_num, rand_stream);
            
            % get the corresponding samples from the discrete measure on
            % the quality space by the atom indices
            samps.DiscreteQualities = ...
                obj.Runtime.DualSolution.QualityAtoms(atom_indices, :);

            % the samples from the continuous measure on the quality space 
            % are obtained by solving the minimization of the sum of cost 
            % functions over the quality space
            samps.ContinuousQualities = obj.computeMinCostOverQuality( ...
                horzcat(samps.ContinuousInputs{:}));
        end

        function coup = randSampleFromOptDiscreteCoupling(obj, ...
                marg_id, samp_num, rand_stream)
            % Generate independent random sample from the coupling between
            % the discrete measure on the quality space and a continuous
            % marginal
            % Inputs: 
            %   marg_id: index of the marginal
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            % Output:
            %   coup: struct containing the following fields:
            %       Qualities: matrix containing the coordinates of samples
            %       from the discrete measure on the quality space
            %       AgentTypes: vector containing the coordinates of 
            %       samples from the continuous marginal

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            [samps, atom_indices] ...
                = obj.doRandSampleFromPartialReassembly( ...
                marg_id, samp_num, rand_stream);

            coup = struct;
            coup.AgentTypes = samps.ContinuousInputs{marg_id};

            % get the corresponding samples from the discrete measure on
            % the quality space by the atom indices
            coup.Qualities = ...
                obj.Runtime.DualSolution.QualityAtoms(atom_indices, :);
        end

        function ql_samps = randSampleFromOptContinuousQualityDistr( ...
                obj, samp_num, rand_stream)
            % Generate independent random sample from the continuous
            % measure on the quality space
            % Inputs: 
            %   marg_id: index of the marginal
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            % Output:
            %   ql_samps: matrix containing the coordinates of samples from
            %   the discrete measure on the quality space

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            samps = obj.randSampleFromOptJointDistr(samp_num, ...
                rand_stream);
            ql_samps = samps.ContinuousQualities;
        end

        function coup = randSampleFromOptContinuousCoupling(obj, ...
                marg_id, samp_num, rand_stream)
            % Generate independent random sample from the coupling between
            % the continuous measure on the quality space and a continuous
            % marginal
            % Inputs: 
            %   marg_id: index of the marginal
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            % Output:
            %   coup: struct containing the following fields:
            %       Qualities: matrix containing the coordinates of samples
            %       from the continuous measure on the quality space
            %       AgentTypes: vector containing the coordinates of 
            %       samples from the continuous marginal

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            samps = obj.randSampleFromOptJointDistr(samp_num, ...
                rand_stream);
            coup = struct;
            coup.Qualities = samps.ContinuousQualities;
            coup.AgentTypes = samps.ContinuousInputs{marg_id};
        end

        function [UB_disc, samps] = getMTUpperBoundDiscrete(obj, ...
                samp_num, rand_stream, batch_size)
            % Compute the upper bound for the matching for teams problem
            % based on the discrete measure on the quality space. The 
            % bounds are approximated by Monte Carlo integration.
            % Computation is done in batches if necessary.
            % Inputs: 
            %   samp_num: number of samples for Monte Carlo integration
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: number of samples to be handled simultaneously
            %   during the computation of the sum of the cost functions
            %   (default is 1e4)
            % Output:
            %   UB_disc: upper bound for the matching for teams problem
            %   based on the discrete measure on the quality space
            %   samps: the Monte Carlo samples used in the approximation of
            %   the bounds

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') ...
                    || isempty(batch_size)
                batch_size = 1e4;
            end

            marg_num = length(obj.Marginals);
            [samps, atom_indices] ...
                = obj.doRandSampleFromPartialReassembly( ...
                (1:marg_num)', samp_num, rand_stream);
            
            % get the corresponding samples from the discrete measure on
            % the quality space by the atom indices
            samps.DiscreteQualities = ...
                obj.Runtime.DualSolution.QualityAtoms(atom_indices, :);

            UB_disc_list = obj.evaluateSumOfCostFuncs( ...
                horzcat(samps.ContinuousInputs{:}), ...
                samps.DiscreteQualities, batch_size);

            UB_disc = mean(UB_disc_list);

            samps.UBDiscrete = UB_disc_list;
        end

        function [UB_disc, UB_cont, samps] = getMTUpperBounds(obj, ...
                samp_num, rand_stream, batch_size)
            % Compute two upper bounds for the matching for teams problem
            % based on the discrete measure on the quality space and based
            % on the continuous measure on the quality space. The bounds
            % are approximated by Monte Carlo integration. Computation is 
            % done in batches if necessary.
            % Inputs: 
            %   samp_num: number of samples for Monte Carlo integration
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: number of samples to be handled simultaneously
            %   during the computation of the sum of the cost functions 
            %   (default is 1e4)
            % Output:
            %   UB_disc: upper bound for the matching for teams problem
            %   based on the discrete measure on the quality space
            %   UB_cont: upper bound for the matching for teams problem
            %   based on the continuous measure on the quality space
            %   samps: the Monte Carlo samples used in the approximation of
            %   the bounds

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') ...
                    || isempty(batch_size)
                batch_size = 1e4;
            end

            samps = obj.randSampleFromOptJointDistr(samp_num, ...
                rand_stream);

            UB_disc_list = obj.evaluateSumOfCostFuncs( ...
                horzcat(samps.ContinuousInputs{:}), ...
                samps.DiscreteQualities, batch_size);
            UB_cont_list = obj.evaluateSumOfCostFuncs( ...
                horzcat(samps.ContinuousInputs{:}), ...
                samps.ContinuousQualities, batch_size);

            UB_disc = mean(UB_disc_list);
            UB_cont = mean(UB_cont_list);

            samps.UBDiscrete = UB_disc_list;
            samps.UBContinuous = UB_cont_list;
        end

        function EB = getMTErrorBoundBasedOnOT(obj)
            % Compute the error bound for the objective value of the 
            % matching for teams problem based on the Lipschitz constant of 
            % the cost functions and the optimal transport distances 
            % between the marginals and their discretizations
            % Output:
            %   EB: the computed error bound

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            marg_num = length(obj.Marginals);
            EB = obj.Runtime.LSIP_UB - obj.Runtime.LSIP_LB;

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_atoms = obj.Runtime.DualSolution.InputAtoms(:, ...
                    marg_id);
                marg_probs = obj.Runtime.DualSolution.Probabilities;

                % combine atoms that are very close to each other
                [~, uind, mapping] = unique(round(marg_atoms, 6), ...
                    'stable');
                marg_atoms = marg_atoms(uind);
                marg_probs = accumarray(mapping, marg_probs, ...
                    [length(marg_atoms), 1]);

                marg.setCoupledDiscreteMeasure(marg_atoms, marg_probs);
                EB = EB + obj.CostFuncs{marg_id}.LipschitzConstMarg ...
                    * marg.computeOTCost();
            end
        end

        function EB = getMTTheoreticalErrorBound(obj, LSIP_tolerance)
            % Compute the theoretical error bound for the objective value 
            % of the matching for teams problem via the Lipschitz constant 
            % of the cost functions and the mesh sizes of the simplicial 
            % covers for the marginals
            % Input:
            %   LSIP_tolerance: tolerance value used in the computation of
            %   the LSIP problem
            % Output:
            %   EB: the computed error bound

            marg_num = length(obj.Marginals);
            EB = LSIP_tolerance;

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                EB = EB + obj.CostFuncs{marg_id}.LipschitzConstMarg ...
                    * 2 * marg.SimplicialTestFuncs.MeshSize;
            end
        end
    end

    methods(Access = protected)

        function prepareRuntime(obj)
            % Prepare the runtime environment by initializing some
            % variables
            prepareRuntime@LSIPMinCuttingPlaneAlgo(obj);

            % initialize the cuts to be empty
            obj.Runtime.CutInputs = zeros(0, length(obj.Marginals));
            obj.Runtime.CutQualities = zeros(0, obj.Quality.Dim);
        end

        function model = generateInitialMinModel(obj)
            % Generate the initial linear programming model for gurobi
            % Output:
            %   model: struct containing the linear programming model in
            %   gurobi

            model = struct;
            marg_num = length(obj.Marginals);
            decivar_num = length(obj.Storage.DeciVarIndicesInTestFuncs) ...
                + 1;

            integrals_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                integrals_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;
            end

            integrals_vec = vertcat(integrals_cell{:});
            
            % since the cutting plane algorithm assumes that the problem is
            % a minimization problem, we need to transform our maximization
            % problem into a minimization problem
            model.modelsense = 'min';
            model.objcon = 0;

            % the coefficients corresponding to the first test function of
            % each marginal is not included in the decision variables for
            % identification purposes
            model.obj = [-1; ...
                -integrals_vec(obj.Storage.DeciVarIndicesInTestFuncs)];

            model.sense = '>';
            model.lb = -inf(decivar_num, 1);
            model.ub = inf(decivar_num, 1);
            model.A = sparse(0, decivar_num);
            model.rhs = zeros(0, 1);
        end

        function [min_lb, optimizers] = callGlobalMinOracle(obj, vec)
            % Given a decision vector, call the global minimization oracle
            % to approximately determine the "most violated" constraints
            % and return a lower bound for the minimal value
            % Input:
            %   vec: vector corresponding to the current LP solution
            % Outputs:
            %   min_lb: lower bound for the global minimization problem
            %   optimizers: struct containing approximate optimizers of the
            %   global minimization problem in the form of inputs in the 
            %   support of the marginals as well as the the corresponding 
            %   points in the quality space

            % disassemble the vector resulted from solving the LP problem

            % the first component corresponds to the constant term; the
            % second component onwards will be filled into the
            % corresponding entries of the vector storing the coefficients
            % corresponding to the test functions
            objective_const = vec(1);
            knots_vals = zeros(obj.Storage.TotalKnotNum, 1);
            knots_vals(obj.Storage.DeciVarIndicesInTestFuncs) = vec(2:end);

            % call the mixed-integer solver
            [pool_inputs, pool_qualities, pool_marg_testfuncs, ...
                pool_vals, pool_costs, LB] ...
                = obj.computeGlobalMinMILP(knots_vals, objective_const);

            min_lb = LB;

            % ignore those approximate minimizers whose objective values
            % are non-negative since they do not generate cuts
            pool_neg_list = pool_vals < 0;

            optimizers = struct;
            optimizers.inputs = pool_inputs(pool_neg_list, :);
            optimizers.qualities = pool_qualities(pool_neg_list, :);
            optimizers.marg_testfunc_vals = ...
                pool_marg_testfuncs(pool_neg_list, :);
            optimizers.costs = pool_costs(pool_neg_list, :);
        end

        function [pool_inputs, pool_qualities, pool_marg_testfuncs, ...
                pool_vals, pool_costs, LB] ...
                = computeGlobalMinMILP(obj, knot_vals, objective_const)
            % Approximately solve the global minimization problem using the
            % gurobi mixed-integer solver. The algorithm will retain up to 
            % obj.GlobalOptions.pool_size approximate optimizers.
            % Input:
            %   vert_vals: values of the simplicial test functions at the
            %   knots concatenated into a vector
            %   objective_const: constant part of the objective function
            % Outputs:
            %   pool_inputs: matrix containing the inputs in the support of 
            %   the marginals from the approximate optimizers
            %   pool_qualities: matrix containing the quality part of the
            %   approximate optimizers
            %   pool_marg_testfuncs: sparse matrix containing the values of
            %   the simplicial test functions for the marginals evaluated 
            %   at the resulting input points
            %   pool_vals: the corresponding objective values of the
            %   approximate optimizers
            %   pool_costs: the corresponding values of the sum of the cost
            %   functions evaluated at the approximate optimizers on the
            %   quality space
            %   LB: the best lower bound found throughout the algorithm

            model = obj.Storage.GlobalMinGurobiModel;

            % this vector has the coefficients of all test functions set to
            % zero; thus it can be used to evaluate the sum of the cost
            % functions at the approximate optimizers on the quality space
            costfunc_obj = model.obj;
            costfunc_obj_const = model.objcon;

            model.obj(model.obj_coef_indices) = model.obj_coef_mat ...
                * knot_vals;
            model.objcon = model.objcon - objective_const;

            % options of the mixed-integer solver
            gl_options = struct;
            gl_options.OutputFlag = 0;
            gl_options.IntFeasTol = 1e-6;
            gl_options.FeasibilityTol = 1e-8;
            gl_options.OptimalityTol = 1e-8;
            gl_options.PoolSolutions = 100;
            gl_options.PoolGap = 0.8;
            gl_options.NodefileStart = 2;
            gl_options.BestBdStop = 1e-6;
            gl_options.BestObjStop = -inf;
            gl_options.MIPGap = 1e-4;
            gl_options.MIPGapAbs = 1e-10;

            % set the additional options for the mixed-integer solver
            gl_options_fields = fieldnames(obj.GlobalOptions);
            gl_options_values = struct2cell(obj.GlobalOptions);

            for fid = 1:length(gl_options_fields)
                gl_options.(gl_options_fields{fid}) ...
                    = gl_options_values{fid};
            end

            result = gurobi(model, gl_options);

            if ~strcmp(result.status, 'OPTIMAL') ...
                    && ~strcmp(result.status, 'USER_OBJ_LIMIT')
                error('error in the mixed-integer solver');
            end

            LB = result.objbound;

            % get a set of approximate optimizers
            if isfield(result, 'pool')
                pool_cuts = horzcat(result.pool.xn)';
                pool_vals = [result.pool.objval]';
            else
                pool_cuts = result.x';
                pool_vals = result.objval;
            end

            pool_inputs = pool_cuts(:, model.input_indices);
            pool_qualities = pool_cuts(:, model.quality_indices);
            pool_marg_testfuncs = pool_cuts * model.tfcoef_extract_mat' ...
                + model.tfcoef_extract_intercept';

            % make sure that the inputs are within their respective bounds
            % (they may be slightly out of bound due to small numerical
            % inaccuracies)
            pool_inputs = min(max(pool_inputs, ...
                obj.Storage.MarginalLowerBounds'), ...
                obj.Storage.MarginalUpperBounds');

            pool_costs = pool_cuts * costfunc_obj + costfunc_obj_const;
        end

        function addConstraints(obj, optimizers)
            % Given a collection of approximate optimizers from the global
            % minimization oracle, generate and add the corresponding
            % linear constraints
            % Inputs:
            %   optimizers: output of the method callGlobalOracle

            marg_num = length(obj.Marginals);

            % generate the constraint matrix by evaluating the simplicial
            % test functions with the inputs
            A_cell = cell(marg_num, 1);

            % there are two possibilities: either the inputs are provided 
            % as coordinates (when the field inputs is set), or the inputs 
            % are provided as knot indices (when the field knot_indices is 
            % set)
            if isfield(optimizers, 'marg_knot_indices')
                constr_num = size(optimizers.marg_knot_indices, 1);
                optimizers.inputs = zeros(constr_num, marg_num);
            else
                constr_num = size(optimizers.inputs, 1);
            end

            if isfield(optimizers, 'marg_testfunc_vals')
                A_full = optimizers.marg_testfunc_vals;
            else
                for marg_id = 1:marg_num
                    marg = obj.Marginals{marg_id};
    
                    if isfield(optimizers, 'marg_knot_indices')
                        A_cell{marg_id} = sparse((1:constr_num)', ...
                            optimizers.marg_knot_indices(:, ...
                            marg_id), 1, ...
                            constr_num, ...
                            obj.Storage.MargKnotNumList(marg_id));
                        optimizers.inputs(:, marg_id) = ...
                            marg.SimplicialTestFuncs.Knots( ...
                            optimizers.marg_knot_indices(:, marg_id));
                    else
                        A_cell{marg_id} = ...
                            marg.evaluateSimplicialTestFuncs( ...
                            optimizers.inputs(:, marg_id));
                    end
                end

                % first generate a matrix containing all test functions 
                % (each row corresponds to an approximate optimizer, each 
                % column corresponds to a test function)
                A_full = horzcat(A_cell{:});
            end

            col_counter = 0;

            for marg_id = 1:marg_num
                col_list = col_counter ...
                    + (1:obj.Storage.MargKnotNumList(marg_id));
                marg_testfunc_vals = A_full(:, col_list);

                % sanitize the output to avoid potential numerical issues;
                % first, check for entries that are very close to 1, set
                % them to 1 and set the remaining entries in those rows to 
                % 0; then, set entries that are close to 0 to 0
                [m_r, m_c, m_v] = find(A_full(:, col_list));
                closeto1_list = abs(m_v - 1) ...
                    < obj.Options.sanitation_threshold;
                m_v(closeto1_list) = 1;
                keep_list = ~(ismember(m_r, m_r(closeto1_list)) ...
                    & ~closeto1_list) & m_v ...
                    >= obj.Options.sanitation_threshold;
                A_full(:, col_list) = sparse(m_r(keep_list), ...
                    m_c(keep_list), m_v(keep_list), ...
                    size(marg_testfunc_vals, 1), ...
                    size(marg_testfunc_vals, 2));

                col_counter = col_counter ...
                    + obj.Storage.MargKnotNumList(marg_id);
            end

            % filter out those test functions whose coefficients are not
            % included in the decision variable, then prepend a column of
            % 1
            A_new = [sparse(ones(constr_num, 1)), ...
                A_full(:, obj.Storage.DeciVarIndicesInTestFuncs)];

            rhs_new = optimizers.costs;

            % add the newly generated constraints to the end
            obj.Runtime.CurrentLPModel.A = ...
                [obj.Runtime.CurrentLPModel.A; -A_new];
            obj.Runtime.CurrentLPModel.rhs = ...
                [obj.Runtime.CurrentLPModel.rhs; -rhs_new];

            if ~isempty(obj.Runtime.vbasis) && ~isempty(obj.Runtime.cbasis)
                obj.Runtime.CurrentLPModel.vbasis = obj.Runtime.vbasis;
                obj.Runtime.CurrentLPModel.cbasis = ...
                    [obj.Runtime.cbasis; zeros(constr_num, 1)];
            else
                if isfield(obj.Runtime.CurrentLPModel, 'vbasis')
                    obj.Runtime.CurrentLPModel ...
                        = rmfield(obj.Runtime.CurrentLPModel, 'vbasis');
                end

                if isfield(obj.Runtime.CurrentLPModel, 'cbasis')
                    obj.Runtime.CurrentLPModel ...
                        = rmfield(obj.Runtime.CurrentLPModel, 'cbasis');
                end
            end

            % add the inputs and qualities to the runtime environment
            obj.Runtime.CutInputs = [obj.Runtime.CutInputs; ...
                optimizers.inputs];
            obj.Runtime.CutQualities = [obj.Runtime.CutQualities; ...
                optimizers.qualities];

            % if obj.Runtime.NumOfInitialConstraints is not set, it means
            % that this is the first call to obj.addConstraints which
            % generates the initial constraints; this number is stored in
            % the runtime environment
            if ~isfield(obj.Runtime, 'NumOfInitialConstraints') ...
                    || isempty(obj.Runtime.NumOfInitialConstraints)
                obj.Runtime.NumOfInitialConstraints = constr_num;
            end
        end

        function reduceConstraints(obj, result)
            % Remove some of the constraints to speed up the LP solver
            % Input:
            %   result: the output from the gurobi LP solver

            if ~isinf(obj.Options.reduce.thres) ...
                    && obj.Runtime.iter > 0 ...
                    && obj.Runtime.iter <= obj.Options.reduce.max_iter ...
                    && mod(obj.Runtime.iter, obj.Options.reduce.freq) == 0

                % the list of constraints to be kept (here since the
                % directions of the inequalities are all >=, the slackness
                % is non-positive; the threshold specifies the maximum
                % absolute value of slackness)
                keep_list = result.slack >= -obj.Options.reduce.thres;

                % always keep the initial constraints
                keep_list(1:obj.Runtime.NumOfInitialConstraints) = true;

                % update all variables
                obj.Runtime.CutInputs ...
                    = obj.Runtime.CutInputs(keep_list, :);
                obj.Runtime.CutQualities ...
                    = obj.Runtime.CutQualities(keep_list, :);
                obj.Runtime.CurrentLPModel.A ...
                    = obj.Runtime.CurrentLPModel.A(keep_list, :);
                obj.Runtime.CurrentLPModel.rhs ...
                    = obj.Runtime.CurrentLPModel.rhs(keep_list);

                if ~isempty(obj.Runtime.cbasis)
                    obj.Runtime.cbasis = obj.Runtime.cbasis(keep_list);
                end
            end
        end

        function primal_sol = buildPrimalSolution(obj, result, violation)
            % Given the output from gurobi and a lower bound for the
            % optimal value of the global minimization oracle, build the
            % corresponding primal solution
            % Inputs:
            %   result: output of the gurobi LP solver
            %   violation: a lower bound for the global minimization
            %   problem
            % Output:
            %   primal_sol: the constructed potential functions on the
            %   support of the input measures

            vec = result.x;
            primal_sol = struct;
            primal_sol.Constant = vec(1) + violation;

            vert_coefs = zeros(obj.Storage.TotalKnotNum, 1);
            vert_coefs(obj.Storage.DeciVarIndicesInTestFuncs) = vec(2:end);
            marg_num = length(obj.Marginals);
            primal_sol.Coefficients = cell(marg_num, 1);

            for marg_id = 1:marg_num
                primal_sol.Coefficients{marg_id} = vert_coefs( ...
                    obj.Storage.MargKnotNumOffsets(marg_id) ...
                    + (1:obj.Storage.MargKnotNumList(marg_id)));
            end
        end

        function dual_sol = buildDualSolution(obj, result)
            % Given the output from gurobi, build the corresponding dual
            % solution
            % Input:
            %   result: output of the gurobi LP solver
            % Output:
            %   dual_sol: the constructed discrete probability measure for
            %   the relaxed MMOT problem coupled with the corresponding
            %   discrete probability measure on the quality space

            dual_sol = struct;
            pos_list = result.pi > 0;
            dual_sol.Probabilities = result.pi(pos_list);
            dual_sol.InputAtoms = obj.Runtime.CutInputs(pos_list, :);
            dual_sol.QualityAtoms = obj.Runtime.CutQualities(pos_list, :);
        end

        function [samps, disc_atom_index_samps] ...
                = doRandSampleFromPartialReassembly(obj, ...
                marg_to_reassemble, samp_num, rand_stream)
            % Generate independent random samples from a partial reassembly
            % where only couplings of certain marginals are sampled
            % Inputs: 
            %   marg_to_reassemble: vector containing the indices of the
            %   marginals whose corresponding couplings will be sampled
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling
            % Outputs:
            %   samps: struct with the following fields:
            %       DiscreteInputs: cell array containing the coordinates
            %       of samples from the coupling of the discretized
            %       marginals
            %       ContinuousInputs: cell array containing the coordinates
            %       of samples from the continuous marginals (only those
            %       included in marg_to_reassmeble will be sampled)
            %   disc_atom_index_samps: vector containing the indices of
            %   atoms in the discrete coupling

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            dual_sol = obj.Runtime.DualSolution;
            marg_num = length(obj.Marginals);

            if ~isfield(obj.Runtime, 'CouplingDone') ...
                || ~obj.Runtime.CouplingDone
                % set the discrete measures that the marginals are coupled 
                % with
                for marg_id = 1:marg_num
                    obj.Marginals{marg_id}.setCoupledDiscreteMeasure( ...
                        dual_sol.InputAtoms(:, marg_id), ...
                        dual_sol.Probabilities);
                end

                obj.Runtime.CouplingDone = true;
            end

            n = length(obj.Runtime.DualSolution.Probabilities);
            
            % generate random indices of the atoms according to the
            % probabilities
            disc_atom_index_samps = randsample(rand_stream, n, ...
                samp_num, true, obj.Runtime.DualSolution.Probabilities);

            samps = struct;
            samps.DiscreteInputs = cell(marg_num, 1);
            samps.ContinuousInputs = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};

                % store the coordinates of the samples with discrete
                % marginals
                samps.DiscreteInputs{marg_id} ...
                    = obj.Runtime.DualSolution.InputAtoms( ...
                    disc_atom_index_samps, marg_id);

                if ~ismember(marg_id, marg_to_reassemble)
                    % skip those continuous marginals that are not required
                    % to be sampled
                    continue;
                end

                samps.ContinuousInputs{marg_id} = zeros(samp_num, 1);

                % count the number of samples coupled with each of the
                % atoms in the discretized marginal
                atom_num_list = accumarray(disc_atom_index_samps, 1, ...
                    [n, 1]);

                % generate from the conditional distributions
                cont_samp_cell = marg.conditionalRandSample( ...
                    atom_num_list, rand_stream);

                % fill in the coupled samples from the continuous marginals
                for atom_id = 1:length(atom_num_list)
                    samps.ContinuousInputs{marg_id}( ...
                        disc_atom_index_samps == atom_id, :) = ...
                        cont_samp_cell{atom_id};
                end
            end
        end
    end
end

