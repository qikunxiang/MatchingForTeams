classdef OTDiscrete < LSIPMinCuttingPlaneAlgo
    % Class for computing discrete optimal transport with large numbers of atoms using the cutting plane algorithm (effectively a 
    % constraint generate scheme)
    
    properties(GetAccess = public, SetAccess = protected)
        % probabilities in the first discrete probability measure
        Probs1;

        % probabilities in the second discrete probability measure
        Probs2;

        % cost matrix
        Costs;
    end
    
    methods(Access = public)
        function obj = OTDiscrete(probs1, probs2, costs, varargin)
            % Constructor
            % Inputs:
            %   probs1: vector containing probabilities in the first discrete probability measure
            %   probs2: vector containing probabilities in the second discrete probability measure
            %   costs: matrix containing the costs

            obj@LSIPMinCuttingPlaneAlgo(varargin{:});

            % set the default options for reducing constraints
            if ~isfield(obj.Options, 'reduce') || isempty(obj.Options.reduce)
                obj.Options.reduce = struct;
            end

            if ~isfield(obj.Options.reduce, 'thres') || isempty(obj.Options.reduce.thres)
                obj.Options.reduce.thres = inf;
            end

            if ~isfield(obj.Options.reduce, 'freq') || isempty(obj.Options.reduce.freq)
                obj.Options.reduce.freq = 20;
            end

            if ~isfield(obj.Options.reduce, 'max_iter') || isempty(obj.Options.reduce.max_iter)
                obj.Options.reduce.max_iter = inf;
            end

            if ~isfield(obj.GlobalOptions, 'pool_size') || isempty(obj.GlobalOptions.pool_size)
                obj.GlobalOptions.pool_size = 100;
            end

            atom_num1 = length(probs1);
            atom_num2 = length(probs2);

            assert(size(costs, 1) == atom_num1 && size(costs, 2) == atom_num2, 'cost matrix mis-specified');
            assert(abs(sum(probs1) - 1) < 1e-10, 'the first probability measure is not normalized');
            assert(abs(sum(probs2) - 1) < 1e-10, 'the second probability measure is not normalized');

            obj.Probs1 = probs1 / sum(probs1);
            obj.Probs2 = probs2 / sum(probs2);
            obj.Costs = costs;
        end

        function [coup_indices, coup_probs] = generateHeuristicCoupling(obj)
            % Compute a heuristic coupling of the marginals via a greedy algorithm
            % Output:
            %   coup_indices: two-column matrix containing the computed coupled indices
            %   coup_probs: vector containing the computed coupled probabilities

            probs1 = obj.Probs1;
            probs2 = obj.Probs2;
            atom_num1 = length(probs1);
            atom_num2 = length(probs2);

            coup_indices = zeros(atom_num1 + atom_num2, 2);
            coup_probs = zeros(atom_num1 + atom_num2, 1);
            coup_counter = 0;

            % sort the atoms in the first probability measure into an order with descending probabilities to guarantee that atoms with
            % high probabilities are coupled first
            probs1_atom_indices = (1:atom_num1)';
            [probs1_res, sorted_order] = sort(probs1, 'descend');
            probs1_atom_indices = probs1_atom_indices(sorted_order);
            probs2_atom_indices = (1:atom_num2)';
            costs_sorted = obj.Costs(sorted_order, :);
            probs2_available = true(atom_num2, 1);
            probs2_res = probs2;
            
            for i = 1:atom_num1
                costs_list = costs_sorted(i, probs2_available)';

                if isempty(costs_list)
                    % if all the atoms in the second probability measure has been coupled, terminate the loop
                    break;
                end

                probs2_atom_list = probs2_atom_indices(probs2_available);

                % sort the available atoms based on costs
                [costs_list, sorted_order] = sort(costs_list, 'ascend');
                probs2_atom_list = probs2_atom_list(sorted_order);

                for j = 1:length(costs_list)
                    atom_id = probs2_atom_list(j);

                    % couple the two atoms
                    coup_prob = min(probs2_res(atom_id), probs1_res(i));
                    coup_counter = coup_counter + 1;
                    coup_probs(coup_counter) = coup_prob;
                    coup_indices(coup_counter, :) = [probs1_atom_indices(i), probs2_atom_list(j)];

                    probs2_res(atom_id) = probs2_res(atom_id) - coup_prob;
                    probs1_res(i) = probs1_res(i) - coup_prob;

                    if probs2_res(atom_id) >= 1e-14
                        % if the current atom in the first probability measure has been fully coupled, skip to the next one
                        break;
                    else
                        probs2_available(atom_id) = false;
                    end

                    if probs1_res(i) < 1e-14
                        % if the current atom in the first probability measure has been fully coupled, skip to the next one
                        break;
                    end
                end
            end

            coup_indices = coup_indices(1:coup_counter, :);
            coup_probs = coup_probs(1:coup_counter);
        end
    end

    methods(Access = protected)
        function initializeBeforeRun(obj) %#ok<MANU>
            % Initialize the algorithm by preparing some constant quantities

            % do nothing
        end

        function prepareRuntime(obj)
            % Prepare the runtime environment by initializing some variables
            prepareRuntime@LSIPMinCuttingPlaneAlgo(obj);

            % initialize the cuts to be empty
            obj.Runtime.CoupIndices = zeros(0, 2);
        end

        function model = generateInitialMinModel(obj)
            % Generate the initial linear programming model for gurobi
            % Input:
            %   model: gurobi model representing the initial LP problem
            
            model = struct;
            model.modelsense = 'min';
            model.objcon = 0;
            model.obj = -[obj.Probs1; obj.Probs2];
            decivar_num = length(obj.Probs1) + length(obj.Probs2);
            model.A = sparse(0, decivar_num);
            model.rhs = zeros(0, 1);
            model.sense = '>';
            model.lb = -inf(decivar_num, 1);
            model.ub = inf(decivar_num, 1);
        end
        
        function [min_lb, optimizers] = callGlobalMinOracle(obj, vec)
            % Given a decision vector, call the global minimization oracle to approximately determine the "most violated" constraints
            % and return a lower bound for the minimal value
            % Input:
            %   vec: vector representing a sub-optimal solution
            % Outputs: 
            %   min_lb: a lower bound for the global minimization problem
            %   optimizers: two-column matrix containing the atom indices of a collection of approximate minimizers

            atom_num1 = length(obj.Probs1);
            pot1 = vec(1:atom_num1);
            pot2 = vec((atom_num1 + 1):end);
            violation_mat = obj.Costs - pot1 - pot2';

            [vio_r, vio_c, vio_v] = find(min(violation_mat, 0));

            if isempty(vio_v)
                % if all violations are non-negative, just return a minimizer with minimum value 0
                [min_over_col, min_col_indices] = min(violation_mat, [], 2);
                [min_lb, min_row_index] = min(min_over_col);
                optimizers = [min_row_index, min_col_indices(min_row_index)];
            else
                [~, sorted_order] = sort(vio_v, 'ascend');
                sorted_order = sorted_order(1:min(length(sorted_order), obj.GlobalOptions.pool_size));
                min_lb = vio_v(sorted_order(1));
                optimizers = [vio_r(sorted_order), vio_c(sorted_order)];
            end
        end

        function addConstraints(obj, optimizers)
            % Given a collection of approximate optimizers from the global minimization oracle, generate and add the corresponding
            % linear constraints
            % Input:
            %   optimizers: two-column matrix containing the atom indices of a collection of approximate minimizers of the global
            %   minimization problem
            
            atom_num1 = length(obj.Probs1);
            atom_num2 = length(obj.Probs2);

            constr_num = size(optimizers, 1);
            A_new = sparse([(1:constr_num)'; (1:constr_num)'], [optimizers(:, 1); atom_num1 + optimizers(:, 2)], 1, ...
                constr_num, atom_num1 + atom_num2);
            rhs_new = obj.Costs(sub2ind([atom_num1, atom_num2], optimizers(:, 1), optimizers(:, 2)));

            % add the newly generated constraints to the end
            obj.Runtime.CurrentLPModel.A = [obj.Runtime.CurrentLPModel.A; -A_new];
            obj.Runtime.CurrentLPModel.rhs = [obj.Runtime.CurrentLPModel.rhs; -rhs_new];

            if ~isempty(obj.Runtime.vbasis) && ~isempty(obj.Runtime.cbasis)
                obj.Runtime.CurrentLPModel.vbasis = obj.Runtime.vbasis;
                obj.Runtime.CurrentLPModel.cbasis = [obj.Runtime.cbasis; zeros(constr_num, 1)];
            else
                if isfield(obj.Runtime.CurrentLPModel, 'vbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'vbasis');
                end

                if isfield(obj.Runtime.CurrentLPModel, 'cbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'cbasis');
                end
            end

            % add the added indices to the runtime environment
            obj.Runtime.CoupIndices = [obj.Runtime.CoupIndices; optimizers];

            % if obj.Runtime.NumOfInitialConstraints is not set, it means that this is the first call to obj.addConstraints which
            % generates the initial constraints; this number is stored in the runtime environment
            if ~isfield(obj.Runtime, 'NumOfInitialConstraints') || isempty(obj.Runtime.NumOfInitialConstraints)
                obj.Runtime.NumOfInitialConstraints = constr_num;
            end
        end

        function reduceConstraints(obj, result)
            % Remove some of the constraints to speed up the LP solver
            % Input:
            %   result: the output from the gurobi LP solver

            if ~isinf(obj.Options.reduce.thres) && obj.Runtime.iter > 0 && obj.Runtime.iter <= obj.Options.reduce.max_iter ...
                    && mod(obj.Runtime.iter, obj.Options.reduce.freq) == 0

                % the list of constraints to be kept (here since the directions of the inequalities are all >=, the slackness is 
                % non-positive; the threshold specifies the maximum absolute value of slackness)
                keep_list = result.slack >= -obj.Options.reduce.thres;

                % always keep the initial constraints
                keep_list(1:obj.Runtime.NumOfInitialConstraints) = true;

                % update all variables
                obj.Runtime.CoupIndices = obj.Runtime.CoupIndices(keep_list, :);
                obj.Runtime.CurrentLPModel.A = obj.Runtime.CurrentLPModel.A(keep_list, :);
                obj.Runtime.CurrentLPModel.rhs = obj.Runtime.CurrentLPModel.rhs(keep_list);

                if ~isempty(obj.Runtime.cbasis)
                    obj.Runtime.cbasis = obj.Runtime.cbasis(keep_list);
                end
            end
        end

        function updateLSIPUB(obj, min_lb, optimizers) %#ok<INUSD> 
            % Update the LSIP upper bound after each call to the global minimization oracle
            % Inputs:
            %   min_lb: the lower bound for the global minimization problem
            %   optimizers: a set of approximate optimizers of the global minimization problem

            obj.Runtime.LSIP_UB = min(obj.Runtime.LSIP_UB, obj.Runtime.LSIP_LB - min_lb);
        end

        function primal_sol = buildPrimalSolution(obj, result, violation)
            % Given the output from gurobi and a lower bound for the optimal value of the global minimization oracle, build the
            % corresponding primal solution
            % Inputs:
            %   result: output of the gurobi LP solver
            %   violation: a lower bound for the global minimization problem
            % Output:
            %   primal_sol: the constructed potentials at each atom

            atom_num1 = length(obj.Probs1);
            primal_sol = struct;
            primal_sol.Potential1 = result.x(1:atom_num1) + violation;
            primal_sol.Potential2 = result.x((atom_num1 + 1):end);
        end

        function dual_sol = buildDualSolution(obj, result)
            % Given the output from gurobi, build the corresponding dual solution
            % Input:
            %   result: output of the gurobi LP solver
            % Output:
            %   dual_sol: the constructed coupling represented by coupled atom indices

            dual_sol = struct;
            pos_list = result.pi > 0;
            dual_sol.Probabilities = result.pi(pos_list);
            dual_sol.CoupIndices = obj.Runtime.CoupIndices(pos_list, :);
        end
    end
end