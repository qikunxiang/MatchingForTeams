classdef (Abstract) MT2DQuad < LSIPMinCuttingPlaneAlgo
    % Abstract class for matching for teams problems with two-dimensional 
    % marginals, two-dimensional quality space, and quadratic cost 
    % functions. 

    properties(Constant)
        % numerical tolerance for deciding whether a point is inside a
        % polytope
        INSIDE_TOLERANCE = 1e-12;
    end

    properties(GetAccess = public, SetAccess = protected)
        % cell array containing the marginals
        Marginals;

        % vector containing the weights
        MarginalWeights;

        % struct containing information about the quality space
        Quality;

        % constant part of the quadratic cost functions that does not
        % affect the matching
        QuadraticConstant = 0;
    end

    methods(Access = public)
        function obj = MT2DQuad(marginals, weights, quality_cell, ...
                varargin)
            % Constructor method
            % Inputs:
            %   marginals: cell array containing objects of
            %   ProbMeas2D_ConvexPolytope
            %   weights: vector containing the weights corresponding to the
            %   marginals, the weights will sum up to 1
            %   quality_cell: cell array where each cell contains a
            %   two-column matrix indicating the vertices of a polytope

            obj@LSIPMinCuttingPlaneAlgo(varargin{:});

            % set the default options for semi-discrete optimal transport
            if ~isfield(obj.Options, 'OT') ...
                    || isempty(obj.Options.OT)
                obj.Options.OT = struct;
            end

            if ~isfield(obj.Options.OT, 'angle_num') ...
                    || isempty(obj.Options.OT.angle_num)
                obj.Options.OT.angle_num = [];
            elseif isscalar(obj.Options.OT.angle_num)
                obj.Options.OT.angle_num = ones(length(weights), 1) ...
                    * obj.Options.OT.angle_num;
            end

            if ~isfield(obj.Options.OT, 'pp_angle_indices') ...
                    || isempty(obj.Options.OT.pp_angle_indices)
                obj.Options.OT.pp_angle_indices = repmat({[]}, ...
                    length(weights), 1);
            elseif ~iscell(obj.Options.OT.pp_angle_indices)
                obj.Options.OT.pp_angle_indices = ...
                    repmat({obj.Options.OT.pp_angle_indices}, ...
                    length(weights), 1);
            end

            if ~isfield(obj.Options.OT, 'optimization_options')
                obj.Options.OT.optimization_options ...
                    = cell(length(weights), 1);
            elseif ~iscell(obj.Options.OT.optimization_options)
                optim_options = obj.Options.OT.optimization_options;
                obj.Options.OT.optimization_options ...
                    = cell(length(weights), 1);
                obj.Options.OT.optimization_options(:, :) ...
                    = {optim_options};
            end

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

            marg_num = length(weights);
            assert(length(marginals) == marg_num, ...
                'input mis-specified');
            assert(abs(sum(weights) - 1) < 1e-12, ...
                'weights do not sum up to 1');
            
            % normalize the weights to remove numerical inaccuracies
            weights = weights / sum(weights);

            obj.Marginals = marginals;
            obj.MarginalWeights = weights;

            % compute the constant terms in the cost functions that are
            % related to the quadratic expectation with respect to the
            % marginals
            quad_consts = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                quad_consts(marg_id) = sum(diag( ...
                    obj.Marginals{marg_id}.SecondMomentMat));
            end

            obj.QuadraticConstant = quad_consts' * obj.MarginalWeights;

            obj.Quality = struct;
            poly_num = length(quality_cell);
            obj.Quality.Polytopes = cell(poly_num, 1);

            % number of vertices in each polytope
            vert_num_list = zeros(poly_num, 1);

            % store all vertices of the polytopes which will be used in the
            % first iteration of the global minimization algorithm
            vert_cell = cell(poly_num, 1);

            % store one vertex from each polytope in order to compute the
            % rectangles that overlap with the polytopes; this is to
            % resolve an edge case where the entire polytope is contained
            % in a rectangle
            first_vert = zeros(poly_num, 2);

            for poly_id = 1:poly_num
                vertices = quality_cell{poly_id};

                % compute the convex hull of the vertices
                ccorder = convhull(vertices(:, 1), vertices(:, 2), ...
                    'Simplify', true);

                vert_num_list(poly_id) = length(ccorder) - 1;

                % rearrange the vertices in counterclockwise order
                vertices_cc = vertices(ccorder(1:end - 1), :);

                obj.Quality.Polytopes{poly_id} = vertices_cc;
                vert_cell{poly_id} = vertices_cc;
                first_vert(poly_id, :) = vertices_cc(1, :);
            end

            % store all vertices of the polytopes
            obj.Quality.Vertices = unique(vertcat(vert_cell{:}), 'rows');

            % store the maximum norm of points in the quality space
            obj.Quality.MaxNorm = ...
                sqrt(max(sum(obj.Quality.Vertices .^ 2, 2)));

            % compute all hyperplanes characterizing the polytopes in the
            % form of inequalities w' * z <= b
            hp_w_cell = cell(poly_num, 1);
            hp_b_cell = cell(poly_num, 1);

            % hyperplanes with redundant inequality as padding such that
            % each polytope has the same number of inequalities
            hp_w_cell_rd = cell(poly_num, 1);
            hp_b_cell_rd = cell(poly_num, 1);

            % compute the maximum number of vertices among the polytopes;
            % this will be used to add redundant inequalities such that all
            % polytopes are characterized by the same number of
            % inequalities
            vert_num_max = max(vert_num_list);

            for poly_id = 1:poly_num
                % retrieve the vertices in the polytopes and their next
                % counterclockwise neighbor
                vertices_cc1 = obj.Quality.Polytopes{poly_id};
                vertices_cc2 = vertices_cc1([2:end, 1], :);
                vert_num = vert_num_list(poly_id);

                % compute the weights and the intercepts of the hyperplanes
                % characterizing the polytope
                hp_w_cell{poly_id} = ...
                    [vertices_cc2(:, 2) - vertices_cc1(:, 2), ...
                    vertices_cc1(:, 1) - vertices_cc2(:, 1)];
                hp_b_cell{poly_id} = ...
                    vertices_cc2(:, 2) .* vertices_cc1(:, 1) ...
                    - vertices_cc1(:, 2) .* vertices_cc2(:, 1);

                % correct numerical inaccuracies and make sure that the
                % vertices themselves are contained in all of the
                % half-spaces
                hp_b_cell{poly_id} = max(hp_b_cell{poly_id}, ...
                    max(hp_w_cell{poly_id} * vertices_cc1', [], 2));

                % add redundant inequalities 0' * z <= 1 below such that
                % the number of inequalities for each polytope is the same
                hp_w_cell_rd{poly_id} = [hp_w_cell{poly_id}; 
                    zeros(vert_num_max - vert_num, 2)];
                hp_b_cell_rd{poly_id} = [hp_b_cell{poly_id}; 
                    ones(vert_num_max - vert_num, 1)];
            end

            % store the information about the hyperplanes characterizing
            % the polytopes
            obj.Quality.Hyperplane = struct;
            obj.Quality.Hyperplane.num = vert_num_max;
            obj.Quality.Hyperplane.w = vertcat(hp_w_cell_rd{:});
            obj.Quality.Hyperplane.b = vertcat(hp_b_cell_rd{:});

            % prepare quantities for the minimization of quadratic function
            % over the quality space
            obj.Storage.QuadMin = struct;

            convcomb_coef_cell = cell(poly_num, 1);
            convcomb_intercept_cell = cell(poly_num, 1);
            candidate_cell = cell(poly_num, 1);

            for poly_id = 1:poly_num
                % retrieve the vertices in the polytopes and their next
                % counterclockwise neighbor
                vertices_cc1 = obj.Quality.Polytopes{poly_id};
                vertices_cc2 = vertices_cc1([2:end, 1], :);

                % difference of the end points of each edge in x and y
                % coordinates
                vertices_diff_x = vertices_cc2(:, 1) - vertices_cc1(:, 1);
                vertices_diff_y = vertices_cc2(:, 2) - vertices_cc1(:, 2);

                % retrieve the inequalities characterizing the polytope
                hp_w = hp_w_cell{poly_id};
                hp_w_ss = sum(hp_w .^ 2, 2);
                hp_b = hp_b_cell{poly_id};

                % compute the coefficient and the intercept when computing
                % the weight of a point as the convex combination of two
                % end points of an edge of the polytope
                convcomb_coef_cell{poly_id} ...
                    = [hp_w(:, 2) .^ 2, -hp_w(:, 1) .* hp_w(:, 2)] ...
                    ./ hp_w_ss ./ vertices_diff_x;
                convcomb_intercept_cell{poly_id} ...
                    = (hp_w(:, 1) .* hp_b ./ hp_w_ss ...
                    - vertices_cc1(:, 1)) ./ vertices_diff_x;

                % if the edge is vertical, the method above will fail, and
                % we thus determine the weight using the y-coordinates
                convcomb_coef_y = [-hp_w(:, 1) .* hp_w(:, 2), ...
                    hp_w(:, 1) .^ 2] ./ hp_w_ss ./ vertices_diff_y;
                convcomb_intercept_y = (hp_w(:, 2) .* hp_b ./ hp_w_ss ...
                    - vertices_cc1(:, 2)) ./ vertices_diff_y;

                vertical_list = abs(vertices_diff_x) < 1e-12;
                convcomb_coef_cell{poly_id}(vertical_list, :) ...
                    = convcomb_coef_y(vertical_list, :);
                convcomb_intercept_cell{poly_id}(vertical_list) ...
                    = convcomb_intercept_y(vertical_list);

                % the candidate optimizer before scaling
                candidate_cell{poly_id} = hp_w ./ hp_w_ss;
            end

            obj.Storage.QuadMin.IneqWeights = vertcat(hp_w_cell{:});
            obj.Storage.QuadMin.IneqRHS = vertcat(hp_b_cell{:});
            obj.Storage.QuadMin.ConvCombCoefs ...
                = vertcat(convcomb_coef_cell{:});
            obj.Storage.QuadMin.ConvCombIntercepts ...
                = vertcat(convcomb_intercept_cell{:});
            obj.Storage.QuadMin.UnscaledEdgeCandidates ...
                = vertcat(candidate_cell{:});
            obj.Storage.QuadMin.VerticesSS ...
                = sum(obj.Quality.Vertices .^ 2, 2);
        end

        function inside = checkIfInsideQualitySpace(obj, pts, batch_size)
            % Check if the points are inside the quality space, which is
            % the union of polytopes. Compute in batches if necessary to
            % avoid excessive memory use
            % Inputs:
            %   pts: two-column matrix containing the input points
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
        
        function setSimplicialTestFuncs(obj, args_cell)
            % Set the simplicial test functions for all marginals at the
            % same time
            % Input:
            %   args_cell: cell array where each cell is a cell array
            %   containing all inputs to the method setSimplicialTestFuncs
            %   of each marginal

            for marg_id = 1:length(obj.MarginalWeights)
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

        function OT_info_cell = performReassembly(obj)
            % Perform reassembly by computing semi-discrete optimal
            % transport
            % Output:
            %   OT_info_cell: cell array where each cell is a cell array
            %   containing optimal transport-related information that can
            %   be saved and loaded later

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, ...
                    '--- semi-discrete OT starts ---\n');
            end
            
            marg_num = length(obj.MarginalWeights);

            for marg_id = 1:marg_num
                if isempty(obj.Options.OT.angle_num)
                    marg_angle_num = [];
                else
                    marg_angle_num = obj.Options.OT.angle_num(marg_id);
                end

                marg_pp_angle_indices = ...
                    obj.Options.OT.pp_angle_indices{marg_id};

                marg_options ...
                    = obj.Options.OT.optimization_options{marg_id};

                % the atoms and the corresponding probabilities of the
                % discretized marginal are exactly given by the test
                % functions and their respective integrals; this is only
                % valid due to the quadratic structure of the cost function
                marg = obj.Marginals{marg_id};
                marg_atoms = marg.SimplicialTestFuncs.Vertices;
                marg_probs = marg.SimplicialTestFuncs.Integrals;

                marg.computeOptimalTransport(marg_atoms, marg_probs, ...
                    marg_angle_num, [], marg_pp_angle_indices, ...
                    marg_options);

                if obj.Options.display
                    fprintf('%s: marginal %d done\n', ...
                        class(obj), marg_id);
                end

                % logging
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, ...
                        '%s: marginal %d done\n', ...
                        class(obj), marg_id);
                end

                % clear the temporary data to save memory
                marg.clearPreparation();
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, ...
                    '--- semi-discrete OT ends ---\n\n');
                fclose(log_file);
            end

            obj.Storage.OTComputed = true;

            OT_info_cell = obj.saveOptimalTransportInfo();
        end

        function OT_info_cell = saveOptimalTransportInfo(obj)
            % Retrieve computed semi-discrete optimal transport-related
            % information of each marginal
            % Output:
            %   OT_info_cell: cell array where each cell is a cell array
            %   containing optimal transport-related information that can
            %   be saved and loaded later

            marg_num = length(obj.MarginalWeights);
            OT_info_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                OT_info_cell{marg_id} = marg.saveOptimalTransportInfo();
            end
        end

        function loadOptimalTransportInfo(obj, OT_info_cell)
            % Load semi-discrete optimal transport-related information into
            % each marginal
            % Input:
            %   OT_info_cell: cell array where each cell is a cell array
            %   containing optimal transport-related information that can
            %   be saved and loaded later

            marg_num = length(obj.MarginalWeights);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg.loadOptimalTransportInfo(OT_info_cell{marg_id});
            end

            obj.Storage.OTComputed = true;
        end
    
        function [UB_disc_list, UB_cont_list, samps, samp_histpdf] ...
                = getMTUpperBoundsWRepetition(obj, samp_num, rep_num, ...
                rand_stream, batch_size, hist_edge_x, hist_edge_y)
            % Compute two upper bounds for the matching for teams problem 
            % with repetition of Monte Carlo integration.
            % Inputs:
            %   samp_num: number of samples for Monte Carlo integration
            %   rep_num: number of repetitions
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the
            %   continuous distribution on the quality space (default is
            %   1e4)
            %   hist_edge_x: edges of bins on the x-axis for 2D pdf
            %   estimation (default is [])
            %   hist_edge_y: edges of bins on the y-axis for 2D pdf
            %   estimation (default is [])
            % Output:
            %   UB_disc_list: vector containing the computed upper bounds 
            %   for the matching for teams problem based on the discrete 
            %   measure on the quality space
            %   UB_cont_list: vector containing the computed upper bounds 
            %   for the matching for teams problem based on the continuous 
            %   measure on the quality space
            %   samps: one particular set of Monte Carlo samples used in 
            %   the approximation of the bounds
            %   samp_histpdf: 2D pdf estimation via the histogram of the
            %   generated samples from the continuous quality distribution;
            %   only computed if both hist_edge_x and hist_edge_y are set

            if ~exist('rand_stream', 'var') ...
                    || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') ...
                    || isempty(batch_size)
                batch_size = 1e4;
            end

            if ~exist('hist_edge_x', 'var') || ~exist('hist_edge_y', 'var')
                hist_edge_x = [];
                hist_edge_y = [];
                samp_histpdf = [];
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

            % initialize the 2D histogram if it needs to be computed
            if ~isempty(hist_edge_x) && ~isempty(hist_edge_y)
                pdf_computed = true;

                samp_histpdf = zeros(length(hist_edge_x) - 1, ...
                    length(hist_edge_y) - 1);
            else
                pdf_computed = false;
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

                % update the 2D histogram
                if pdf_computed
                    samp_histpdf = samp_histpdf ...
                        + histcounts2(samps.ContinuousQualities(:, 1), ...
                        samps.ContinuousQualities(:, 2), ...
                        hist_edge_x, hist_edge_y, ...
                        'Normalization', 'pdf');
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, ...
                    '--- Monte Carlo sampling ends ---\n\n');
                fclose(log_file);
            end

            % normalize the 2D histogram
            if pdf_computed
                samp_histpdf = samp_histpdf / rep_num;
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
            samp_num, rand_stream, batch_size);
        
        % Generate independent random sample from the coupling between the 
        % continuous measure on the quality space and a continuous marginal
        coup = randSampleFromOptContinuousCoupling(obj, marg_id, ...
            samp_num, rand_stream, batch_size);
        
        % Compute two upper bounds for the matching for teams problem based 
        % on the discrete measure on the quality space and based on the 
        % continuous measure on the quality space. The bounds are 
        % approximated by Monte Carlo integration.
        [UB_disc, UB_cont, samps] = getMTUpperBounds(obj, samp_num, ...
            rand_stream, batch_size);

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
            % Check if the points are inside the quality space, which is
            % the union of polytopes
            % Input:
            %   pts: two-column matrix containing the input points
            % Output:
            %   inside: boolean vector indicating whether each of the
            %   points is inside the quality space

            pt_num = size(pts, 1);
            poly_num = length(obj.Quality.Polytopes);
            inside = any(reshape(all(reshape(obj.Quality.Hyperplane.w ...
                * pts' - obj.Quality.Hyperplane.b <= ...
                MT2DQuad.INSIDE_TOLERANCE, ...
                obj.Quality.Hyperplane.num, poly_num * pt_num), 1), ...
                poly_num, pt_num), 1)'; 
        end

        function initializeBeforeRun(obj)
            % Initialize the algorithm by computing some static quantities

            if ~obj.Storage.SimplicialTestFuncsInitialized
                obj.initializeSimplicialTestFuncs();
            end
        end

        function [min_pts, min_vals] = computeQuadMin(obj, weights, ...
                batch_size)
            % Solve the minimization of z' * z - 2 * weights' * z where z
            % is in the quality space. Computation is done in batches if
            % necessary to avoid excessive memory usage.
            % Inputs:
            %   weights: two-column matrix where each row contains a weight
            %   batch_size: the maximum number of inputs to be computed in
            %   a batch (default is 1e4)
            % Outputs:
            %   min_pts: the computed minimizers
            %   min_vals: the computed minimal values

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            weight_num = size(weights, 1);
            inside = obj.checkIfInsideQualitySpace(weights, batch_size);

            min_pts = zeros(weight_num, 2);
            min_vals = zeros(weight_num, 1);

            weights_inside = weights(inside, :);
            min_pts(inside, :) = weights_inside;
            min_vals(inside) = -sum(weights_inside .^ 2, 2);

            remain_indices = find(~inside);
            weight_remain = weights(remain_indices, :);
            remain_num = size(weight_remain, 1);
            
            if remain_num > 0
                batch_num = ceil(remain_num / batch_size);

                for batch_id = 1:batch_num
                    batch_weight_indices = ((batch_id - 1) ...
                        * batch_size + 1):min(batch_id ...
                        * batch_size, remain_num);
                    [min_pts(remain_indices(batch_weight_indices), :), ...
                        min_vals(remain_indices(batch_weight_indices))] ...
                        = obj.doComputeQuadMinOutside( ...
                        weight_remain(batch_weight_indices, :));
                end
            end
        end

        function [min_pts, min_vals] = doComputeQuadMinOutside(obj, ...
                weights)
            % Solve the minimization of z' * z - 2 * weights' * z where z
            % is in the quality space and weight is outside the quality
            % space (thus some projection is necessary)
            % Input:
            %   weights: two-column matrix where each row contains a weight
            % Outputs:
            %   min_pts: the computed minimizers
            %   min_vals: the computed minimal values

            weight_num = size(weights, 1);
            min_pts = zeros(weight_num, 2);
            
            % shifting the right-hand side values of the constraints
            % characterizing the polytopes and transform the problem into a
            % Euclidean norm minimization problem
            ineq_rhs_mat = obj.Storage.QuadMin.IneqRHS' ...
                - weights * obj.Storage.QuadMin.IneqWeights';
            edgecand_num = size(ineq_rhs_mat, 2);

            % express each candidate point as the affine combination of two
            % end points of an edge of a polytope (if weight is in [0, 1]
            % then the candidate is valid, else it is invalid)
            convcomb_mat = weights * obj.Storage.QuadMin.ConvCombCoefs' ...
                + obj.Storage.QuadMin.ConvCombIntercepts';

            % the x- and y-coordinates of the candidate points stored in
            % two matrices
            edgecand_x_mat_shifted = ineq_rhs_mat ...
                .* obj.Storage.QuadMin.UnscaledEdgeCandidates(:, 1)';
            edgecand_y_mat_shifted = ineq_rhs_mat ...
                .* obj.Storage.QuadMin.UnscaledEdgeCandidates(:, 2)';
            edgecand_x_mat = edgecand_x_mat_shifted + weights(:, 1);
            edgecand_y_mat = edgecand_y_mat_shifted + weights(:, 2);
            
            % the corresponding objective values
            edgecand_val_mat = edgecand_x_mat_shifted .^ 2 ...
                + edgecand_y_mat_shifted .^ 2 ...
                - sum(weights .^ 2, 2);

            % remove those invalid candidates
            edgecand_val_mat(convcomb_mat < 0 | convcomb_mat > 1) = inf;

            cand_val_mat = [edgecand_val_mat, ...
                obj.Storage.QuadMin.VerticesSS' ...
                - weights * (2 * obj.Quality.Vertices')];

            % choose a valid candidate with minimum objective value
            [min_vals, min_inds] = min(cand_val_mat, [], 2);
            edge_min_attained_list = min_inds <= edgecand_num;
            edge_min_attained_indices = find(edge_min_attained_list);
            min_attained_linind = sub2ind(size(edgecand_val_mat), ...
                edge_min_attained_indices, ...
                min_inds(edge_min_attained_indices));
            min_pts(edge_min_attained_indices, 1) ...
                = edgecand_x_mat(min_attained_linind);
            min_pts(edge_min_attained_indices, 2) ...
                = edgecand_y_mat(min_attained_linind);
            min_pts(~edge_min_attained_list, :) ...
                = obj.Quality.Vertices( ...
                min_inds(~edge_min_attained_list) - edgecand_num, :);
        end

        function weighted_sum = computeWeightedSumOfVertices(obj, ...
                indices, batch_size)
            % Compute weighted sum of indices in the triangulation of the
            % marginals, where the weights are specified in
            % obj.MarginalWeights. Computation is done in batches if
            % necessary. 
            % Inputs:
            %   indices: matrix of indices where each row corresponds to an
            %   input and each column corresponds to a marginal; each index
            %   corresponds to the index of a vertex in the triangulation
            %   of a marginal
            %   batch_size: the maximum number of inputs to be handled in a
            %   vectorized procedure (default is 1e4)
            % Output:
            %   weighted_sum: two-column matrix where each row represents a
            %   weighted sum

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            input_num = size(indices, 1);
            marg_num = length(obj.MarginalWeights);

            if input_num <= batch_size
                % if the number of input is less than the batch size, do
                % the computation directly

                % first add the offsets to the indices so that they become
                % row indices in the matrix containing all vertices
                vertex_indices = indices' + obj.Storage.MargVertNumOffsets;

                % sparse combination matrix to add up the coordinates of
                % the vertices
                comb_mat = sparse(repelem((1:input_num)', marg_num, 1), ...
                    (1:input_num * marg_num)', 1, input_num, ...
                    input_num * marg_num);

                weighted_sum = full(comb_mat ...
                    * obj.Storage.WeightedMargVertices(vertex_indices, :));
            else
                % if more than one batch is needed, pre-compute the
                % combination matrix
                batch_num = ceil(input_num / batch_size);
                weighted_sum_cell = cell(batch_num, 1);
                comb_mat = sparse(repelem((1:batch_size)', ...
                    marg_num, 1), ...
                    (1:batch_size * marg_num)', 1, batch_size, ...
                    batch_size * marg_num);

                for batch_id = 1:batch_num
                    if batch_id < batch_num
                        vertex_indices = indices(((batch_id - 1) ...
                            * batch_size + 1):batch_id ...
                            * batch_size, :)' ...
                            + obj.Storage.MargVertNumOffsets;
                        weighted_sum_cell{batch_id} ...
                            = full(comb_mat ...
                            * obj.Storage.WeightedMargVertices( ...
                            vertex_indices, :));
                    else
                        last_batch_size = input_num - (batch_num - 1) ...
                            * batch_size;
                        vertex_indices = indices(((batch_num - 1) ...
                            * batch_size + 1):end, :)' ...
                            + obj.Storage.MargVertNumOffsets;
                        weighted_sum_cell{batch_id} ...
                            = full(comb_mat(1:last_batch_size, ...
                            1:(last_batch_size * marg_num)) ...
                            * obj.Storage.WeightedMargVertices( ...
                            vertex_indices, :));
                    end
                end

                weighted_sum = vertcat(weighted_sum_cell{:});
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

