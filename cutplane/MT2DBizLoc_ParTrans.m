classdef MT2DBizLoc_ParTrans < LSIPMinCuttingPlaneAlgo
    %MT2DBIZLOC_PARTRANS Summary of this class goes here
    %   Detailed explanation goes here

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
        function obj = MT2DBizLoc_ParTrans(marginals, costfuncs, ...
                quality_vertices, quality_testfuncs, varargin)
            % Constructor method
            % Inputs:
            %   marginals: cell array containing objects of 
            %   ProbMeas1D_ConvexPolytope
            %   costfuncs: cell array containing structs with the following
            %   fields:
            %       weights: coefficients in front of the costs
            %       stations: two-column matrix indicating the train
            %       stations
            %       inter_station_costs: matrix containing the pairwise
            %       transportation costs between stations
            %       availabilities: logical vector indicating whether each
            %       agent population can use the train
            %   quality_vertices: two-column matrix containing the vertices
            %   of the quality space
            %   quality_testfuncs: cell array containing the inputs to the
            %   function obj.setQualitySimplicialTestFuncs

            obj@LSIPMinCuttingPlaneAlgo(varargin{:});

            marg_num = length(marginals);

            % set the default option for the formulation of global
            % minimization problem (options are 'DLOG' and 'CC')
            if ~isfield(obj.Options, 'global_formulation') ...
                    || isempty(obj.Options.global_formulation)
                obj.Options.global_formulation = 'DLOG';
            end

            % set the default option for the index of the discrete measure
            % candidate on the quality space
            if ~isfield(obj.Options, 'discmeas_cand_index') ...
                    || isempty(obj.Options.discmeas_cand_index)
                obj.Options.discmeas_cand_index = length(marginals);
            end

            % set the default option for sanitizing the entries in the
            % inequality constraints
            if ~isfield(obj.Options, 'sanitation_threshold') ...
                    || isempty(obj.Options.sanitation_threshold)
                obj.Options.sanitation_threshold = 0;
            end

            % set the default option for combining atoms on the support of
            % the marginals that only differ by numerical errors
            if ~isfield(obj.Options, 'marginal_magnetic_grids') ...
                    || isempty(obj.Options.marginal_magnetic_grids)
                obj.Options.marginal_magnetic_grids = ...
                    1e4 * ones(marg_num, 2);
            end

            % set the default options for semi-discrete optimal transport
            if ~isfield(obj.Options, 'OT') ...
                    || isempty(obj.Options.OT)
                obj.Options.OT = struct;
            end

            if ~isfield(obj.Options.OT, 'angle_num') ...
                    || isempty(obj.Options.OT.angle_num)
                obj.Options.OT.angle_num = [];
            elseif isscalar(obj.Options.OT.angle_num)
                obj.Options.OT.angle_num = ones(marg_num, 1) ...
                    * obj.Options.OT.angle_num;
            end

            if ~isfield(obj.Options.OT, 'pp_angle_indices') ...
                    || isempty(obj.Options.OT.pp_angle_indices)
                obj.Options.OT.pp_angle_indices = repmat({[]}, ...
                    marg_num, 1);
            elseif ~iscell(obj.Options.OT.pp_angle_indices)
                obj.Options.OT.pp_angle_indices = ...
                    repmat({obj.Options.OT.pp_angle_indices}, ...
                    marg_num, 1);
            end

            if ~isfield(obj.Options.OT, 'optimization_options')
                obj.Options.OT.optimization_options ...
                    = cell(marg_num, 1);
            elseif ~iscell(obj.Options.OT.optimization_options)
                optim_options = obj.Options.OT.optimization_options;
                obj.Options.OT.optimization_options ...
                    = cell(marg_num, 1);
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

            % set the marginals
            obj.Marginals = marginals;

            q_vert_ccorder = convhull(quality_vertices(:, 1), ...
                quality_vertices(:, 2), 'Simplify', true);
            q_vert_cc1 = quality_vertices(q_vert_ccorder(1:end - 1), :);
            q_vert_cc2 = quality_vertices(q_vert_ccorder(2:end), :);

            % compute the weights and the intercepts of the hyperplanes
            % characterizing the quality space
            hp_w = [q_vert_cc2(:, 2) - q_vert_cc1(:, 2), ...
                q_vert_cc1(:, 1) - q_vert_cc2(:, 1)];
            hp_b = q_vert_cc2(:, 2) .* q_vert_cc1(:, 1) ...
                - q_vert_cc1(:, 2) .* q_vert_cc2(:, 1);

            % correct numerical inaccuracies and make sure that the
            % vertices themselves are contained in all of the
            % half-spaces
            hp_b = max(hp_b, max(hp_w * q_vert_cc1', [], 2));

            % set the quality space
            obj.Quality = struct;

            obj.Quality.GurobiModel = struct;
            obj.Quality.GurobiModel.A = hp_w;
            obj.Quality.GurobiModel.rhs = hp_b;
            obj.Quality.GurobiModel.sense = repmat('<', length(hp_b), 1);
            obj.Quality.GurobiModel.lb = -inf(2, 1);
            obj.Quality.GurobiModel.ub = inf(2, 1);

            obj.Quality.Vertices = q_vert_cc1;

            % if the test functions are specified, set the test functions
            if exist('quality_testfuncs', 'var') ...
                    && ~isempty(quality_testfuncs)
                obj.setQualitySimplicialTestFuncs(quality_testfuncs{:});
            end

            % store a reference point from the quality space
            obj.Quality.SamplePoint = obj.Quality.Vertices(1, :)';

            % initialize the costs functions
            assert(length(costfuncs.weights) == marg_num, ...
                'weights in the cost function mis-specified');
            station_num = size(costfuncs.stations, 1);
            assert(size(costfuncs.stations, 2) == 2, ...
                'train station locations mis-specified');

            assert(size(costfuncs.inter_station_costs, 1) ...
                == station_num ...
                && size(costfuncs.inter_station_costs, 2) ...
                == station_num, ...
                'inter-station costs mis-specified');
            assert(all(all(costfuncs.inter_station_costs >= 0)) ...
                && issymmetric(costfuncs.inter_station_costs) ...
                && all(diag(costfuncs.inter_station_costs) == 0), ...
                ['inter-station costs must be non-negative, ' ...
                'symmetric, and the diagonal must all be 0']);

            assert(length(costfuncs.availabilities) == marg_num, ...
                'train availabilities mis-specified');

            obj.CostFuncs = struct;
            obj.CostFuncs.Weights = costfuncs.weights;
            obj.CostFuncs.StationNum = station_num;
            obj.CostFuncs.TrainStations = costfuncs.stations;
            obj.CostFuncs.InterStationCosts = ...
                costfuncs.inter_station_costs;
            obj.CostFuncs.TrainAvailabilities = costfuncs.availabilities;

            % save the Lipschitz constants of the cost functions
            obj.CostFuncs.LipschitzConstMarg = ...
                sqrt(2) * obj.CostFuncs.Weights;
            obj.CostFuncs.LipschitzConstQuality = ...
                sqrt(2) * obj.CostFuncs.Weights;

            % build the gurobi model used for computing the minimal cost
            % business outlet location given the locations of the employees
            % and the supplier
            lb_cell = cell(marg_num, 1);
            ub_cell = cell(marg_num, 1);
            vtype_cell = cell(marg_num, 1);
            obj_con = 0;

            % indices of variables representing the coordinates of the
            % business location
            quality_indices = [1; 2];
            var_counter = 2;
            quality_vtype = repmat('C', 2, 1);

            % indices of variables representing the costs contributed by
            % each marginal
            v_indices = var_counter + (1:marg_num)';
            var_counter = var_counter + marg_num;
            v_lb = -inf(marg_num, 1);
            v_ub = inf(marg_num, 1);
            v_vtype = repmat('C', marg_num, 1);

            % indices of variables representing the L1 distance between the
            % business location and each train station
            r_indices = reshape((1:station_num * 2)', 2, station_num)' ...
                + var_counter;
            var_counter = var_counter + 2 * station_num;
            r_lb = -inf(station_num * 2, 1);
            r_ub = inf(station_num * 2, 1);
            r_vtype = repmat('C', station_num * 2, 1);

            % indices of variables representing the L1 distance between the
            % business location and each of the input locations
            q_indices = zeros(marg_num, 2);

            % indices of variables representing the slack variables when
            % taking the minimum among distances
            s_indices = cell(marg_num, 1);

            % indices of binary variables connecting the slack variables
            % when taking the minimum among distances
            iota_indices = cell(marg_num, 1);

            for marg_id = 1:marg_num
                q_indices(marg_id, :) = var_counter + [1, 2];
                var_counter = var_counter + 2;
                lb_cell{marg_id} = -inf(2, 1);
                ub_cell{marg_id} = inf(2, 1);
                vtype_cell{marg_id} = repmat('C', 2, 1);

                if obj.CostFuncs.TrainAvailabilities(marg_id)
                    s_indices{marg_id} = var_counter ...
                        + (1:(station_num + 1))';
                    var_counter = var_counter + station_num + 1;
    
                    iota_indices{marg_id} = var_counter ...
                        + (1:(station_num + 1))';
                    var_counter = var_counter + station_num + 1;

                    lb_cell{marg_id} = [lb_cell{marg_id}; ...
                        zeros(station_num + 1, 1); ...
                        -inf(station_num + 1, 1)];
                    ub_cell{marg_id} = [ub_cell{marg_id}; ...
                        inf(station_num + 1, 1); ...
                        inf(station_num + 1, 1)];
                    vtype_cell{marg_id} = [vtype_cell{marg_id}; ...
                        repmat('C', station_num + 1, 1); ...
                        repmat('B', station_num + 1, 1)];
                end
            end
            
            % total number of decision variables;
            var_length = var_counter;

            % start to assemble the constraints
            constr_counter = 0;

            % the constraints on the business outlet location
            quality_constr_num = size(obj.Quality.GurobiModel.A, 1);
            A_quality = ...
                [sparse(quality_constr_num, quality_indices(1) - 1), ...
                obj.Quality.GurobiModel.A, ...
                sparse(quality_constr_num, ...
                var_length - quality_indices(end))];
            rhs_quality = obj.Quality.GurobiModel.rhs;
            sense_quality = obj.Quality.GurobiModel.sense;
            constr_counter = constr_counter + quality_constr_num;

            % the constraints that relate the business outlet location and
            % the r variables
            A_rcon_cell = cell(station_num, 1);
            rhs_rcon_cell = cell(station_num, 1);
            sense_rcon_cell = cell(station_num, 1);

            for station_id = 1:station_num
                A_rcon_cell{station_id} = sparse(repmat((1:4)', 2, 1), ...
                    [r_indices(station_id, :)'; ...
                    r_indices(station_id, :)'; ...
                    quality_indices; quality_indices], ...
                    [1; 1; 1; 1; -1; -1; 1; 1], ...
                    4, var_length);
                rhs_rcon_cell{station_id} = ...
                    [-obj.CostFuncs.TrainStations(station_id, :)'; ...
                    obj.CostFuncs.TrainStations(station_id, :)'];
                sense_rcon_cell{station_id} = repmat('>', 4, 1);

                % update the constraint counter
                constr_counter = constr_counter + 4;
            end

            A_rcon = vertcat(A_rcon_cell{:});
            rhs_rcon = vertcat(rhs_rcon_cell{:});
            sense_rcon = vertcat(sense_rcon_cell{:});

            % the constraints that link the business outlet location and
            % the q variables
            A_qcon_cell = cell(marg_num, 1);
            rhs_qcon_cell = cell(marg_num, 1);
            sense_qcon_cell = cell(marg_num, 1);

            % store the row indices of these constraints so that the input
            % locations can be filled into the right-hand side of these
            % constraints at runtime
            qcon_row_indices_cell = cell(marg_num, 1);

            % store the sparse matrix to link the right-hand side of these
            % constraints and the input locations; the input locations are
            % vectorized into [x1(1); x1(2); x2(1); x2(2); ... xN(1);
            % xN(2)]
            qcon_link_mat_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                A_qcon_cell{marg_id} = sparse(repmat((1:4)', 2, 1), ...
                    [q_indices(marg_id, :)'; q_indices(marg_id, :)'; ...
                    quality_indices; quality_indices], ...
                    [1; 1; 1; 1; -1; -1; 1; 1], ...
                    4, var_length);
                rhs_qcon_cell{marg_id} = zeros(4, 1);
                sense_qcon_cell{marg_id} = repmat('>', 4, 1);

                qcon_row_indices_cell{marg_id} = constr_counter ...
                    + (1:4)';
                constr_counter = constr_counter + 4;
                qcon_link_mat_cell{marg_id} = sparse((1:4)', ...
                    [1; 2; 1; 2], [-1; -1; 1; 1], 4, 2);
            end

            A_qcon = vertcat(A_qcon_cell{:});
            rhs_qcon = vertcat(rhs_qcon_cell{:});
            sense_qcon = vertcat(sense_qcon_cell{:});
            qcon_row_indices = vertcat(qcon_row_indices_cell{:});
            qcon_link_mat = blkdiag(qcon_link_mat_cell{:});

            % the constraints that enforces the v variables to be equal to
            % the minimum of distances
            A_mincon_cell = cell(marg_num, 1);
            rhs_mincon_cell = cell(marg_num, 1);
            sense_mincon_cell = cell(marg_num, 1);

            % store the row indices at which the distances between the
            % stations and the input locations need to be filled
            mincon_row_indices_cell = cell(marg_num, 1);

            % compute the big-M constant for linking the s variables and
            % the iota variables
            biz_max_L1_norm = max(sum(abs(obj.Quality.Vertices), 2));
            station_L1_norm = sum(abs(obj.CostFuncs.TrainStations), 2);
            inter_station_L1_dist = squareform(pdist( ...
                obj.CostFuncs.TrainStations, 'cityblock'));
            inter_station_L1_dist_vec = inter_station_L1_dist(:);
            mincon_bigM_cand1 = max(max(inter_station_L1_dist_vec ...
                + inter_station_L1_dist_vec' + abs( ...
                repmat(obj.CostFuncs.InterStationCosts, ...
                station_num, station_num) ...
                - repelem(obj.CostFuncs.InterStationCosts, ...
                station_num, station_num))));
            mincon_bigM_cand2 = 2 * biz_max_L1_norm ...
                + max(max(station_L1_norm + station_L1_norm' ...
                + obj.CostFuncs.InterStationCosts));
            mincon_bigM = max(mincon_bigM_cand1, mincon_bigM_cand2);

            for marg_id = 1:marg_num
                if ~obj.CostFuncs.TrainAvailabilities(marg_id)
                    A_mincon_cell{marg_id} = sparse(ones(3, 1), ...
                        [v_indices(marg_id); q_indices(marg_id, :)'], ...
                        [1, -1, -1], ...
                        1, var_length);
                    rhs_mincon_cell{marg_id} = zeros(1, 1);
                    sense_mincon_cell{marg_id} = repmat('=', 1, 1);
                    mincon_row_indices_cell{marg_id} = zeros(0, 1);
                    constr_counter = constr_counter + 1;

                    continue;
                end

                % part 1.1: linking the v variables, the q variables, and
                % the s variables
                A_mincon_11 = sparse(ones(4, 1), ...
                    [v_indices(marg_id); q_indices(marg_id, :)'; ...
                    s_indices{marg_id}(1)], ...
                    [1; -1; -1; 1], 1, var_length);
                rhs_mincon_11 = 0;
                sense_mincon_11 = '=';
                constr_counter = constr_counter + 1;

                % part 1.2: linking the v variables, the r variables, and
                % the s variables
                A_mincon_12 = sparse(repmat((1:station_num)', 4, 1), ...
                    [v_indices(marg_id) * ones(station_num, 1); ...
                    r_indices(:, 1); r_indices(:, 2); ...
                    s_indices{marg_id}(2:end)], ...
                    [ones(station_num, 1); ...
                    -ones(station_num, 1); ...
                    -ones(station_num, 1); ...
                    ones(station_num, 1)], ...
                    station_num, var_length);
                rhs_mincon_12 = zeros(station_num, 1);
                sense_mincon_12 = repmat('=', station_num, 1);

                mincon_row_indices_cell{marg_id} = constr_counter ...
                    + (1:station_num)';
                constr_counter = constr_counter + station_num;

                % part 2.1: linking the s variable and the iota variable
                % via the big-M constraints for the direct distance
                A_mincon_21 = sparse(ones(2, 1), ...
                    [s_indices{marg_id}(1); iota_indices{marg_id}(1)], ...
                    [1; mincon_bigM], 1, var_length);
                rhs_mincon_21 = mincon_bigM;
                sense_mincon_21 = '<';
                constr_counter = constr_counter + 1;

                % part 2.2: linking the s variables and the iota variables
                % via the big-M constraints for the indirect distances
                A_mincon_22 = sparse(repmat((1:station_num)', 2, 1), ...
                    [s_indices{marg_id}(2:end); ...
                    iota_indices{marg_id}(2:end)], ...
                    [ones(station_num, 1); ...
                    mincon_bigM * ones(station_num, 1)], ...
                    station_num, var_length);
                rhs_mincon_22 = mincon_bigM * ones(station_num, 1);
                sense_mincon_22 = repmat('<', station_num, 1);
                constr_counter = constr_counter + station_num;

                % part 3: requiring the iota variables to sum to 1
                A_mincon_3 = sparse(ones(station_num + 1, 1), ...
                    [iota_indices{marg_id}], ...
                    ones(station_num + 1, 1), ...
                    1, var_length);
                rhs_mincon_3 = 1;
                sense_mincon_3 = '=';
                constr_counter = constr_counter + 1;

                A_mincon_cell{marg_id} = ...
                    [A_mincon_11; A_mincon_12; ...
                    A_mincon_21; A_mincon_22;
                    A_mincon_3];
                rhs_mincon_cell{marg_id} = ...
                    [rhs_mincon_11; rhs_mincon_12; ...
                    rhs_mincon_21; rhs_mincon_22; ...
                    rhs_mincon_3];
                sense_mincon_cell{marg_id} = ...
                    [sense_mincon_11; sense_mincon_12; ...
                    sense_mincon_21; sense_mincon_22; ...
                    sense_mincon_3];
            end

            A_mincon = vertcat(A_mincon_cell{:});
            rhs_mincon = vertcat(rhs_mincon_cell{:});
            sense_mincon = vertcat(sense_mincon_cell{:});
            mincon_row_indices = vertcat(mincon_row_indices_cell{:});
            
            % store the linear indices for extracting the non-diagonal
            % extries in a matrix
            nondiag_indices = setdiff((1:station_num^2)', ...
                (1:(station_num + 1):station_num^2)');
            mincon_nondiag_indices_mat = nondiag_indices ...
                + ((0:(marg_num - 1)) * station_num^2);
            mincon_nondiag_indices_mat = mincon_nondiag_indices_mat( ...
                :, obj.CostFuncs.TrainAvailabilities);
            mincon_nondiag_indices = mincon_nondiag_indices_mat(:);


            % assemble the gurobi model for cost minimization
            min_cost_model = struct;
            min_cost_model.modelsense = 'min';
            min_cost_model.objcon = obj_con;
            min_cost_model.obj = zeros(var_length, 1);
            min_cost_model.obj(v_indices) = obj.CostFuncs.Weights;
            min_cost_model.lb = [obj.Quality.GurobiModel.lb; ...
                v_lb; r_lb; vertcat(lb_cell{:})];
            min_cost_model.ub = [obj.Quality.GurobiModel.ub; ...
                v_ub; r_ub; vertcat(ub_cell{:})];
            min_cost_model.vtype = [quality_vtype; v_vtype; r_vtype; ...
                vertcat(vtype_cell{:})];
            min_cost_model.A = ...
                [A_quality; A_rcon; A_qcon; A_mincon];
            min_cost_model.rhs = ...
                [rhs_quality; rhs_rcon; rhs_qcon; rhs_mincon];
            min_cost_model.sense = ...
                [sense_quality; sense_rcon; sense_qcon; sense_mincon];

            obj.Storage.MinCost = struct;
            obj.Storage.MinCost.GurobiModel = min_cost_model;
            obj.Storage.MinCost.MILPBigM = mincon_bigM;
            obj.Storage.MinCost.QConstraintsRowIndices = qcon_row_indices;
            obj.Storage.MinCost.QConstraintsLinkMat = qcon_link_mat;
            obj.Storage.MinCost.MinConstraintsRowIndices = ...
                mincon_row_indices;
            obj.Storage.MinCost.MinConstraintsNonDiagonalIndices = ...
                mincon_nondiag_indices;
            obj.Storage.MinCost.QualityIndices = quality_indices;

            % this flag is used to track if the function
            % obj.initializeSimplicialTestFuncs has been called
            obj.Storage.SimplicialTestFuncsInitialized = false;

            marg_testfunc_set = false(marg_num, 1);

            for marg_id = 1:marg_num
                marg_testfunc_set(marg_id) = ~isempty( ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs);
            end

            if isfield(obj.Quality, 'SimplicialTestFuncs') ...
                    && ~isempty(obj.Quality.SimplicialTestFuncs) ...
                    && all(marg_testfunc_set)
                % initialize the simplicial test functions at the end of
                % the constructor if they have already been set
                obj.initializeSimplicialTestFuncs();
            end
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
            %   marg_inputs: two-column matrix containing inputs in the
            %   support of the marginal
            %   quality_inputs: two-column matrix where each row is an
            %   input in the quality space
            %   batch_size: number of inputs to be handled together
            %   (default is 1e4)
            % Output:
            %   vals: vector containing the computed cost values

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            input_num = size(marg_inputs, 1);
            assert(size(quality_inputs, 1) == input_num, ...
                'input mis-specified');

            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);

            station_num = obj.CostFuncs.StationNum;
            inter_station_dist = obj.CostFuncs.InterStationCosts;

            for batch_id = 1:batch_num
                batch_list = ((batch_id - 1) * batch_size ...
                    + 1):min(batch_id * batch_size, input_num);
                batch_inputs = marg_inputs(batch_list, :);
                batch_qualities = quality_inputs(batch_list, :);
                batch_input_num = size(batch_inputs, 1);

                dist_direct = sum(abs(batch_inputs - batch_qualities), 2);
                
                if ~obj.CostFuncs.TrainAvailabilities(marg_id)
                    vals_cell{batch_id} = dist_direct;
                    continue;
                end

                dist_input2station = ...
                    pdist2(obj.CostFuncs.TrainStations, ...
                    batch_inputs, 'cityblock');
                dist_biz2station = pdist2(obj.CostFuncs.TrainStations, ...
                    batch_qualities, 'cityblock');

                dist_input2biz = reshape(dist_input2station, ...
                    [station_num, 1, batch_input_num]) ...
                    + reshape(dist_biz2station, ...
                    [1, station_num, batch_input_num]) ...
                    + inter_station_dist;

                vals_cell{batch_id} = min(dist_direct, reshape( ...
                    min(min(dist_input2biz, [], 1), [], 2), ...
                    [batch_input_num, 1, 1]));
            end

            vals = vertcat(vals_cell{:}) * obj.CostFuncs.Weights(marg_id);
        end

        function vals = evaluateSumOfCostFuncs(obj, marg_input_mat, ...
                quality_inputs, batch_size)
            % Evaluate the sum of cost functions. Computation is done in 
            % batches if necessary. Warning: this function does not check 
            % if the inputs are inside their corresponding domains. 
            % Inputs: 
            %   marg_input_mat: matrix where each row is an input and each
            %   two columns correspond to an input from a marginal
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
            assert(size(marg_input_mat, 2) == marg_num * 2, ...
                'input mis-specified');

            vals = zeros(input_num, 1);

            for marg_id = 1:marg_num
                vals = vals + obj.evaluateCostFunc(marg_id, ...
                    marg_input_mat(:, (marg_id - 1) * 2 + [1, 2]'), ...
                    quality_inputs, batch_size);
            end
        end

        function [min_pts, min_vals] = computeMinCostOverQuality(obj, ...
                input_mat)
            % Solve the minimization of the sum of cost functions where the
            % arguments in the support of marginals are given.
            % Inputs:
            %   input_mat: matrix where each row corresponds to an input 
            %   and every two columns correspond to an input from a
            %   marginal
            % Outputs:
            %   min_pts: the computed minimizers stored in a matrix where
            %   each row contains a minimizer
            %   min_vals: the computed minimal values

            marg_num = length(obj.Marginals);
            input_num = size(input_mat, 1);
            station_num = obj.CostFuncs.StationNum;
            available_num = sum(obj.CostFuncs.TrainAvailabilities);
            min_pts = zeros(input_num, 2);
            min_vals = zeros(input_num, 1);

            MIP_options = struct;
            MIP_options.OutputFlag = 0;
            MIP_options.IntFeasTol = 1e-6;
            MIP_options.FeasibilityTol = 1e-8;
            MIP_options.OptimalityTol = 1e-8;
            MIP_options.NodefileStart = 2;
            MIP_options.MIPGap = 1e-4;
            MIP_options.MIPGapAbs = 1e-10;

            qcon_row_indices = obj.Storage.MinCost.QConstraintsRowIndices;
            qcon_link_mat = obj.Storage.MinCost.QConstraintsLinkMat;
            mincon_row_indices = ...
                obj.Storage.MinCost.MinConstraintsRowIndices;
            mincon_nodiag_indices = ...
                obj.Storage.MinCost.MinConstraintsNonDiagonalIndices;
            quality_indices = obj.Storage.MinCost.QualityIndices;

            for input_id = 1:input_num
                input_vec = input_mat(input_id, :)';
                input_pts_mat = reshape(input_vec, 2, marg_num)';
                L1_dist_mat = pdist2(obj.CostFuncs.TrainStations, ...
                    input_pts_mat, 'cityblock');
                dist_mat = repelem(L1_dist_mat, 1, station_num) ...
                    + repmat(obj.CostFuncs.InterStationCosts, 1, marg_num);
                dist_mat_nodiag = ...
                    reshape(dist_mat(mincon_nodiag_indices), ...
                    station_num - 1, available_num * station_num);
                min_dist_mat = reshape(min(dist_mat_nodiag, [], 1), ...
                    station_num, available_num);
                model = obj.Storage.MinCost.GurobiModel;

                % fill in the right-hand side of the constraints related to
                % the q variables
                model.rhs(qcon_row_indices) = qcon_link_mat * input_vec;

                % complete the right-hand side of the constraints related
                % to the s variables
                model.rhs(mincon_row_indices) = min_dist_mat(:);
                    
                result = gurobi(model, MIP_options);

                if ~strcmp(result.status, 'OPTIMAL')
                    error('MIP solver failed');
                end

                min_pts(input_id, :) = result.x(quality_indices)';
                min_vals(input_id) = result.objval;
            end
        end

        function setSimplicialTestFuncs(obj, quality_args, args_cell)
            % Set the simplicial test functions on the quality space as
            % well as the simplicial test functions for all marginals at
            % the same time
            % Input:
            %   quality_arts: cell array containing all inputs to the
            %   method setQualitySimplicialTestFuncs
            %   args_cell: cell array where each cell is a cell array
            %   containing all inputs to the method setSimplicialTestFuncs
            %   of each marginal

            obj.setQualitySimplicialTestFuncs(quality_args{:});

            for marg_id = 1:length(obj.Marginals)
                obj.Marginals{marg_id}.setSimplicialTestFuncs( ...
                    args_cell{marg_id}{:});
            end

            % after setting the simplicial functions, initialize the
            % quantities for the cutting-plane algorithm
            obj.initializeSimplicialTestFuncs();
        end

        function setQualitySimplicialTestFuncs(obj, ...
                vertices, triangles, check_simpcover)
            % Initialize the test functions with respect to a simplicial
            % cover and compute the integrals of the test functions with
            % respect to the probability measure.
            % Warning: the function does not check if the union of the
            % triangles in the simplicial cover equals the quality space
            % given by the union of polytopes.
            % Inputs:
            %   vertices: two-column matrix containing the vertices used in
            %   the triangulation
            %   triangles: three-column matrix containing the triangulation
            %   check_simpcover: check whether the given triangulation
            %   forms a simplicial cover (default is true)

            if ~exist('check_simpcover', 'var') || isempty(check_simpcover)
                check_simpcover = true;
            end

            % check if all vertices are contained in the quality space
            assert(all(obj.checkIfInsideQualitySpace(vertices)), ...
                'there are vertices outside the quality space');

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

            obj.Quality.SimplicialTestFuncs = struct;
            obj.Quality.SimplicialTestFuncs.Vertices = vertices;
            obj.Quality.SimplicialTestFuncs.Triangles = triangles;

            tri_num = size(triangles, 1);

            % compute the inverse transformation matrix which transform a
            % coordinate into weights in the affine combination
            invtrans_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert = vertices(triangles(tri_id, :), :);

                % compute the inverse transformation matrix
                invtrans_cell{tri_id} = [tri_vert'; ones(1, 3)] ...
                    \ eye(3);
            end

            obj.Quality.SimplicialTestFuncs.InvTransMat ...
                = vertcat(invtrans_cell{:});

            % compute the mesh size
            v1 = vertices(triangles(:, 1), :);
            v2 = vertices(triangles(:, 2), :);
            v3 = vertices(triangles(:, 3), :);

            edge1sq = sum((v1 - v2) .^ 2, 2);
            edge2sq = sum((v2 - v3) .^ 2, 2);
            edge3sq = sum((v3 - v1) .^ 2, 2);

            obj.Quality.SimplicialTestFuncs.MeshSize = ...
                sqrt(max(max(max(edge1sq, edge2sq), edge3sq)));
        end

        function [vals, inside] = evaluateQualitySimplicialTestFuncs( ...
                obj, pts, batch_size)
            % Evaluate the test functions on the quality space at given
            % locations; the input locations must be inside the quality
            % space
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

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);
            inside_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                [vals_cell{batch_id}, inside_cell{batch_id}] ...
                    = obj.doEvaluateQualitySimplicialTestFuncs( ...
                    pts(((batch_id - 1) * batch_size + 1):min(batch_id ...
                    * batch_size, input_num), :));
            end

            vals = vertcat(vals_cell{:});
            inside = vertcat(inside_cell{:});
        end
        
        function coup_cell = updateSimplicialTestFuncs(obj, ...
                quality_args, args_cell, fixed_index)
            % Update the simplicial test functions after an execution of
            % the cutting-plane algorithm. Besides setting the new
            % simplicial test functions, a new set of couplings of 
            % discretized marginals are generated via reassembly of the 
            % dual solution from the cutting-plane algorithm with the new 
            % discretized marginals. These couplings can be used to 
            % generate initial constraints for the new LSIP problem with 
            % the updated test functions.
            % Input:
            %   quality_args: cell array containing all inputs to the
            %   method setQualitySimplicialTestFuncs
            %   args_cell: cell array where each cell is a cell array
            %   containing all inputs to the methods setSimplicialTestFuncs
            %   of each marginal
            %   fixed_index: index of the marginal whose corresponding
            %   measure on the quality space is used in the updated
            %   couplings (default is 1)
            % Outputs:
            %   coup_cell: cell array containing structs with fields
            %   encoding the couplings between a marginal and the measure
            %   on the quality space

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('fixed_index', 'var') || isempty(fixed_index)
                fixed_index = 1;
            end

            marg_num = length(obj.Marginals);

            % retrieve the dual solution resulted from the cutting-plane
            % algorithm
            dual_sol = obj.Runtime.DualSolution;
            old_coup_probs_cell = cell(marg_num, 1);
            old_coup_indices_cell = cell(marg_num, 1);
            old_coup_marg_atoms_cell = cell(marg_num, 1);
            old_coup_marg_probs_cell = cell(marg_num, 1);
            old_coup_q_atoms_cell = cell(marg_num, 1);
            old_coup_q_probs_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                input_atoms = dual_sol{marg_id}.InputAtoms;
                quality_atoms = dual_sol{marg_id}.QualityAtoms;
                old_coup_probs_cell{marg_id} = ...
                    dual_sol{marg_id}.Probabilities;

                % store only the unique atoms in the discretized marginals
                % and use their indices in the coupling with the quality
                % space
                [~, ulist_marg, uind_marg] ...
                    = unique(round(input_atoms, 6), 'rows', 'stable');
                old_coup_marg_atoms_cell{marg_id} = ...
                    input_atoms(ulist_marg, :);
                old_coup_marg_probs_cell{marg_id} = accumarray( ...
                    uind_marg, old_coup_probs_cell{marg_id});

                % store only the unique atoms in the quality space and
                % use their indices in the coupling with the marginal
                [~, ulist_quality, uind_quality] ...
                    = unique(round(quality_atoms, 6), 'rows', 'stable');
                old_coup_q_atoms_cell{marg_id} = ...
                    quality_atoms(ulist_quality, :);
                old_coup_q_probs_cell{marg_id} = accumarray( ...
                    uind_quality, old_coup_probs_cell{marg_id});
                
                old_coup_indices_cell{marg_id} = ...
                    [uind_marg, uind_quality];

                if marg_id == fixed_index
                    keep_list = old_coup_q_probs_cell{marg_id} >= 1e-12;
                    fixed_cand_q_probs = ...
                        old_coup_q_probs_cell{marg_id}(keep_list);
                    fixed_cand_q_probs = fixed_cand_q_probs ...
                        / sum(fixed_cand_q_probs);
                    fixed_cand_q_atoms = ...
                        old_coup_q_atoms_cell{marg_id}(keep_list, :);
                end
            end

            % set the new simplicial test functions
            obj.setSimplicialTestFuncs(quality_args, args_cell);

            % move the atoms in the fixed candidate on the quality space to
            % the closest vertex in the updated test functions; then,
            % take a discrete measure on the quality space that is the
            % mixture of the fixed candidate and a discrete measure that is
            % uniform on all vertices in the new triangulation
            fixed_meas_interp_const = 0.05;
            fixed_dist_mat = pdist2(fixed_cand_q_atoms, ...
                obj.Quality.SimplicialTestFuncs.Vertices, 'euclidean');
            [~, min_dist_ind] = min(fixed_dist_mat, [], 2);
            fixed_meas_q_atoms = obj.Quality.SimplicialTestFuncs.Vertices;
            fixed_meas_q_atom_num = size(fixed_meas_q_atoms, 1);
            fixed_meas_q_probs = ones(fixed_meas_q_atom_num, 1) ...
                / fixed_meas_q_atom_num * fixed_meas_interp_const;
            fixed_meas_q_probs_additional = accumarray(min_dist_ind, ...
                fixed_cand_q_probs, [fixed_meas_q_atom_num, 1]);
            fixed_meas_q_probs = fixed_meas_q_probs ...
                + fixed_meas_q_probs_additional ...
                * (1 - fixed_meas_interp_const);
            fixed_meas_q_probs = fixed_meas_q_probs ...
                / sum(fixed_meas_q_probs);

            coup_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                new_atoms = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices;
                new_probs = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;

                % compute an optimal coupling between the original discrete
                % marginal and the new discrete marginal discrete optimal
                % transport

                % the cost function is the Euclidean distance
                dist_mat = pdist2(old_coup_marg_atoms_cell{marg_id}, ...
                    new_atoms, 'euclidean');
                if size(dist_mat, 1) * size(dist_mat, 2) > 1e6
                    % if there are too many atoms in the discrete measures,
                    % a direct computation of discrete OT may cause memory
                    % throttling; thus, we resort to a constraint
                    % generation scheme
                    cp_options = struct('display', false);
                    OT = OTDiscrete(old_coup_marg_probs_cell{marg_id}, ...
                        new_probs, dist_mat, cp_options);
                    [hcoup_indices, ~] = ...
                        OT.generateHeuristicCoupling();
                    OT.run(hcoup_indices, 1e-6);
                    coup = OT.Runtime.DualSolution;
                    marg_coup_atom_indices = coup.CoupIndices;
                    marg_coup_probs = coup.Probabilities;
                else
                    [marg_coup_atom_indices, marg_coup_probs] ...
                        = discrete_OT( ...
                        old_coup_marg_probs_cell{marg_id}, ...
                        new_probs, dist_mat);
                end

                % compute an optimal coupling between the measure on
                % the quality space and the fixed measure on the
                % quality space; this is to ensure that all measures on
                % the quality space satisfy the constraints with
                % respect to the test functions
                q_dist_mat = pdist2(old_coup_q_atoms_cell{marg_id}, ...
                    fixed_meas_q_atoms, 'euclidean');

                if size(q_dist_mat, 1) * size(q_dist_mat, 2) > 1e6
                    % if there are too many atoms in the discrete measures,
                    % a direct computation of discrete OT may cause memory
                    % throttling; thus, we resort to a constraint
                    % generation scheme
                    cp_options = struct('display', false);
                    OT = OTDiscrete(old_coup_q_probs_cell{marg_id}, ...
                        fixed_meas_q_probs, q_dist_mat, cp_options);
                    [hcoup_indices, ~] = ...
                        OT.generateHeuristicCoupling();
                    OT.run(hcoup_indices, 1e-6);
                    coup = OT.Runtime.DualSolution;
                    q_coup_atom_indices = coup.CoupIndices;
                    q_coup_probs = coup.Probabilities;
                else
                    [q_coup_atom_indices, q_coup_probs] = discrete_OT( ...
                        old_coup_q_probs_cell{marg_id}, ...
                        fixed_meas_q_probs, q_dist_mat);
                end

                % perform discrete reassembly to get the new coupling;
                % in this case, we need to reassemble both marginals in
                % the coupling
                [coup_indices, coup_probs] = discrete_reassembly( ...
                    old_coup_indices_cell{marg_id}, ...
                    old_coup_probs_cell{marg_id}, ...
                    {marg_coup_atom_indices; q_coup_atom_indices}, ...
                    {marg_coup_probs; q_coup_probs}, 1:2);
                coup_qualities = fixed_meas_q_atoms( ...
                    coup_indices(:, 2), :);
                coup_inputs = new_atoms(coup_indices(:, 1), :);
                coup_atom_num = length(coup_probs);
                marg_testfunc_vals = sparse(1:coup_atom_num, ...
                    coup_indices(:, 1), 1, coup_atom_num, ...
                    length(new_atoms));
                quality_testfunc_vals = sparse(1:coup_atom_num, ...
                    coup_indices(:, 2), 1, coup_atom_num, ...
                    fixed_meas_q_atom_num);

                % compute the corresponding cost function values
                coup_costfunc_vals = obj.evaluateCostFunc(marg_id, ...
                    coup_inputs, coup_qualities);

                coup_cell{marg_id} = struct( ...
                    'probabilities', coup_probs, ...
                    'inputs', coup_inputs, ...
                    'qualities', coup_qualities, ...
                    'marg_testfunc_vals', marg_testfunc_vals, ...
                    'quality_testfunc_vals', quality_testfunc_vals, ...
                    'costfunc_vals', coup_costfunc_vals);
            end
        end

        function initializeSimplicialTestFuncs(obj)
            % Initialize some quantities related to the simplicial test
            % functions of the marginals and the simplicial test functions
            % on the quality space

            marg_num = length(obj.Marginals);

            quality_testfuncs = obj.Quality.SimplicialTestFuncs;
            quality_vert_num = size(quality_testfuncs.Vertices, 1);
            quality_tri_num = size(quality_testfuncs.Triangles, 1);
            obj.Storage.QualityTriNum = quality_tri_num;

            % store the number of vertices in the test functions for each
            % marginal
            obj.Storage.MargVertNumList = zeros(marg_num, 1);

            % store the indices of the decision variables that correspond
            % to the intercepts, the coefficients of the test functions for
            % the marginals, and the coefficients of the test functions on
            % the quality space
            obj.Storage.DeciVarInterceptIndices = zeros(marg_num, 1);
            obj.Storage.DeciVarMargTestFuncIndices = cell(marg_num, 1);
            obj.Storage.DeciVarQualityTestFuncIndices = cell(marg_num, 1);

            ind_counter = 0;

            for marg_id = 1:marg_num
                obj.Storage.DeciVarInterceptIndices(marg_id) = ...
                    ind_counter + 1;
                ind_counter = ind_counter + 1;

                marg = obj.Marginals{marg_id};
                obj.Storage.MargVertNumList(marg_id) = length( ...
                    marg.SimplicialTestFuncs.Vertices);

                % note that the first test function for the marginal is
                % removed for identification purposes
                obj.Storage.DeciVarMargTestFuncIndices{marg_id} = ...
                    ind_counter ...
                    + (1:obj.Storage.MargVertNumList(marg_id) - 1)';
                ind_counter = ind_counter ...
                    + obj.Storage.MargVertNumList(marg_id) - 1;

                % note that the first test function on the quality space is
                % removed for identification purposes
                obj.Storage.DeciVarQualityTestFuncIndices{marg_id} = ...
                    ind_counter + (1:quality_vert_num - 1)';
                ind_counter = ind_counter + quality_vert_num - 1;
            end

            obj.Storage.DeciVarLength = ind_counter;
            obj.Storage.TotalMargVertNum ...
                = sum(obj.Storage.MargVertNumList);
            
            % initialize the global minimization models
            obj.initializeGlobalMinModels();

            % initialize the minimization models for evaluating the
            % transfer functions
            obj.initializeTransFuncModels();

            obj.Storage.SimplicialTestFuncsInitialized = true;

            % updating the simplicial test functions will invalidate all
            % quantities in the runtime environment, thus all variables in
            % the runtime environment need to be flushed
            obj.Runtime = [];
        end

        function coup_cell = generateDiscreteCoupling(obj)
            % Generate a feasible dual solution by solving the discretized
            % version of the problem
            % Output:
            %   coup_cell: cell array containing structs with fields
            %   encoding the couplings between a marginal and the measure
            %   on the quality space

            marg_num = length(obj.Marginals);
            quality_vertices = obj.Quality.SimplicialTestFuncs.Vertices;
            quality_vert_num = size(quality_vertices, 1);
            marg_vert_num_list = obj.Storage.MargVertNumList;
            decivar_num = sum(marg_vert_num_list) ...
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
                marg_vertices = marg.SimplicialTestFuncs.Vertices;
                marg_integral = marg.SimplicialTestFuncs.Integrals;
                marg_vert_num = marg_vert_num_list(marg_id);
                
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
                    marg_vert_num), speye(quality_vert_num - 1)];

                % matrix for the inequality constraints will have all
                % possible combinations of the vertices in the test
                % functions for the marginal and the vertices in the test
                % functions on the quality space
                [marg_grid, q_grid] = meshgrid(1:marg_vert_num, ...
                    1:quality_vert_num);
                marg_grid_indices = marg_grid(:);
                q_grid_indices = q_grid(:);
                A_marg_full = sparse( ...
                    1:marg_vert_num * quality_vert_num, ...
                    marg_grid_indices, 1, ...
                    marg_vert_num * quality_vert_num, marg_vert_num);
                A_quality_full = sparse( ...
                    1:marg_vert_num * quality_vert_num, ...
                    q_grid_indices, 1, ...
                    marg_vert_num * quality_vert_num, quality_vert_num);

                % remove the first test function for the marginal and the
                % first test function on the quality space; also prepend a
                % column of 1s
                A_ineq_cell{marg_id} = [sparse(ones(marg_vert_num ...
                    * quality_vert_num, 1)), ...
                    A_marg_full(:, 2:end), ...
                    A_quality_full(:, 2:end)];

                 % the right-hand side of the inequality constraints are
                % computed from the corresponding coordinates of the
                % vertices
                marg_grid_pts = marg_vertices(marg_grid_indices, :);
                quality_grid_pts = quality_vertices(q_grid_indices, :);
                rhs_ineq_cell{marg_id} = ...
                    obj.evaluateCostFunc(marg_id, marg_grid_pts, ...
                    quality_grid_pts);

                constr_num_list(marg_id) = marg_vert_num ...
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

            % the equality constraints
            A_eq = horzcat(A_eq_cell{:});
            rhs_eq = zeros(quality_vert_num - 1, 1);

            % assemble the constraints in the model
            model.A = [A_eq; blkdiag(A_ineq_cell{:})];
            model.rhs = [rhs_eq; vertcat(rhs_ineq_cell{:})];
            model.sense = [repmat('=', length(rhs_eq), 1); ...
                repmat('<', length(model.rhs) - length(rhs_eq), 1)];

            % parameters of the LP solver
            LP_options = struct;
            LP_options.OutputFlag = 0;

            result = gurobi(model, LP_options);

            if ~strcmp(result.status, 'OPTIMAL')
                error('unexpected error in LP');
            end

            probs = result.pi(quality_vert_num:end);
            coup_cell = cell(marg_num, 1);
            constr_counter = 0;

            % unwrap the dual optimal solution
            for marg_id = 1:marg_num
                coup_probs = probs(constr_counter ...
                    + (1:constr_num_list(marg_id)));

                % retain only those atoms with positive probabilities
                pos_list = coup_probs > 0;
                coup_probs = coup_probs(pos_list);
                coup_atom_num = length(coup_probs);
                coup_marg_points = ...
                    marg_points_cell{marg_id}(pos_list, :);
                coup_quality_indices = ...
                    quality_indices_cell{marg_id}(pos_list);
                coup_quality_points = ...
                    quality_points_cell{marg_id}(pos_list, :);
                marg_testfunc_vals = sparse(1:coup_atom_num, ...
                    marg_indices_cell{marg_id}(pos_list), 1, ...
                    coup_atom_num, marg_vert_num_list(marg_id));
                quality_testfunc_vals = sparse(1:coup_atom_num, ...
                    coup_quality_indices, 1, coup_atom_num, ...
                    quality_vert_num);
                coup_cost = rhs_ineq_cell{marg_id}(pos_list, :);

                coup_cell{marg_id} = struct( ...
                    'probabilities', coup_probs, ...
                    'inputs', coup_marg_points, ...
                    'qualities', coup_quality_points, ...
                    'marg_testfunc_vals', marg_testfunc_vals, ...
                    'quality_testfunc_vals', quality_testfunc_vals, ...
                    'costfunc_vals', coup_cost);

                % update the constraint counter
                constr_counter = constr_counter + constr_num_list(marg_id);
            end
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

        function prepareDiscreteOT(obj)
            % Compute discrete to discrete optimal transport between the
            % candidate discrete measure on the quality space and the rest
            % of the discrete measures on the quality space

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            cand_id = obj.Options.discmeas_cand_index;
            cand = obj.Runtime.DualSolution{cand_id};
            [~, cand_ulist, cand_uind] = unique( ...
                round(cand.QualityAtoms, 4), 'rows', 'stable');
            cand_atoms = cand.QualityAtoms(cand_ulist, :);
            cand_atom_num = size(cand_atoms, 1);
            cand_probs = accumarray(cand_uind, cand.Probabilities, ...
                [cand_atom_num, 1]);

            marg_num = length(obj.Marginals);

            marg_atoms_cell = cell(marg_num, 1);
            coup_cell = cell(marg_num, 1);
            discreteOT_cost_list = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                dual_sol = obj.Runtime.DualSolution{marg_id};
                [~, m_ulist, m_uind] = unique( ...
                    round(dual_sol.InputAtoms ...
                    .* obj.Options.marginal_magnetic_grids(marg_id, ...
                    :)), 'rows', 'stable');
                marg_atoms = dual_sol.InputAtoms(m_ulist, :);
                marg_atoms = marg.sanitizePoints(marg_atoms, 1e-6);

                if any(pdist(marg_atoms) < 1e-3)
                    error(['some atoms are too close, which will ' ...
                        'cause numerical issues; consider changing ' ...
                        'the option marginal_magnetic_grids']);
                end

                marg_atoms_cell{marg_id} = marg_atoms;
                marg_atom_num = size(marg_atoms, 1);

                if marg_id == cand_id
                    % simply store the coupling between the candidate
                    % measure and the discretized marginal
                    coup_cell{marg_id} = sparse(cand_uind, ...
                        m_uind, cand.Probabilities, ...
                        cand_atom_num, marg_atom_num);
                    discreteOT_cost_list(marg_id) = 0;

                    continue;
                end

                [~, q_ulist, q_uind] = unique( ...
                    round(dual_sol.QualityAtoms, 4), 'rows', 'stable');
                quality_atoms = dual_sol.QualityAtoms(q_ulist, :);
                quality_atom_num = size(quality_atoms, 1);

                % construct the sparse matrix representing the coupling
                % between the discrete measure on the quality space and the
                % discretized marginal
                q2m_coup_mat = sparse(q_uind, m_uind, ...
                    dual_sol.Probabilities, quality_atom_num, ...
                    marg_atom_num);

                % compute the optimal coupling between the candidate
                % discrete measure and this discrete measure on the quality
                % space
                dist_mat = pdist2(cand_atoms, quality_atoms, 'euclidean');
                quality_probs = accumarray(q_uind, ...
                    dual_sol.Probabilities, [quality_atom_num, 1]);
                
                if size(dist_mat, 1) * size(dist_mat, 2) > 1e6
                    % if there are too many atoms in the discrete measures,
                    % a direct computation of discrete OT may cause memory
                    % throttling; thus, we resort to a constraint
                    % generation scheme
                    cp_options = struct('display', false);
                    OT = OTDiscrete(cand_probs, quality_probs, ...
                        dist_mat, cp_options);
                    [hcoup_indices, ~] = ...
                        OT.generateHeuristicCoupling();
                    OT.run(hcoup_indices, 1e-6);
                    coup = OT.Runtime.DualSolution;
                    coup_atoms = coup.CoupIndices;
                    coup_probs = coup.Probabilities;
                    discreteOT_cost_list(marg_id) = -OT.Runtime.LSIP_LB;
                else
                    [coup_atoms, coup_probs, ...
                        discreteOT_cost_list(marg_id)] ...
                        = discrete_OT(cand_probs, quality_probs, dist_mat);
                end
                
                % construct the sparse matrix representing the coupling
                % between the candidate discrete measure and this discrete
                % measure on the quality space
                q2q_coup_mat = sparse(coup_atoms(:, 1), ...
                    coup_atoms(:, 2), coup_probs, cand_atom_num, ...
                    quality_atom_num);

                % the coupling between the discrete candidate measure and
                % the discretized marginal is formed by composing the two
                % couplings
                q2q_cond_mat_sum = sum(q2q_coup_mat, 1);
                q2q_cond_mat_sum(q2q_cond_mat_sum == 0) = 1; %#ok<SPRIX>
                q2q_cond_mat = q2q_coup_mat ./ q2q_cond_mat_sum;
                coup_cell{marg_id} = q2q_cond_mat * q2m_coup_mat;
            end

            % store the candidate discrete measure on the quality space
            obj.Runtime.DiscreteQualityCandidate = struct( ...
                'Probabilities', cand_probs, ...
                'Atoms', cand_atoms);

            % store the computed discrete OT costs in the runtime
            % environment
            obj.Runtime.DiscreteOTCosts = discreteOT_cost_list;

            % store the computed couplings in the runtime environment
            obj.Runtime.DiscretizedMarginalAtoms = marg_atoms_cell;
            obj.Runtime.DiscreteCouplings = coup_cell;

            obj.Runtime.DiscreteCouplingsComputed = true;
        end

        function OT_info = prepareSemiDiscreteOT(obj)
            % Compute semi-discrete optimal transport between the
            % discretized marginals and the continuous marginals
            % Output:
            %   OT_info: optimal transport-related information returned by
            %   saveOptimalTransportInfo()

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~isfield(obj.Runtime, 'DiscreteCouplingsComputed') ...
                    || ~obj.Runtime.DiscreteCouplingsComputed
                error(['need to first compute the discrete couplings ' ...
                    'by calling prepareDiscreteOT()']);
            end

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, ...
                    '--- semi-discrete OT starts ---\n');
            end
            
            marg_num = length(obj.Marginals);

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
                % discretized marginal are taken from the computed
                % approximately optimal solution of the dual of the LSIP
                % problem (duplicate atoms are already combined)
                marg = obj.Marginals{marg_id};
                marg_atoms = obj.Runtime.DiscretizedMarginalAtoms{marg_id};
                marg_probs = full(sum( ...
                    obj.Runtime.DiscreteCouplings{marg_id}, 1)');

                marg.computeOptimalTransport(marg_atoms, marg_probs, ...
                    marg_angle_num, [], marg_pp_angle_indices, ...
                    marg_options);

                if obj.Options.display
                    fprintf('%s: marginal %d done\n', class(obj), marg_id);
                end

                % logging
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, ...
                        '%s: marginal %d done\n', class(obj), marg_id);
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, ...
                    '--- semi-discrete OT ends ---\n\n');
                fclose(log_file);
            end

            obj.Runtime.OTComputed = true;

            OT_info = obj.saveOptimalTransportInfo();
        end

        function OT_info = saveOptimalTransportInfo(obj)
            % Retrieve computed semi-discrete optimal transport-related
            % information of each marginal as well as the couplings between
            % the discretized marginals and the candidate measure on the
            % quality space
            % Output:
            %   OT_info: struct with fields semi_discrete and discrete

            if ~isfield(obj.Runtime, 'OTComputed') ...
                    || ~obj.Runtime.OTComputed
                error(['must first compute discrete and semi-discrete ' ...
                    'optimal transport']);
            end

            marg_num = length(obj.Marginals);
            semidisc_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                semidisc_cell{marg_id} = marg.saveOptimalTransportInfo();
            end

            OT_info = struct;
            OT_info.semi_discrete = semidisc_cell;
            OT_info.discrete = struct;
            OT_info.discrete.candidate = ...
                obj.Runtime.DiscreteQualityCandidate;
            OT_info.discrete.costs = obj.Runtime.DiscreteOTCosts;
            OT_info.discrete.marginal_atoms = ...
                obj.Runtime.DiscretizedMarginalAtoms;
            OT_info.discrete.couplings = obj.Runtime.DiscreteCouplings;
        end

        function loadOptimalTransportInfo(obj, OT_info)
            % Load semi-discrete optimal transport-related information into
            % each marginal and load the couplings between the discretized
            % marginals and the candidate measure on the quality space
            % Input:
            %   OT_info: struct produced by saveOptimalTransportInfo()

            marg_num = length(obj.Marginals);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg.loadOptimalTransportInfo( ...
                    OT_info.semi_discrete{marg_id});
            end

            obj.Runtime.DiscreteQualityCandidate = ...
                OT_info.discrete.candidate;
            obj.Runtime.DiscreteOTCosts = OT_info.discrete.costs;
            obj.Runtime.DiscretizedMarginalAtoms = ...
                OT_info.discrete.marginal_atoms;
            obj.Runtime.DiscreteCouplings = OT_info.discrete.couplings;

            obj.Runtime.DiscreteCouplingsComputed = true;
            obj.Runtime.OTComputed = true;
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
            %   batch_size: number of samples to be handled simultaneously
            %   during the computation of the sum of the cost functions 
            %   (default is 1e4)
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

        function vals = evaluateOptParametricFunc(obj, ...
                marg_id, pts, batch_size)
            % Evaluate the parametric function on the support of a marginal
            % parametrized by simplicial test functions that is optimized
            % using the cutting plane algorithm. This corresponds to the
            % primal solution of the LSIP problem. The computation is done
            % is batches if necessary.
            % Inputs: 
            %   marg_id: the index of the marginal
            %   pts: two-column matrix containing the input points as rows
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

            marg = obj.Marginals{marg_id};
            primal_sol = obj.Runtime.PrimalSolution{marg_id};
            vals = marg.evaluateWeightedSumOfSimplicialTestFuncs( ...
                pts, primal_sol.Coefficients, batch_size) ...
                + primal_sol.Constant;
        end

        function vals_mat = evaluateOptTransferFuncs(obj, pts, ref_pt)
            % Evaluate the transfer functions resulted from optimized
            % parametric functions. These transfer functions is part of
            % approximate matching equilibria. 
            % Inputs:
            %   pts: two-column matrix containing the input points as rows
            %   ref_pt: vector indicating a reference point where the 
            %   transfer function will evaluate to 0 (default is the sample
            %   point from the quality space)
            % Output:
            %   vals_mat: matrix where each column contains the computed
            %   transfer function corresponding to a marginal at the input
            %   points

            if ~isfield(obj.Runtime, 'PrimalSolution') ...
                    || isempty(obj.Runtime.PrimalSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('ref_pt', 'var') || isempty(ref_pt)
                ref_pt = obj.Quality.SamplePoint;
            else
                assert(obj.checkIfInsideQualitySpace(ref_pt'), ...
                    'the reference point is not in the quality space');
            end

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, ...
                    '--- Computation of transfer functions starts ---\n');
            end

            marg_num = length(obj.Marginals);
            stations = obj.CostFuncs.TrainStations;
            station_num = size(stations, 1);
            inter_station_costs = obj.CostFuncs.InterStationCosts;
            tr_avail = obj.CostFuncs.TrainAvailabilities;

            nondiag_list = setdiff((1:station_num^2)', ...
                (1:station_num + 1:station_num^2)');

            MIP_options = struct;
            MIP_options.OutputFlag = 0;
            MIP_options.IntFeasTol = 1e-6;
            MIP_options.FeasibilityTol = 1e-8;
            MIP_options.OptimalityTol = 1e-8;
            MIP_options.NodefileStart = 2;
            MIP_options.MIPGap = 1e-4;
            MIP_options.MIPGapAbs = 1e-10;

            % add the reference point
            pts = [ref_pt'; pts];
            pt_num = size(pts, 1);

            assert(all(obj.checkIfInsideQualitySpace(pts)), ...
                'some points are outside the quality space');

            vals_mat = zeros(pt_num, marg_num);

            for marg_id = 1:marg_num - 1
                tf = obj.Storage.TransFunc{marg_id};
                primal_sol = obj.Runtime.PrimalSolution{marg_id};

                for pt_id = 1:pt_num
                    quality_pt = pts(pt_id, :)';
                    model = tf.GurobiModel;
                    model.objcon = -primal_sol.Constant;
                    model.obj(tf.MargTFCoefsIndices) = ...
                        tf.MargTFCoefsLinkMat * primal_sol.Coefficients;
                    model.rhs(tf.QConstraintsRowIndices) = ...
                        tf.QConstraintsLinkMat * quality_pt;

                    dist_sta = (sum(abs(quality_pt' - stations), 2)' ...
                        + inter_station_costs)';
                    dist_sta = reshape(dist_sta(nondiag_list), ...
                        station_num - 1, station_num);

                    if tr_avail(marg_id)
                        model.rhs(tf.MinConstraintsRowIndices) = ...
                            min(dist_sta, [], 1);
                    end

                    result = gurobi(model, MIP_options);

                    if ~strcmp(result.status, 'OPTIMAL')
                        error('MIP solver failed');
                    end

                    vals_mat(pt_id, marg_id) = result.objval;

                    if mod(pt_id, 1e3) == 0
                        % display output
                        if obj.Options.display
                            fprintf(['%s: ' ...
                                'Transfer function for marginal %2d, ' ...
                                '%5d points computed\n'], ...
                                class(obj), marg_id, pt_id);
                        end
    
                        % write log
                        if ~isempty(obj.Options.log_file)
                            fprintf(log_file, ['%s: ' ...
                                'Transfer function for marginal %2d, ' ...
                                '%5d points computed\n'], ...
                                class(obj), marg_id, pt_id);
                        end
                    end
                end

                % display output
                if obj.Options.display
                    fprintf(['%s: ' ...
                        'Transfer function for marginal %2d ' ...
                        'computed\n'], ...
                        class(obj), marg_id);
                end

                % write log
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, ['%s: ' ...
                        'Transfer function for marginal %2d ' ...
                        'computed\n'], ...
                        class(obj), marg_id);
                end
            end

            % subtract the values evaluated at the reference point to
            % guarantee that all the transfer functions will evaluate to 0
            % at the reference point
            vals_mat = vals_mat - vals_mat(1, :);

            % the last transfer function is used to guarantee that the sum
            % of all transfer functions is equal to 0
            vals_mat(:, marg_num) = -sum(vals_mat(:, 1:marg_num - 1), 2);

            % remove the reference point
            vals_mat = vals_mat(2:end, :);

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, ...
                    '--- Computation of transfer functions  ends ---\n\n');
                fclose(log_file);
            end
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

            if ~isfield(obj.Runtime, 'DiscreteCouplingsComputed') ...
                    || ~obj.Runtime.DiscreteCouplingsComputed
                obj.prepareDiscreteOT();
            end

            distri = obj.Runtime.DiscreteQualityCandidate;
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
                obj.Runtime.DiscreteQualityCandidate.Atoms( ...
                atom_indices, :);

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
            %       Qualities: two-column matrix containing the coordinates
            %       of samples from the discrete measure on the quality
            %       space 
            %       AgentTypes: two-column matrix containing the
            %       coordinates of samples from the continuous marginal

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
                obj.Runtime.DiscreteQualityCandidate.Atoms( ...
                atom_indices, :);
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
            %   ql_samps: two-column matrix containing the coordinates of
            %   samples from the discrete measure on the quality space

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
            %       Qualities: two-column matrix containing the coordinates
            %       of samples from the continuous measure on the quality
            %       space 
            %       AgentTypes: two-column matrix containing the
            %       coordinates of samples from the continuous marginal

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
                obj.Runtime.DiscreteQualityCandidate.Atoms( ...
                atom_indices, :);

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

            % make sure that the couplings between the candidate discrete 
            % measure and the discretized marginals are computed
            if ~isfield(obj.Runtime, 'DiscreteCouplingsComputed') ...
                    || ~obj.Runtime.DiscreteCouplingsComputed
                obj.prepareDiscreteOT();
            end

            marg_num = length(obj.Marginals);
            EB = obj.Runtime.LSIP_UB - obj.Runtime.LSIP_LB;

            for marg_id = 1:marg_num
                EB = EB + sqrt(2) * obj.CostFuncs.Weights(marg_id) ...
                    * (obj.Runtime.DiscreteOTCosts(marg_id) ...
                    + obj.Marginals{marg_id}.OT.Cost);
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

            cf = obj.CostFuncs;

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                Lip = sqrt(2) * cf.Weights(marg_id);
                EB = EB + Lip * 2 * marg.SimplicialTestFuncs.MeshSize;
                
                if marg_id ~= obj.Options.discmeas_cand_index
                    EB = EB + Lip ...
                        * 2 * obj.Quality.SimplicialTestFuncs.MeshSize;
                end
            end
        end
    end

    methods(Access = protected)

        function initializeGlobalMinModels(obj)
            % Initialize the gurobi models for solving the global
            % minimization problems in the cutting-plane algorithm

            marg_num = length(obj.Marginals);

            quality_testfuncs = obj.Quality.SimplicialTestFuncs;
            quality_vert_num = size(quality_testfuncs.Vertices, 1);
            quality_tri_num = size(quality_testfuncs.Triangles, 1);

            cf_weights = obj.CostFuncs.Weights;
            station_num = obj.CostFuncs.StationNum;
            station_comb_num = station_num * (station_num - 1);
            stations = obj.CostFuncs.TrainStations;
            intsta_cost_mat = obj.CostFuncs.InterStationCosts;
            tr_avail = obj.CostFuncs.TrainAvailabilities;
            nondiag_indices = setdiff((1:station_num^2)', ...
                (1:(station_num + 1):station_num^2)');
            nondiag_row_indices = repmat((1:station_num)', ...
                1, station_num);
            nondiag_row_indices = nondiag_row_indices(nondiag_indices);
            nondiag_col_indices = repmat(1:station_num, ...
                station_num, 1);
            nondiag_col_indices = nondiag_col_indices(nondiag_indices);
            mincon_bigM = obj.Storage.MinCost.MILPBigM;

            obj.Storage.GlobalMin = cell(marg_num, 1);

            if strcmp(obj.Options.global_formulation, 'CC')
                % build the gurobi models for the global minimization
                % problems; the two-dimensional CPWA functions are modeled 
                % by the convex combination (CC) formulation
                global_form_CC = true;
            elseif strcmp(obj.Options.global_formulation, 'DLOG')
                % build the gurobi models for the global minimization 
                % problems; the two-dimensional CPWA functions are modeled 
                % by the logarithmic disaggregated convex combination
                % (DLOG) formulation
                global_form_CC = false;

                qtf_bit_length = ceil(log2(quality_tri_num));
                qtf_bisect_cell = cell(qtf_bit_length, 2);

                for bit_id = 1:qtf_bit_length
                    logi0 = bitget((0:quality_tri_num - 1)', bit_id) == 0;
                    qtf_bisect_cell{bit_id, 1} = find(logi0);
                    qtf_bisect_cell{bit_id, 2} = find(~logi0);
                end
            else
                error(['unknown formulation for the global ' ...
                    'minimization problem']);
            end

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_tf = marg.SimplicialTestFuncs;
                marg_tf_vert_num = size(marg_tf.Vertices, 1);
                marg_tf_tri_num = size(marg_tf.Triangles, 1);

                % first compute the indices of the variables
                input_indices = (1:2)';
                input_lb = -inf(2, 1);
                input_ub = inf(2, 1);
                input_vtype = repmat('C', 2, 1);
                quality_indices = (3:4)';
                quality_lb = -inf(2, 1);
                quality_ub = inf(2, 1);
                quality_vtype = repmat('C', 2, 1);

                % variable that represents the value of the cost
                % function
                xi_index = 5;
                var_counter = 5;
                xi_lb = -inf;
                xi_ub = inf;
                xi_vtype = 'C';

                if tr_avail(marg_id)
                    % variables that represent the L1 distance between
                    % the business location and the train stations
                    r_indices = reshape((1:2 * station_num)', ...
                        2, station_num)' + var_counter;
                    var_counter = var_counter + 2 * station_num;
                    r_lb = -inf(2 * station_num, 1);
                    r_ub = inf(2 * station_num, 1);
                    r_vtype = repmat('C', 2 * station_num, 1);

                    % variables that represent the L1 distance between
                    % the employee location and the train stations
                    p_indices = reshape((1:2 * station_num)', ...
                        2, station_num)' + var_counter;
                    var_counter = var_counter + 2 * station_num;
                    p_lb = -inf(2 * station_num, 1);
                    p_ub = inf(2 * station_num, 1);
                    p_vtype = repmat('C', 2 * station_num, 1);
                end

                % variables that represent the L1 distance between the
                % employee and the business location
                q_indices = var_counter + (1:2)';
                var_counter = var_counter + 2;
                q_lb = -inf(2, 1);
                q_ub = inf(2, 1);
                q_vtype = repmat('C', 2, 1);

                if tr_avail(marg_id)
                    % variable that represents the residual distance
                    % with respect to the direct distance
                    s_indices0 = var_counter + 1;
                    var_counter = var_counter + 1;

                    % variables that represent the residual distances
                    % with respect to the indirect distances
                    s_indices = var_counter + (1:station_comb_num)';
                    var_counter = var_counter + station_comb_num;

                    s_lb = zeros(station_comb_num + 1, 1);
                    s_ub = inf(station_comb_num + 1, 1);
                    s_vtype = repmat('C', station_comb_num + 1, 1);

                    % binary variable that enforces the minimum
                    % relation between the xi variable and the q
                    % variables
                    iota_indices0 = var_counter + 1;
                    var_counter = var_counter + 1;

                    % binary variables that enforces the minimum
                    % relations between the xi variable, the r
                    % variables, and the p variables
                    iota_indices = var_counter + (1:station_comb_num)';
                    var_counter = var_counter + station_comb_num;

                    iota_lb = -inf(station_comb_num + 1, 1);
                    iota_ub = inf(station_comb_num + 1, 1);
                    iota_vtype = repmat('B', station_comb_num + 1, 1);

                    cf_lb = [xi_lb; r_lb; p_lb; q_lb; s_lb; iota_lb];
                    cf_ub = [xi_ub; r_ub; p_ub; q_ub; s_ub; iota_ub];
                    cf_vtype = [xi_vtype; r_vtype; p_vtype; q_vtype; ...
                        s_vtype; iota_vtype];
                else
                    cf_lb = [xi_lb; q_lb];
                    cf_ub = [xi_ub; q_ub];
                    cf_vtype = [xi_vtype; q_vtype];
                end

                if global_form_CC
                    % variables between 0 and 1 that are used to represent 
                    % a point in the marginal support as a convex
                    % combination 
                    beta_indices = var_counter + (1:marg_tf_vert_num)';
                    var_counter = var_counter + marg_tf_vert_num;
                    beta_lb = zeros(marg_tf_vert_num, 1);
                    beta_ub = inf(marg_tf_vert_num, 1);
                    beta_vtype = repmat('C', marg_tf_vert_num, 1);

                    % transformation matrix applied to the vector
                    % containing coefficients of the test functions for the
                    % marginal before inserting the values into the
                    % objective vector of the mixed-integer program
                    marg_tf_coef_link_mat = -speye(marg_tf_vert_num);

                    % indices of locations in the objective vector of the
                    % mixed-integer program where the transformed
                    % coefficients of the test functions for the marginal
                    % will be inserted
                    marg_tf_coef_indices = beta_indices;

                    % binary-valued variables for indicating the triangle 
                    % for the parametric function on the marginal support
                    chi_indices = var_counter + (1:marg_tf_tri_num)';
                    var_counter = var_counter + marg_tf_tri_num;
                    chi_lb = -inf(marg_tf_tri_num, 1);
                    chi_ub = inf(marg_tf_tri_num, 1);
                    chi_vtype = repmat('B', marg_tf_tri_num, 1);

                    % variables between 0 and 1 that are used to represent 
                    % a point in the quality space as a convex combination
                    gamma_indices = var_counter + (1:quality_vert_num)';
                    var_counter = var_counter + quality_vert_num;
                    gamma_lb = zeros(quality_vert_num, 1);
                    gamma_ub = inf(quality_vert_num, 1);
                    gamma_vtype = repmat('C', quality_vert_num, 1);

                    % transformation matrix applied to the vector
                    % containing coefficients of the test functions on the
                    % quality space before inserting the values into the 
                    % objective vector of the mixed-integer program
                    quality_tf_coef_link_mat = -speye(quality_vert_num);

                    % indices of locations in the objective vector of the
                    % mixed-integer program where the transformed
                    % coefficients of the test functions on the quality
                    % space will be inserted
                    quality_tf_coef_indices = gamma_indices;

                    % binary-valued variables for indicating the triangle 
                    % for the parametric transfer function on the quality 
                    % space
                    psi_indices = var_counter + (1:quality_tri_num)';
                    var_counter = var_counter + quality_tri_num;
                    psi_lb = -inf(quality_tri_num, 1);
                    psi_ub = inf(quality_tri_num, 1);
                    psi_vtype = repmat('B', quality_tri_num, 1);

                    tf_lb = [beta_lb; chi_lb; gamma_lb; psi_lb];
                    tf_ub = [beta_ub; chi_ub; gamma_ub; psi_ub];
                    tf_vtype = [beta_vtype; chi_vtype; ...
                        gamma_vtype; psi_vtype];
                else
                    mtf_bit_length = ceil(log2(marg_tf_tri_num));
                    mtf_bisect_cell = cell(mtf_bit_length, 2);

                    for bit_id = 1:mtf_bit_length
                        logi0 = bitget((0:marg_tf_tri_num - 1)', ...
                            bit_id) == 0;
                        mtf_bisect_cell{bit_id, 1} = find(logi0);
                        mtf_bisect_cell{bit_id, 2} = find(~logi0);
                    end

                    % variables between 0 and 1 that are used to represent
                    % a point in a triangle in the marginal support as a
                    % convex combination
                    beta_indices = var_counter ...
                        + reshape((1:3 * marg_tf_tri_num)', ...
                        3, marg_tf_tri_num)';
                    var_counter = var_counter + 3 * marg_tf_tri_num;
                    beta_lb = zeros(3 * marg_tf_tri_num, 1);
                    beta_ub = inf(3 * marg_tf_tri_num, 1);
                    beta_vtype = repmat('C', 3 * marg_tf_tri_num, 1);

                    % transformation matrix applied to the vector
                    % containing coefficients of the test functions for the
                    % marginal before inserting the values into the
                    % objective vector of the mixed-integer program
                    marg_tf_coef_link_mat = sparse( ...
                        (1:3 * marg_tf_tri_num)', ...
                        marg_tf.Triangles(:), ...
                        -1, 3 * marg_tf_tri_num, marg_tf_vert_num);

                    % indices of locations in the objective vector of the
                    % mixed-integer program where the transformed
                    % coefficients of the test functions for the marginal
                    % will be inserted
                    marg_tf_coef_indices = beta_indices(:);

                    % binary-valued variables for indicating the triangle
                    % for the parametric transfer function on the marginal
                    % support
                    chi_indices = var_counter + (1:mtf_bit_length)';
                    var_counter = var_counter + mtf_bit_length;
                    chi_lb = -inf(mtf_bit_length, 1);
                    chi_ub = inf(mtf_bit_length, 1);
                    chi_vtype = repmat('B', mtf_bit_length, 1);

                    % variables between 0 and 1 that are used to represent
                    % a point in a triangle in the quality space as a
                    % convex combination
                    gamma_indices = var_counter ...
                        + reshape((1:3 * quality_tri_num)', ...
                        3, quality_tri_num)';
                    var_counter = var_counter + 3 * quality_tri_num;
                    gamma_lb = zeros(3 * quality_tri_num, 1);
                    gamma_ub = inf(3 * quality_tri_num, 1);
                    gamma_vtype = repmat('C', 3 * quality_tri_num, 1);

                    % transformation matrix applied to the vector
                    % containing coefficients of the test functions on the
                    % quality space before inserting the values into the 
                    % objective vector of the mixed-integer program
                    quality_tf_coef_link_mat = sparse( ...
                        (1:3 * quality_tri_num)', ...
                        quality_testfuncs.Triangles(:), ...
                        -1, 3 * quality_tri_num, quality_vert_num);

                    % indices of locations in the objective vector of the
                    % mixed-integer program where the transformed
                    % coefficients of the test functions on the quality
                    % space will be inserted
                    quality_tf_coef_indices = gamma_indices(:);

                    % binary-valued variables for indicating the triangle
                    % for the parametric transfer function on the quality
                    % space
                    psi_indices = var_counter + (1:qtf_bit_length)';
                    var_counter = var_counter + qtf_bit_length;
                    psi_lb = -inf(qtf_bit_length, 1);
                    psi_ub = inf(qtf_bit_length, 1);
                    psi_vtype = repmat('B', qtf_bit_length, 1);

                    tf_lb = [beta_lb; chi_lb; gamma_lb; psi_lb];
                    tf_ub = [beta_ub; chi_ub; gamma_ub; psi_ub];
                    tf_vtype = [beta_vtype; chi_vtype; ...
                        gamma_vtype; psi_vtype];
                end

                var_length = var_counter;

                % assemble the constraint matrix
                if tr_avail(marg_id)
                    % the constraints that relate the business outlet
                    % location and the r variables
                    A_rcon_cell = cell(station_num, 1);
                    rhs_rcon_cell = cell(station_num, 1);
                    sense_rcon_cell = cell(station_num, 1);

                    for station_id = 1:station_num
                        A_rcon_cell{station_id} = ...
                            sparse(repmat((1:4)', 2, 1), ...
                            [r_indices(station_id, :)'; ...
                            r_indices(station_id, :)'; ...
                            quality_indices; quality_indices], ...
                            [1; 1; 1; 1; -1; -1; 1; 1], ...
                            4, var_length);
                        rhs_rcon_cell{station_id} = ...
                            [-stations(station_id, :)'; ...
                            stations(station_id, :)'];
                        sense_rcon_cell{station_id} = ...
                            repmat('>', 4, 1);
                    end

                    A_rcon = vertcat(A_rcon_cell{:});
                    rhs_rcon = vertcat(rhs_rcon_cell{:});
                    sense_rcon = vertcat(sense_rcon_cell{:});

                    % the constraints that relate the employee
                    % location and the p variables
                    A_pcon_cell = cell(station_num, 1);
                    rhs_pcon_cell = cell(station_num, 1);
                    sense_pcon_cell = cell(station_num, 1);

                    for station_id = 1:station_num
                        A_pcon_cell{station_id} = ...
                            sparse(repmat((1:4)', 2, 1), ...
                            [p_indices(station_id, :)'; ...
                            p_indices(station_id, :)'; ...
                            input_indices; input_indices], ...
                            [1; 1; 1; 1; -1; -1; 1; 1], ...
                            4, var_length);
                        rhs_pcon_cell{station_id} = ...
                            [-stations(station_id, :)'; ...
                            stations(station_id, :)'];
                        sense_pcon_cell{station_id} = ...
                            repmat('>', 4, 1);
                    end

                    A_pcon = vertcat(A_pcon_cell{:});
                    rhs_pcon = vertcat(rhs_pcon_cell{:});
                    sense_pcon = vertcat(sense_pcon_cell{:});
                end

                % the constraints that relate the employee location,
                % the business outlet location, and the q variables
                A_qcon = sparse(repmat((1:4)', 3, 1), ...
                    [q_indices; q_indices; ...
                    quality_indices; quality_indices; ...
                    input_indices; input_indices], ...
                    [1; 1; 1; 1; -1; -1; 1; 1; 1; 1; -1; -1], ...
                    4, var_length);
                rhs_qcon = zeros(4, 1);
                sense_qcon = repmat('>', 4, 1);

                if tr_avail(marg_id)
                    % the following part will encode the constraints
                    % that enforce the minimum relation in the cost
                    % function

                    % part 1.1: linking the xi variable, the q
                    % variables, and the s variable for the direct
                    % distance
                    A_cfmin_11 = sparse(ones(4, 1), ...
                        [xi_index; q_indices; s_indices0], ...
                        [1; -1; -1; 1], 1, var_length);
                    rhs_cfmin_11 = 0;
                    sense_cfmin_11 = '=';

                    % part 1.2: linking the xi variable, the r
                    % variables, the p variables, and the s variables
                    % for the indirect distances
                    A_cfmin_12 = sparse( ...
                        repmat((1:station_comb_num)', 6, 1), ...
                        [xi_index * ones(station_comb_num, 1); ...
                        r_indices(nondiag_col_indices, 1); ...
                        r_indices(nondiag_col_indices, 2); ...
                        p_indices(nondiag_row_indices, 1); ...
                        p_indices(nondiag_row_indices, 2); ...
                        s_indices], ...
                        [ones(station_comb_num, 1); ...
                        -ones(station_comb_num, 1); ...
                        -ones(station_comb_num, 1); ...
                        -ones(station_comb_num, 1); ...
                        -ones(station_comb_num, 1); ...
                        ones(station_comb_num, 1)], ...
                        station_comb_num, var_length);
                    rhs_cfmin_12 = intsta_cost_mat(nondiag_indices);
                    sense_cfmin_12 = repmat('=', station_comb_num, 1);

                    % part 2.1: linking the s variable and the iota
                    % variable in the direct distance case via the
                    % big-M constraint
                    A_cfmin_21 = sparse(ones(2, 1), ...
                        [s_indices0; iota_indices0], ...
                        [1; mincon_bigM], ...
                        1, var_length);
                    rhs_cfmin_21 = mincon_bigM;
                    sense_cfmin_21 = '<';

                    % part 2.2: linking the s variables and the iota
                    % variables in the indirect distance case via the
                    % big-M constraints
                    A_cfmin_22 = sparse( ...
                        repmat((1:station_comb_num)', 2, 1), ...
                        [s_indices; iota_indices], ...
                        [ones(station_comb_num, 1); ...
                        mincon_bigM * ones(station_comb_num, 1)], ...
                        station_comb_num, var_length);
                    rhs_cfmin_22 = mincon_bigM ...
                        * ones(station_comb_num, 1);
                    sense_cfmin_22 = repmat('<', station_comb_num, 1);

                    % part 3: requiring all iota variables to sum to 1
                    A_cfmin_3 = sparse( ...
                        ones(station_comb_num + 1, 1), ...
                        [iota_indices0; iota_indices], ...
                        ones(station_comb_num + 1, 1), ...
                        1, var_length);
                    rhs_cfmin_3 = 1;
                    sense_cfmin_3 = '=';

                    A_cfmin = [A_cfmin_11; A_cfmin_12; ...
                        A_cfmin_21; A_cfmin_22; ...
                        A_cfmin_3];
                    rhs_cfmin = [rhs_cfmin_11; rhs_cfmin_12; ...
                        rhs_cfmin_21; rhs_cfmin_22; ...
                        rhs_cfmin_3];
                    sense_cfmin = [sense_cfmin_11; sense_cfmin_12; ...
                        sense_cfmin_21; sense_cfmin_22; ...
                        sense_cfmin_3];

                    A_cf = [A_rcon; A_pcon; A_qcon; A_cfmin];
                    rhs_cf = [rhs_rcon; rhs_pcon; rhs_qcon; rhs_cfmin];
                    sense_cf = [sense_rcon; sense_pcon; sense_qcon; ...
                        sense_cfmin];
                else
                    A_cfmin = sparse(ones(3, 1), ...
                        [xi_index; q_indices], ...
                        [1; -1; -1], ...
                        1, var_length);
                    rhs_cfmin = 0;
                    sense_cfmin = '=';

                    A_cf = [A_qcon; A_cfmin];
                    rhs_cf = [rhs_qcon; rhs_cfmin];
                    sense_cf = [sense_qcon; sense_cfmin];
                end

                if global_form_CC
                    % constraint that the sum of all beta coefficients must 
                    % be equal to 1
                    A_mtf_sum1 = sparse( ...
                        ones(marg_tf_vert_num, 1), beta_indices, ...
                        ones(marg_tf_vert_num, 1), 1, var_length);
                    rhs_mtf_sum1 = 1;
                    sense_mtf_sum1 = '=';

                    % constraint that the sum of all chi variables must be
                    % equal to 1
                    A_mtf_sum2 = sparse( ...
                        ones(marg_tf_tri_num, 1), chi_indices, ...
                        ones(marg_tf_tri_num, 1), 1, var_length);
                    rhs_mtf_sum2 = 1;
                    sense_mtf_sum2 = '=';

                    % constraint identifying non-zero beta's with chi
                    A_mtf_id_cell = cell(marg_tf_vert_num, 1);

                    for vert_id = 1:marg_tf_vert_num
                        % list of triangles that contain this vertex
                        tri_contain_list = find(any( ...
                            marg_tf.Triangles == vert_id, 2));
                        tri_contain_num = length(tri_contain_list);

                        A_mtf_id_cell{vert_id} = sparse( ...
                            ones(tri_contain_num + 1, 1), ...
                            [beta_indices(vert_id); ...
                            chi_indices(tri_contain_list)], ...
                            [1; -ones(tri_contain_num, 1)], ...
                            1, var_length);
                    end

                    A_mtf_id = vertcat(A_mtf_id_cell{:});
                    rhs_mtf_id = zeros(marg_tf_vert_num, 1);
                    sense_mtf_id = repmat('<', marg_tf_vert_num, 1);

                    % constraint linking the beta variables and the point
                    % in the marginal support
                    A_mtf_link = sparse( ...
                        repelem([1; 2], marg_tf_vert_num + 1, 1), ...
                        [beta_indices; input_indices(1); ...
                        beta_indices; input_indices(2)], ...
                        [marg_tf.Vertices(:, 1); -1; ...
                        marg_tf.Vertices(:, 2); -1], ...
                        2, var_length);
                    rhs_mtf_link = zeros(2, 1);
                    sense_mtf_link = repmat('=', 2, 1);

                    % sparse matrix used to directly extract the input
                    % components of the computed (approximate) minimizers;
                    % used to improve numerical stability since the input
                    % components of the computed (approximate) minimizers
                    % might slightly violate the equality constraints
                    input_extract_mat = sparse( ...
                        repelem([1; 2], marg_tf_vert_num, 1), ...
                        [beta_indices; beta_indices], ...
                        [marg_tf.Vertices(:, 1); ...
                        marg_tf.Vertices(:, 2)], ...
                        2, var_length);

                    % assemble all constraints related to beta variables
                    % and chi variables
                    A_mtf = [A_mtf_sum1; A_mtf_sum2; ...
                        A_mtf_id; A_mtf_link];
                    rhs_mtf = [rhs_mtf_sum1; rhs_mtf_sum2; ...
                        rhs_mtf_id; rhs_mtf_link];
                    sense_mtf = [sense_mtf_sum1; sense_mtf_sum2; ...
                        sense_mtf_id; sense_mtf_link];

                    % constraint that the sum of all gamma coefficients
                    % must be equal to 1
                    A_qtf_sum1 = sparse( ...
                        ones(quality_vert_num, 1), gamma_indices, ...
                        ones(quality_vert_num, 1), 1, var_length);
                    rhs_qtf_sum1 = 1;
                    sense_qtf_sum1 = '=';

                    % constraint that the sum of all psi variables must be
                    % equal to 1
                    A_qtf_sum2 = sparse( ...
                        ones(quality_tri_num, 1), psi_indices, ...
                        ones(quality_tri_num, 1), 1, var_length);
                    rhs_qtf_sum2 = 1;
                    sense_qtf_sum2 = '=';

                    % constraint identifying non-zero gamma's with psi
                    A_qtf_id_cell = cell(quality_vert_num, 1);

                    for vert_id = 1:quality_vert_num
                        % list of triangles that contain this vertex
                        tri_contain_list = find(any( ...
                            quality_testfuncs.Triangles == vert_id, 2));
                        tri_contain_num = length(tri_contain_list);

                        A_qtf_id_cell{vert_id} = sparse( ...
                            ones(tri_contain_num + 1, 1), ...
                            [gamma_indices(vert_id); ...
                            psi_indices(tri_contain_list)], ...
                            [1; -ones(tri_contain_num, 1)], ...
                            1, var_length);
                    end

                    A_qtf_id = vertcat(A_qtf_id_cell{:});
                    rhs_qtf_id = zeros(quality_vert_num, 1);
                    sense_qtf_id = repmat('<', quality_vert_num, 1);

                    % constraint linking the gamma variables and the point
                    % in the quality space
                    A_qtf_link = sparse( ...
                        repelem([1; 2], quality_vert_num + 1, 1), ...
                        [gamma_indices; quality_indices(1); ...
                        gamma_indices; quality_indices(2)], ...
                        [quality_testfuncs.Vertices(:, 1); -1; ...
                        quality_testfuncs.Vertices(:, 2); -1], ...
                        2, var_length);
                    rhs_qtf_link = zeros(2, 1);
                    sense_qtf_link = repmat('=', 2, 1);

                    % sparse matrix used to directly extract the quality
                    % components of the computed (approximate) minimizers;
                    % used to improve numerical stability since the quality
                    % components of the computed (approximate) minimizers
                    % might slightly violate the equality constraints
                    quality_extract_mat = sparse( ...
                        repelem([1; 2], quality_vert_num, 1), ...
                        [gamma_indices; gamma_indices], ...
                        [quality_testfuncs.Vertices(:, 1); ...
                        quality_testfuncs.Vertices(:, 2)], ...
                        2, var_length);

                    % assemble all constraints related to gamma variables
                    % and psi variables
                    A_qtf = [A_qtf_sum1; A_qtf_sum2; ...
                        A_qtf_id; A_qtf_link];
                    rhs_qtf = [rhs_qtf_sum1; rhs_qtf_sum2; ...
                        rhs_qtf_id; rhs_qtf_link];
                    sense_qtf = [sense_qtf_sum1; sense_qtf_sum2; ...
                        sense_qtf_id; sense_qtf_link];

                    % matrix for directly extracting the values of the test
                    % functions for the marginals from the result of the
                    % solver
                    marg_tf_extract_mat = sparse( ...
                        (1:marg_tf_vert_num)', beta_indices, 1, ...
                        marg_tf_vert_num, var_length);

                    % matrix for directly extracting the values of the test
                    % functions on the quality space from the result of the
                    % solver
                    quality_tf_extract_mat = sparse( ...
                        (1:quality_vert_num)', gamma_indices, 1, ...
                        quality_vert_num, var_length);
                else
                    % constraint that the sum of all beta coefficients must
                    % be equal to 1
                    A_mtf_sum = sparse( ...
                        ones(3 * marg_tf_tri_num, 1), ...
                        beta_indices(:), ...
                        ones(3 * marg_tf_tri_num, 1), ...
                        1, var_length);
                    rhs_mtf_sum = 1;
                    sense_mtf_sum = '=';

                    % constraints identifying non-zero beta's with chi's
                    A_mtf_id_cell = cell(mtf_bit_length, 1);

                    for bit_id = 1:mtf_bit_length
                        % list of beta variables that are representing
                        % triangles with that specific bit being 0 or 1
                        mtf_bit0_beta_list = ...
                            reshape(beta_indices( ...
                            mtf_bisect_cell{bit_id, 1}, :)', [], 1);
                        mtf_bit0_beta_num = length(mtf_bit0_beta_list);
                        mtf_bit1_beta_list = ...
                            reshape(beta_indices( ...
                            mtf_bisect_cell{bit_id, 2}, :)', [], 1);
                        mtf_bit1_beta_num = length(mtf_bit1_beta_list);

                        A_mtf_id_cell{bit_id} = sparse( ...
                            [ones(mtf_bit1_beta_num + 1, 1); ...
                            2 * ones(mtf_bit0_beta_num + 1, 1)], ...
                            [mtf_bit1_beta_list; chi_indices(bit_id); ...
                            mtf_bit0_beta_list; chi_indices(bit_id)], ...
                            [ones(mtf_bit1_beta_num, 1); -1; ...
                            ones(mtf_bit0_beta_num, 1); 1], ...
                            2, var_length);
                    end

                    A_mtf_id = vertcat(A_mtf_id_cell{:});
                    rhs_mtf_id = repmat([0; 1], mtf_bit_length, 1);
                    sense_mtf_id = repmat('<', 2 * mtf_bit_length, 1);

                    % constraint linking the beta variables and the point
                    % in the marginal support
                    A_mtf_link = sparse( ...
                        repelem([1; 2], 3 * marg_tf_tri_num + 1, 1), ...
                        [beta_indices(:); input_indices(1); ...
                        beta_indices(:); input_indices(2)], ...
                        [marg_tf.Vertices(marg_tf.Triangles(:), 1); -1; ...
                        marg_tf.Vertices(marg_tf.Triangles(:), 2); -1], ...
                        2, var_length);
                    rhs_mtf_link = [0; 0];
                    sense_mtf_link = repmat('=', 2, 1);

                    % sparse matrix used to directly extract the input
                    % components of the computed (approximate) minimizers;
                    % used to improve numerical stability since the input
                    % components of the computed (approximate) minimizers
                    % might slightly violate the equality constraints
                    input_extract_mat = sparse( ...
                        repelem([1; 2], 3 * marg_tf_tri_num, 1), ...
                        [beta_indices(:); beta_indices(:)], ...
                        [marg_tf.Vertices(marg_tf.Triangles(:), 1); ...
                        marg_tf.Vertices(marg_tf.Triangles(:), 2)], ...
                        2, var_length);

                    % assemble all constraints related to beta variables
                    % and chi variables
                    A_mtf = [A_mtf_sum; A_mtf_id; A_mtf_link];
                    rhs_mtf = [rhs_mtf_sum; rhs_mtf_id; rhs_mtf_link];
                    sense_mtf = [sense_mtf_sum; sense_mtf_id; ...
                        sense_mtf_link];

                    % constraint that the sum of all gamma coefficients
                    % must be equal to 1
                    A_qtf_sum = sparse( ...
                        ones(3 * quality_tri_num, 1), ...
                        gamma_indices(:), ...
                        ones(3 * quality_tri_num, 1), ...
                        1, var_length);
                    rhs_qtf_sum = 1;
                    sense_qtf_sum = '=';

                    % constraints identifying non-zero gamma's with psi's
                    A_qtf_id_cell = cell(qtf_bit_length, 1);

                    for bit_id = 1:qtf_bit_length
                        % list of gamma variables that are representing
                        % triangles with that specific bit being 0 or 1
                        qtf_bit0_gamma_list = ...
                            reshape(gamma_indices( ...
                            qtf_bisect_cell{bit_id, 1}, :)', [], 1);
                        qtf_bit0_gamma_num = length(qtf_bit0_gamma_list);
                        qtf_bit1_gamma_list = ...
                            reshape(gamma_indices( ...
                            qtf_bisect_cell{bit_id, 2}, :)', [], 1);
                        qtf_bit1_gamma_num = length(qtf_bit1_gamma_list);

                        A_qtf_id_cell{bit_id} = sparse( ...
                            [ones(qtf_bit1_gamma_num + 1, 1); ...
                            2 * ones(qtf_bit0_gamma_num + 1, 1)], ...
                            [qtf_bit1_gamma_list; psi_indices(bit_id); ...
                            qtf_bit0_gamma_list; psi_indices(bit_id)], ...
                            [ones(qtf_bit1_gamma_num, 1); -1; ...
                            ones(qtf_bit0_gamma_num, 1); 1], ...
                            2, var_length);
                    end

                    A_qtf_id = vertcat(A_qtf_id_cell{:});
                    rhs_qtf_id = repmat([0; 1], qtf_bit_length, 1);
                    sense_qtf_id = repmat('<', 2 * qtf_bit_length, 1);

                    % constraint linking the gamma variables and the point
                    % in the quality space
                    A_qtf_link = sparse( ...
                        repelem([1; 2], 3 * quality_tri_num + 1, 1), ...
                        [gamma_indices(:); quality_indices(1); ...
                        gamma_indices(:); quality_indices(2)], ...
                        [quality_testfuncs.Vertices( ...
                        quality_testfuncs.Triangles(:), 1); -1; ...
                        quality_testfuncs.Vertices( ...
                        quality_testfuncs.Triangles(:), 2); -1], ...
                        2, var_length);
                    rhs_qtf_link = [0; 0];
                    sense_qtf_link = repmat('=', 2, 1);

                    % sparse matrix used to directly extract the quality
                    % components of the computed (approximate) minimizers;
                    % used to improve numerical stability since the quality
                    % components of the computed (approximate) minimizers
                    % might slightly violate the equality constraints
                    quality_extract_mat = sparse( ...
                        repelem([1; 2], 3 * quality_tri_num, 1), ...
                        [gamma_indices(:); gamma_indices(:)], ...
                        [quality_testfuncs.Vertices( ...
                        quality_testfuncs.Triangles(:), 1); ...
                        quality_testfuncs.Vertices( ...
                        quality_testfuncs.Triangles(:), 2)], ...
                        2, var_length);

                    % assemble all constraints related to gamma variables
                    % and psi variables
                    A_qtf = [A_qtf_sum; A_qtf_id; A_qtf_link];
                    rhs_qtf = [rhs_qtf_sum; rhs_qtf_id; rhs_qtf_link];
                    sense_qtf = [sense_qtf_sum; sense_qtf_id; ...
                        sense_qtf_link];

                    % matrix for directly extracting the values of the test
                    % functions for the marginals from the result of the
                    % solver
                    marg_tf_extract_mat = sparse( ...
                        marg_tf.Triangles(:), ...
                        beta_indices(:), ...
                        1, marg_tf_vert_num, var_length);

                    % matrix for directly extracting the values of the test
                    % functions on the quality space from the result of the
                    % solver
                    quality_tf_extract_mat = sparse( ...
                        quality_testfuncs.Triangles(:), ...
                        gamma_indices(:), ...
                        1, quality_vert_num, var_length);
                end

                % assemble the gurobi model
                gl_model = struct;
                gl_model.modelsense = 'min';
                gl_model.objcon = 0;
                gl_model.obj = zeros(var_length, 1);
                gl_model.obj(xi_index) = cf_weights(marg_id);

                gl_model.lb = [input_lb; quality_lb; cf_lb; tf_lb];
                gl_model.ub = [input_ub; quality_ub; cf_ub; tf_ub];
                gl_model.vtype = [input_vtype; quality_vtype; ...
                    cf_vtype; tf_vtype];

                gl_model.A = [A_cf; A_mtf; A_qtf];
                gl_model.rhs = [rhs_cf; rhs_mtf; rhs_qtf];
                gl_model.sense = [sense_cf; sense_mtf; sense_qtf];

                % store the gurobi model and various indices for
                % inserting inputs and retrieving outputs
                gm = struct;
                gm.GurobiModel = gl_model;
                gm.InputIndices = input_indices;
                gm.InputExtractMat = input_extract_mat;
                gm.QualityIndices = quality_indices;
                gm.QualityExtractMat = quality_extract_mat;
                gm.MargTFCoefsLinkMat = marg_tf_coef_link_mat;
                gm.MargTFCoefsIndices = marg_tf_coef_indices;
                gm.MargTFValuesExtractMat = marg_tf_extract_mat;
                gm.QualityTFCoefsLinkMat = quality_tf_coef_link_mat;
                gm.QualityTFCoefsIndices = quality_tf_coef_indices;
                gm.QualityTFValuesExtractMat = quality_tf_extract_mat;
                obj.Storage.GlobalMin{marg_id} = gm;
            end
        end

        function initializeTransFuncModels(obj)
            % Initialize the gurobi models for solving the minimization
            % problems involved in evaluating the optimized transfer
            % functions

            marg_num = length(obj.Marginals);

            cf_weights = obj.CostFuncs.Weights;
            station_num = obj.CostFuncs.StationNum;
            stations = obj.CostFuncs.TrainStations;
            tr_avail = obj.CostFuncs.TrainAvailabilities;
            mincon_bigM = obj.Storage.MinCost.MILPBigM;

            obj.Storage.TransFunc = cell(marg_num, 1);

            if strcmp(obj.Options.global_formulation, 'CC')
                % build the gurobi models for the global minimization
                % problems; the two-dimensional CPWA functions are modeled 
                % by the convex combination (CC) formulation
                global_form_CC = true;
            elseif strcmp(obj.Options.global_formulation, 'DLOG')
                % build the gurobi models for the global minimization 
                % problems; the two-dimensional CPWA functions are modeled 
                % by the logarithmic disaggregated convex combination
                % (DLOG) formulation
                global_form_CC = false;
            else
                error(['unknown formulation for the global ' ...
                    'minimization problem']);
            end

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_tf = marg.SimplicialTestFuncs;
                marg_tf_vert_num = size(marg_tf.Vertices, 1);
                marg_tf_tri_num = size(marg_tf.Triangles, 1);

                % first compute the indices of the variables
                input_indices = (1:2)';
                input_lb = -inf(2, 1);
                input_ub = inf(2, 1);
                input_vtype = repmat('C', 2, 1);

                % variable that represents the value of the cost
                % function
                xi_index = 3;
                var_counter = 3;
                xi_lb = -inf;
                xi_ub = inf;
                xi_vtype = 'C';

                if tr_avail(marg_id)
                    % variables that represent the L1 distance between
                    % the employee location and the train stations
                    p_indices = reshape((1:2 * station_num)', ...
                        2, station_num)' + var_counter;
                    var_counter = var_counter + 2 * station_num;
                    p_lb = -inf(2 * station_num, 1);
                    p_ub = inf(2 * station_num, 1);
                    p_vtype = repmat('C', 2 * station_num, 1);
                end

                % variables that represent the L1 distance between the
                % employee and the business location
                q_indices = var_counter + (1:2)';
                var_counter = var_counter + 2;
                q_lb = -inf(2, 1);
                q_ub = inf(2, 1);
                q_vtype = repmat('C', 2, 1);

                if tr_avail(marg_id)
                    % variable that represents the residual distance
                    % with respect to the direct distance
                    s_indices0 = var_counter + 1;
                    var_counter = var_counter + 1;

                    % variables that represent the residual distances
                    % with respect to the indirect distances
                    s_indices = var_counter + (1:station_num)';
                    var_counter = var_counter + station_num;

                    s_lb = zeros(station_num + 1, 1);
                    s_ub = inf(station_num + 1, 1);
                    s_vtype = repmat('C', station_num + 1, 1);

                    % binary variable that enforces the minimum
                    % relation between the xi variable and the q
                    % variables
                    iota_indices0 = var_counter + 1;
                    var_counter = var_counter + 1;

                    % binary variables that enforces the minimum
                    % relations between the xi variable, the r
                    % variables, and the p variables
                    iota_indices = var_counter + (1:station_num)';
                    var_counter = var_counter + station_num;

                    iota_lb = -inf(station_num + 1, 1);
                    iota_ub = inf(station_num + 1, 1);
                    iota_vtype = repmat('B', station_num + 1, 1);

                    cf_lb = [xi_lb; p_lb; q_lb; s_lb; iota_lb];
                    cf_ub = [xi_ub; p_ub; q_ub; s_ub; iota_ub];
                    cf_vtype = [xi_vtype; p_vtype; q_vtype; ...
                        s_vtype; iota_vtype];
                else
                    cf_lb = [xi_lb; q_lb];
                    cf_ub = [xi_ub; q_ub];
                    cf_vtype = [xi_vtype; q_vtype];
                end

                if global_form_CC
                    % variables between 0 and 1 that are used to represent 
                    % a point in the marginal support as a convex
                    % combination 
                    beta_indices = var_counter + (1:marg_tf_vert_num)';
                    var_counter = var_counter + marg_tf_vert_num;
                    beta_lb = zeros(marg_tf_vert_num, 1);
                    beta_ub = inf(marg_tf_vert_num, 1);
                    beta_vtype = repmat('C', marg_tf_vert_num, 1);

                    % transformation matrix applied to the vector
                    % containing coefficients of the test functions for the
                    % marginal before inserting the values into the
                    % objective vector of the mixed-integer program
                    marg_tf_coef_link_mat = -speye(marg_tf_vert_num);

                    % indices of locations in the objective vector of the
                    % mixed-integer program where the transformed
                    % coefficients of the test functions for the marginal
                    % will be inserted
                    marg_tf_coef_indices = beta_indices;

                    % binary-valued variables for indicating the triangle 
                    % for the parametric function on the marginal support
                    chi_indices = var_counter + (1:marg_tf_tri_num)';
                    var_counter = var_counter + marg_tf_tri_num;
                    chi_lb = -inf(marg_tf_tri_num, 1);
                    chi_ub = inf(marg_tf_tri_num, 1);
                    chi_vtype = repmat('B', marg_tf_tri_num, 1);

                    tf_lb = [beta_lb; chi_lb];
                    tf_ub = [beta_ub; chi_ub];
                    tf_vtype = [beta_vtype; chi_vtype];
                else
                    mtf_bit_length = ceil(log2(marg_tf_tri_num));
                    mtf_bisect_cell = cell(mtf_bit_length, 2);

                    for bit_id = 1:mtf_bit_length
                        logi0 = bitget((0:marg_tf_tri_num - 1)', ...
                            bit_id) == 0;
                        mtf_bisect_cell{bit_id, 1} = find(logi0);
                        mtf_bisect_cell{bit_id, 2} = find(~logi0);
                    end

                    % variables between 0 and 1 that are used to represent
                    % a point in a triangle in the marginal support as a
                    % convex combination
                    beta_indices = var_counter ...
                        + reshape((1:3 * marg_tf_tri_num)', ...
                        3, marg_tf_tri_num)';
                    var_counter = var_counter + 3 * marg_tf_tri_num;
                    beta_lb = zeros(3 * marg_tf_tri_num, 1);
                    beta_ub = inf(3 * marg_tf_tri_num, 1);
                    beta_vtype = repmat('C', 3 * marg_tf_tri_num, 1);

                    % transformation matrix applied to the vector
                    % containing coefficients of the test functions for the
                    % marginal before inserting the values into the
                    % objective vector of the mixed-integer program
                    marg_tf_coef_link_mat = sparse( ...
                        (1:3 * marg_tf_tri_num)', ...
                        marg_tf.Triangles(:), ...
                        -1, 3 * marg_tf_tri_num, marg_tf_vert_num);

                    % indices of locations in the objective vector of the
                    % mixed-integer program where the transformed
                    % coefficients of the test functions for the marginal
                    % will be inserted
                    marg_tf_coef_indices = beta_indices(:);

                    % binary-valued variables for indicating the triangle
                    % for the parametric transfer function on the marginal
                    % support
                    chi_indices = var_counter + (1:mtf_bit_length)';
                    var_counter = var_counter + mtf_bit_length;
                    chi_lb = -inf(mtf_bit_length, 1);
                    chi_ub = inf(mtf_bit_length, 1);
                    chi_vtype = repmat('B', mtf_bit_length, 1);

                    tf_lb = [beta_lb; chi_lb];
                    tf_ub = [beta_ub; chi_ub];
                    tf_vtype = [beta_vtype; chi_vtype];
                end

                var_length = var_counter;

                % assemble the constraint matrix
                constr_counter = 0;

                if tr_avail(marg_id)
                    % the constraints that relate the employee
                    % location and the p variables
                    A_pcon_cell = cell(station_num, 1);
                    rhs_pcon_cell = cell(station_num, 1);
                    sense_pcon_cell = cell(station_num, 1);

                    for station_id = 1:station_num
                        A_pcon_cell{station_id} = ...
                            sparse(repmat((1:4)', 2, 1), ...
                            [p_indices(station_id, :)'; ...
                            p_indices(station_id, :)'; ...
                            input_indices; input_indices], ...
                            [1; 1; 1; 1; -1; -1; 1; 1], ...
                            4, var_length);
                        rhs_pcon_cell{station_id} = ...
                            [-stations(station_id, :)'; ...
                            stations(station_id, :)'];
                        sense_pcon_cell{station_id} = ...
                            repmat('>', 4, 1);

                        constr_counter = constr_counter + 4;
                    end

                    A_pcon = vertcat(A_pcon_cell{:});
                    rhs_pcon = vertcat(rhs_pcon_cell{:});
                    sense_pcon = vertcat(sense_pcon_cell{:});
                end

                % the constraints that relate the employee location,
                % the business outlet location, and the q variables
                A_qcon = sparse(repmat((1:4)', 2, 1), ...
                    [q_indices; q_indices; ...
                    input_indices; input_indices], ...
                    [1; 1; 1; 1; 1; 1; -1; -1], ...
                    4, var_length);
                rhs_qcon = zeros(4, 1);
                sense_qcon = repmat('>', 4, 1);

                qcon_row_indices = constr_counter + (1:4)';
                constr_counter = constr_counter + 4;
                qcon_link_mat = sparse((1:4)', ...
                    [1; 2; 1; 2], [1; 1; -1; -1], 4, 2);

                if tr_avail(marg_id)
                    % the following part will encode the constraints
                    % that enforce the minimum relation in the cost
                    % function

                    % part 1.1: linking the xi variable, the q
                    % variables, and the s variable for the direct
                    % distance
                    A_cfmin_11 = sparse(ones(4, 1), ...
                        [xi_index; q_indices; s_indices0], ...
                        [1; -1; -1; 1], 1, var_length);
                    rhs_cfmin_11 = 0;
                    sense_cfmin_11 = '=';
                    constr_counter = constr_counter + 1;

                    % part 1.2: linking the xi variable, the r
                    % variables, the p variables, and the s variables
                    % for the indirect distances
                    A_cfmin_12 = sparse( ...
                        repmat((1:station_num)', 4, 1), ...
                        [xi_index * ones(station_num, 1); ...
                        p_indices(:, 1); ...
                        p_indices(:, 2); ...
                        s_indices], ...
                        [ones(station_num, 1); ...
                        -ones(station_num, 1); ...
                        -ones(station_num, 1); ...
                        ones(station_num, 1)], ...
                        station_num, var_length);
                    rhs_cfmin_12 = zeros(station_num, 1);
                    sense_cfmin_12 = repmat('=', station_num, 1);
                    mincon_row_indices = constr_counter ...
                        + (1:station_num)';
                    constr_counter = constr_counter + station_num;

                    % part 2.1: linking the s variable and the iota
                    % variable in the direct distance case via the
                    % big-M constraint
                    A_cfmin_21 = sparse(ones(2, 1), ...
                        [s_indices0; iota_indices0], ...
                        [1; mincon_bigM], ...
                        1, var_length);
                    rhs_cfmin_21 = mincon_bigM;
                    sense_cfmin_21 = '<';
                    constr_counter = constr_counter + 1;

                    % part 2.2: linking the s variables and the iota
                    % variables in the indirect distance case via the
                    % big-M constraints
                    A_cfmin_22 = sparse( ...
                        repmat((1:station_num)', 2, 1), ...
                        [s_indices; iota_indices], ...
                        [ones(station_num, 1); ...
                        mincon_bigM * ones(station_num, 1)], ...
                        station_num, var_length);
                    rhs_cfmin_22 = mincon_bigM * ones(station_num, 1);
                    sense_cfmin_22 = repmat('<', station_num, 1);
                    constr_counter = constr_counter + station_num;

                    % part 3: requiring all iota variables to sum to 1
                    A_cfmin_3 = sparse( ...
                        ones(station_num + 1, 1), ...
                        [iota_indices0; iota_indices], ...
                        ones(station_num + 1, 1), ...
                        1, var_length);
                    rhs_cfmin_3 = 1;
                    sense_cfmin_3 = '=';
                    constr_counter = constr_counter + 1;

                    A_cfmin = [A_cfmin_11; A_cfmin_12; ...
                        A_cfmin_21; A_cfmin_22; ...
                        A_cfmin_3];
                    rhs_cfmin = [rhs_cfmin_11; rhs_cfmin_12; ...
                        rhs_cfmin_21; rhs_cfmin_22; ...
                        rhs_cfmin_3];
                    sense_cfmin = [sense_cfmin_11; sense_cfmin_12; ...
                        sense_cfmin_21; sense_cfmin_22; ...
                        sense_cfmin_3];

                    A_cf = [A_pcon; A_qcon; A_cfmin];
                    rhs_cf = [rhs_pcon; rhs_qcon; rhs_cfmin];
                    sense_cf = [sense_pcon; sense_qcon; sense_cfmin];
                else
                    A_cfmin = sparse(ones(3, 1), ...
                        [xi_index; q_indices], ...
                        [1; -1; -1], ...
                        1, var_length);
                    rhs_cfmin = 0;
                    sense_cfmin = '=';
                    mincon_row_indices = [];
                    constr_counter = constr_counter + 1;

                    A_cf = [A_qcon; A_cfmin];
                    rhs_cf = [rhs_qcon; rhs_cfmin];
                    sense_cf = [sense_qcon; sense_cfmin];
                end

                if global_form_CC
                    % constraint that the sum of all beta coefficients must 
                    % be equal to 1
                    A_mtf_sum1 = sparse( ...
                        ones(marg_tf_vert_num, 1), beta_indices, ...
                        ones(marg_tf_vert_num, 1), 1, var_length);
                    rhs_mtf_sum1 = 1;
                    sense_mtf_sum1 = '=';
                    constr_counter = constr_counter + 1;

                    % constraint that the sum of all chi variables must be
                    % equal to 1
                    A_mtf_sum2 = sparse( ...
                        ones(marg_tf_tri_num, 1), chi_indices, ...
                        ones(marg_tf_tri_num, 1), 1, var_length);
                    rhs_mtf_sum2 = 1;
                    sense_mtf_sum2 = '=';
                    constr_counter = constr_counter + 1;

                    % constraint identifying non-zero beta's with chi
                    A_mtf_id_cell = cell(marg_tf_vert_num, 1);

                    for vert_id = 1:marg_tf_vert_num
                        % list of triangles that contain this vertex
                        tri_contain_list = find(any( ...
                            marg_tf.Triangles == vert_id, 2));
                        tri_contain_num = length(tri_contain_list);

                        A_mtf_id_cell{vert_id} = sparse( ...
                            ones(tri_contain_num + 1, 1), ...
                            [beta_indices(vert_id); ...
                            chi_indices(tri_contain_list)], ...
                            [1; -ones(tri_contain_num, 1)], ...
                            1, var_length);
                    end

                    A_mtf_id = vertcat(A_mtf_id_cell{:});
                    rhs_mtf_id = zeros(marg_tf_vert_num, 1);
                    sense_mtf_id = repmat('<', marg_tf_vert_num, 1);
                    constr_counter = constr_counter + marg_tf_vert_num;

                    % constraint linking the beta variables and the point
                    % in the marginal support
                    A_mtf_link = sparse( ...
                        repelem([1; 2], marg_tf_vert_num + 1, 1), ...
                        [beta_indices; input_indices(1); ...
                        beta_indices; input_indices(2)], ...
                        [marg_tf.Vertices(:, 1); -1; ...
                        marg_tf.Vertices(:, 2); -1], ...
                        2, var_length);
                    rhs_mtf_link = zeros(2, 1);
                    sense_mtf_link = repmat('=', 2, 1);
                    constr_counter = constr_counter + 2; %#ok<NASGU>

                    % assemble all constraints related to beta variables
                    % and chi variables
                    A_mtf = [A_mtf_sum1; A_mtf_sum2; ...
                        A_mtf_id; A_mtf_link];
                    rhs_mtf = [rhs_mtf_sum1; rhs_mtf_sum2; ...
                        rhs_mtf_id; rhs_mtf_link];
                    sense_mtf = [sense_mtf_sum1; sense_mtf_sum2; ...
                        sense_mtf_id; sense_mtf_link];
                else
                    % constraint that the sum of all beta coefficients must
                    % be equal to 1
                    A_mtf_sum = sparse( ...
                        ones(3 * marg_tf_tri_num, 1), ...
                        beta_indices(:), ...
                        ones(3 * marg_tf_tri_num, 1), ...
                        1, var_length);
                    rhs_mtf_sum = 1;
                    sense_mtf_sum = '=';
                    constr_counter = constr_counter + 1;

                    % constraints identifying non-zero beta's with chi's
                    A_mtf_id_cell = cell(mtf_bit_length, 1);

                    for bit_id = 1:mtf_bit_length
                        % list of beta variables that are representing
                        % triangles with that specific bit being 0 or 1
                        mtf_bit0_beta_list = ...
                            reshape(beta_indices( ...
                            mtf_bisect_cell{bit_id, 1}, :)', [], 1);
                        mtf_bit0_beta_num = length(mtf_bit0_beta_list);
                        mtf_bit1_beta_list = ...
                            reshape(beta_indices( ...
                            mtf_bisect_cell{bit_id, 2}, :)', [], 1);
                        mtf_bit1_beta_num = length(mtf_bit1_beta_list);

                        A_mtf_id_cell{bit_id} = sparse( ...
                            [ones(mtf_bit1_beta_num + 1, 1); ...
                            2 * ones(mtf_bit0_beta_num + 1, 1)], ...
                            [mtf_bit1_beta_list; chi_indices(bit_id); ...
                            mtf_bit0_beta_list; chi_indices(bit_id)], ...
                            [ones(mtf_bit1_beta_num, 1); -1; ...
                            ones(mtf_bit0_beta_num, 1); 1], ...
                            2, var_length);
                    end

                    A_mtf_id = vertcat(A_mtf_id_cell{:});
                    rhs_mtf_id = repmat([0; 1], mtf_bit_length, 1);
                    sense_mtf_id = repmat('<', 2 * mtf_bit_length, 1);
                    constr_counter = constr_counter + 2 * mtf_bit_length;

                    % constraint linking the beta variables and the point
                    % in the marginal support
                    A_mtf_link = sparse( ...
                        repelem([1; 2], 3 * marg_tf_tri_num + 1, 1), ...
                        [beta_indices(:); input_indices(1); ...
                        beta_indices(:); input_indices(2)], ...
                        [marg_tf.Vertices(marg_tf.Triangles(:), 1); -1; ...
                        marg_tf.Vertices(marg_tf.Triangles(:), 2); -1], ...
                        2, var_length);
                    rhs_mtf_link = [0; 0];
                    sense_mtf_link = repmat('=', 2, 1);
                    constr_counter = constr_counter + 2; %#ok<NASGU>

                    % assemble all constraints related to beta variables
                    % and chi variables
                    A_mtf = [A_mtf_sum; A_mtf_id; A_mtf_link];
                    rhs_mtf = [rhs_mtf_sum; rhs_mtf_id; rhs_mtf_link];
                    sense_mtf = [sense_mtf_sum; sense_mtf_id; ...
                        sense_mtf_link];
                end

                % assemble the gurobi model
                trans_model = struct;
                trans_model.modelsense = 'min';
                trans_model.objcon = 0;
                trans_model.obj = zeros(var_length, 1);
                trans_model.obj(xi_index) = cf_weights(marg_id);

                trans_model.lb = [input_lb; cf_lb; tf_lb];
                trans_model.ub = [input_ub; cf_ub; tf_ub];
                trans_model.vtype = [input_vtype; cf_vtype; tf_vtype];

                trans_model.A = [A_cf; A_mtf];
                trans_model.rhs = [rhs_cf; rhs_mtf];
                trans_model.sense = [sense_cf; sense_mtf];

                % store the gurobi model and various indices for
                % inserting inputs and retrieving outputs
                transfunc = struct;
                transfunc.GurobiModel = trans_model;
                transfunc.InputIndices = input_indices;
                transfunc.MargTFCoefsLinkMat = marg_tf_coef_link_mat;
                transfunc.MargTFCoefsIndices = marg_tf_coef_indices;
                transfunc.QConstraintsRowIndices = qcon_row_indices;
                transfunc.QConstraintsLinkMat = qcon_link_mat;
                transfunc.MinConstraintsRowIndices = mincon_row_indices;

                obj.Storage.TransFunc{marg_id} = transfunc;
            end
        end

        function initializeBeforeRun(obj)
            % Initialize the algorithm by computing some static quantities

            if ~obj.Storage.SimplicialTestFuncsInitialized
                obj.initializeSimplicialTestFuncs();
            end
        end

        function prepareRuntime(obj)
            % Prepare the runtime environment by initializing some
            % variables
            obj.Runtime = struct;

            marg_num = length(obj.Marginals);

            % the warm-start basis for the constraints are stored in a cell
            % array
            obj.Runtime.cbasis_eq = [];
            obj.Runtime.cbasis_ineq_cell = cell(marg_num, 1);
            obj.Runtime.cbasis = [];
            obj.Runtime.vbasis = [];

            obj.Runtime.CutInputs = cell(marg_num, 1);
            obj.Runtime.CutQualities = cell(marg_num, 1);
            obj.Runtime.CutNumList = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                obj.Runtime.cbasis_ineq_cell{marg_id} = zeros(0, 1);

                % initialize the cuts to be empty
                obj.Runtime.CutInputs{marg_id} = zeros(0, 2);
                obj.Runtime.CutQualities{marg_id} = zeros(0, 2);
            end
        end

        function updateRuntimeAfterLP(obj, result)
            % Update the runtime environment after solving each LP
            % Input:
            %   result: struct produced by gurobi

            % update the warm-start basis used by gurobi
            if isfield(result, 'vbasis') && ~isempty(result.vbasis) ...
                    && isfield(result, 'cbasis') && ~isempty(result.cbasis)
                quality_vert_num = ...
                    size(obj.Quality.SimplicialTestFuncs.Vertices, 1);

                obj.Runtime.vbasis = result.vbasis;
                obj.Runtime.cbasis = result.cbasis;

                % update the equality and inequality parts of the basis
                % separately
                obj.Runtime.cbasis_eq = ...
                    obj.Runtime.cbasis(1:quality_vert_num - 1);
                cbasis_ineq = obj.Runtime.cbasis(quality_vert_num:end);

                marg_num = length(obj.Marginals);
                cut_num_list = obj.Runtime.CutNumList;
                cut_counter = 0;

                for marg_id = 1:marg_num
                    obj.Runtime.cbasis_ineq_cell{marg_id} = ...
                        cbasis_ineq(cut_counter ...
                        + (1:cut_num_list(marg_id)));
                    cut_counter = cut_counter + cut_num_list(marg_id);
                end
            end
        end

        function inside = doCheckIfInsideQualitySpace(obj, pts)
            % Check if the points are inside the quality space
            % Input:
            %   pts: two-column matrix containing the input points, where
            %   each row corresponds to an input
            % Output:
            %   inside: boolean vector indicating whether each of the
            %   points is inside the quality space

            check_model = obj.Quality.GurobiModel;
            lb_list = all(pts >= check_model.lb' ...
                - MT2DBizLoc_ParTrans.INSIDE_TOLERANCE, 2);
            ub_list = all(pts <= check_model.ub' ...
                + MT2DBizLoc_ParTrans.INSIDE_TOLERANCE, 2);

            constr_diff =  pts * check_model.A' - check_model.rhs';

            ineq_list = all(constr_diff(:, check_model.sense ...
                == '<') <= MT2DBizLoc_ParTrans.INSIDE_TOLERANCE, 2);

            inside = lb_list & ub_list & ineq_list;
        end

        function [vals, inside] = doEvaluateQualitySimplicialTestFuncs( ...
                obj, pts)
            % Function that actually evaluates the test functions on the
            % quality space at the given locations
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
            vert_num = size(obj.Quality.SimplicialTestFuncs.Vertices, 1);
            tri_num = size(obj.Quality.SimplicialTestFuncs.Triangles, 1);

            % apply the inverse transformation to recover the coefficients
            % in the affine combination (with respect to the three vertices
            % of each triangle)
            coef_mat = obj.Quality.SimplicialTestFuncs.InvTransMat ...
                * [pts'; ones(1, input_num)];

            % allow a small tolerance to avoid numerical precision issues
            inside = reshape(all(reshape(sparse(coef_mat(:) > -1e-10), ...
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
            triangles = obj.Quality.SimplicialTestFuncs.Triangles';
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

        function model = generateInitialMinModel(obj)
            % Generate the initial linear programming model for gurobi
            % Output:
            %   model: struct containing the linear programming model in
            %   gurobi

            model = struct;
            marg_num = length(obj.Marginals);
            marg_vert_num_list = obj.Storage.MargVertNumList;
            quality_vert_num = ...
                size(obj.Quality.SimplicialTestFuncs.Vertices, 1);
            decivar_num = obj.Storage.DeciVarLength;

            objective_cell = cell(marg_num, 1);
            A_eq_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};

                % the integrals are placed into the objective vector; note
                % that the first test function for the marginal is removed
                % for identification
                marg_integral = marg.SimplicialTestFuncs.Integrals;
                objective_cell{marg_id} = [-1; ...
                    -marg_integral(2:end); ...
                    zeros(quality_vert_num - 1, 1)];

                % matrix for the equality constraints requiring that the
                % transfer functions must sum up to 0; note that the first
                % test function on the quality space is removed for
                % identification
                A_eq_cell{marg_id} = [sparse(quality_vert_num - 1, ...
                    length(marg_integral)), speye(quality_vert_num - 1)];
            end

            % since the cutting plane algorithm assumes that the problem is
            % a minimization problem, we need to transform our maximization
            % problem into a minimization problem
            model.modelsense = 'min';
            model.objcon = 0;

            % the coefficients corresponding to the first test function of
            % each marginal is not included in the decision variables for
            % identification purposes
            model.obj = vertcat(objective_cell{:});

            model.lb = -inf(decivar_num, 1);
            model.ub = inf(decivar_num, 1);

            % store the equality constraints as fields of the model
            model.A_eq = horzcat(A_eq_cell{:});
            model.rhs_eq = zeros(quality_vert_num - 1, 1);

            % store the inequality constraints into cell arrays as fields
            % of the model
            model.A_ineq_cell = cell(marg_num, 1);
            model.rhs_ineq_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                model.A_ineq_cell{marg_id} = sparse(0, ...
                    1 + marg_vert_num_list(marg_id) - 1 ...
                    + quality_vert_num - 1);
                model.rhs_ineq_cell{marg_id} = zeros(0, 1);
            end

            % the constraints in the model need to be assembled from these
            % fields
            model.A = [model.A_eq; blkdiag(model.A_ineq_cell{:})];
            model.rhs = [model.rhs_eq; vertcat(model.rhs_ineq_cell{:})];
            model.sense = [repmat('=', length(model.rhs_eq), 1); ...
                repmat('>', length(model.rhs) - length(model.rhs_eq), 1)];
        end

        function [min_lb, optimizers] = callGlobalMinOracle(obj, vec)
            % Given a decision vector, call the global minimization oracle
            % to approximately determine the "most violated" constraints
            % and return a lower bound for the minimal value
            % Input:
            %   vec: vector corresponding to the current LP solution
            % Outputs:
            %   min_lb: lower bound for the global minimization problem
            %   optimizers: cell array containing structs with fields
            %   inputs and qualites that correspond to the approximate 
            %   optimizers of the global minimization problems

            marg_num = length(obj.Marginals);
            intecept_indices = obj.Storage.DeciVarInterceptIndices;
            marg_tfuncs_indices = obj.Storage.DeciVarMargTestFuncIndices;
            q_tfuncs_indices = obj.Storage.DeciVarQualityTestFuncIndices;

            min_vals_list = zeros(marg_num, 1);
            optimizers = cell(marg_num, 1);

            for marg_id = 1:marg_num
                % disassemble the vector resulted from solving the LP
                % problem; the first test function is added back to the
                % list of test functions for the marginals; similarly the
                % first test function is added back to the list of test
                % functions on the quality space
                objective_const = vec(intecept_indices(marg_id));
                marg_vert_vals = [0; vec(marg_tfuncs_indices{marg_id})];
                quality_vert_vals = [0; vec(q_tfuncs_indices{marg_id})];

                % compute a pool of approximate solutions
                [pool_inputs, pool_qualities, pool_marg_testfuncs, ...
                    pool_quality_testfuncs, pool_vals, pool_costs, ...
                    min_vals_list(marg_id)] ...
                    = obj.computeGlobalMinMILP(marg_id, ...
                    quality_vert_vals, marg_vert_vals, objective_const);
                pool_neg_list = pool_vals < 0;
                optimizers{marg_id} = struct( ...
                    'inputs', pool_inputs(pool_neg_list, :), ...
                    'qualities', pool_qualities(pool_neg_list, :), ...
                    'marg_testfunc_vals', ...
                    pool_marg_testfuncs(pool_neg_list, :), ...
                    'quality_testfunc_vals', ...
                    pool_quality_testfuncs(pool_neg_list, :), ...
                    'costfunc_vals', pool_costs(pool_neg_list), ...
                    'min_val', min_vals_list(marg_id));
            end

            min_lb = sum(min_vals_list);
        end

        function [pool_inputs, pool_qualities, pool_marg_testfuncs, ...
                pool_quality_testfuncs, pool_vals, pool_costs, LB] ...
                = computeGlobalMinMILP(obj, marg_id, ...
                quality_vert_vals, marg_vals, objective_const)
            % Solve the global minimization problem via the mixed-integer
            % linear programming (MILP) formulation
            % Input:
            %   marg_id: the index of the marginal with respect to which
            %   the global minimization problem needs to be solved
            %   quality_vert_vals: values of the simplicial test functions
            %   on the quality spaces at the vertices
            %   marg_vals: values of the simplicial test functions for the 
            %   marginals at the vertices
            %   objective_const: constant part of the objective function
            % Outputs:
            %   pool_inputs: matrix containing the inputs in the support of 
            %   the marginals from the approximate optimizers
            %   pool_qualities: matrix containing the quality part of the
            %   approximate optimizers
            %   pool_marg_testfuncs: sparse matrix containing the values of
            %   the simplicial test functions for the marginal evaluated at
            %   the resulting inputs
            %   pool_quality_testfuncs: sparse matrix containing the values
            %   of the simplicial test functions on the quality space 
            %   evaluated at the resulting quality points
            %   pool_vals: the corresponding objective values of the
            %   approximate optimizers
            %   pool_costs: the corresponding values of the sum of the cost
            %   functions evaluated at the approximate optimizers on the
            %   quality space
            %   LB: the best lower bound found throughout the algorithm

            gm = obj.Storage.GlobalMin{marg_id};
            model = gm.GurobiModel;

            % this vector has the coefficients of all test functions set to
            % zero; thus it can be used to evaluate the sum of the cost
            % functions at the approximate optimizers on the quality space
            costfunc_obj = model.obj;
            costfunc_obj_const = model.objcon;

            % fill in the coefficients in the parametric potential function
            model.obj(gm.MargTFCoefsIndices) = gm.MargTFCoefsLinkMat ...
                * marg_vals;
            model.objcon = model.objcon - objective_const;

            % fill in the coefficients in the parametric transfer function
            model.obj(gm.QualityTFCoefsIndices) = ...
                gm.QualityTFCoefsLinkMat * quality_vert_vals;

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

            pool_inputs = pool_cuts * gm.InputExtractMat';
            pool_qualities = pool_cuts * gm.QualityExtractMat';

            % remove sub-optimal solutions that are duplicates (up to a
            % small tolerance)
            [~, uind] = unique(round([pool_inputs, pool_qualities], 6), ...
                'rows', 'stable');
            pool_inputs = pool_inputs(uind, :);
            pool_qualities = pool_qualities(uind, :);
            pool_cuts = pool_cuts(uind, :);
            pool_vals = pool_vals(uind);

            pool_marg_testfuncs = sparse(pool_cuts ...
                * gm.MargTFValuesExtractMat');
            pool_quality_testfuncs = sparse(pool_cuts ...
                * gm.QualityTFValuesExtractMat');

            pool_costs = pool_cuts * costfunc_obj + costfunc_obj_const;
        end
    
        function updateLSIPUB(obj, min_lb, optimizers) 
            % Update the LSIP upper bound after each call to the global
            % minimization oracle
            % Inputs:
            %   min_lb: the lower bound for the global minimization problem
            %   optimizers: a set of approximate optimizers of the global
            %   minimization problem

            obj.Runtime.LSIP_UB = min(obj.Runtime.LSIP_UB, ...
                obj.Runtime.LSIP_LB - min_lb);

            marg_num = length(obj.Marginals);
            obj.Runtime.GlobalMin = struct;
            obj.Runtime.GlobalMin.MinVals = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                obj.Runtime.GlobalMin.MinVals(marg_id) = ...
                    optimizers{marg_id}.min_val;
            end
        end
        
        function addConstraints(obj, optimizers)
            % Given a collection of approximate optimizers from the global
            % minimization oracle, generate and add the corresponding
            % linear constraints
            % Inputs:
            %   optimizers: output of the method callGlobalOracle

            marg_num = length(obj.Marginals);

            for marg_id = 1:marg_num
                optimizers_i = optimizers{marg_id};
                constr_num = size(optimizers_i.inputs, 1);

                if isfield(optimizers_i, 'marg_testfunc_vals')
                    marg_testfunc_vals = ...
                        optimizers_i.marg_testfunc_vals;
                else
                    marg_testfunc_vals = obj.Marginals{marg_id ...
                        }.evaluateSimplicialTestFuncs(optimizers_i.inputs);
                end

                if isfield(optimizers_i, 'quality_testfunc_vals')
                    quality_testfunc_vals = ...
                        optimizers_i.quality_testfunc_vals;
                else
                    quality_testfunc_vals = ...
                        obj.evaluateQualitySimplicialTestFuncs( ...
                        optimizers_i.qualities);
                end

                % sanitize the output to avoid potential numerical issues;
                % first, check for entries that are very close to 1, set
                % them to 1 and set the remaining entries in those rows to 
                % 0; then, set entries that are close to 0 to 0
                [m_r, m_c, m_v] = find(marg_testfunc_vals);
                closeto1_list = abs(m_v - 1) ...
                    < obj.Options.sanitation_threshold;
                m_v(closeto1_list) = 1;
                keep_list = ~(ismember(m_r, m_r(closeto1_list)) ...
                    & ~closeto1_list) & m_v ...
                    >= obj.Options.sanitation_threshold;
                marg_testfunc_vals = sparse(m_r(keep_list), ...
                    m_c(keep_list), m_v(keep_list), ...
                    size(marg_testfunc_vals, 1), ...
                    size(marg_testfunc_vals, 2));

                [q_r, q_c, q_v] = find(quality_testfunc_vals);
                closeto1_list = abs(q_v - 1) ...
                    < obj.Options.sanitation_threshold;
                q_v(closeto1_list) = 1;
                keep_list = ~(ismember(q_r, q_r(closeto1_list)) ...
                    & ~closeto1_list) & q_v ...
                    >= obj.Options.sanitation_threshold;
                quality_testfunc_vals = sparse(q_r(keep_list), ...
                    q_c(keep_list), q_v(keep_list), ...
                    size(quality_testfunc_vals, 1), ...
                    size(quality_testfunc_vals, 2));

                % remove the first test function for identification, then 
                % prepend a column of 1 and add the values of the
                % simplicial test functions on the quality space on the
                % right
                A_new = [sparse(ones(constr_num, 1)), ...
                    sparse(marg_testfunc_vals(:, 2:end)), ...
                    sparse(quality_testfunc_vals(:, 2:end))];

                rhs_new = optimizers_i.costfunc_vals;
                closeto0_list = abs(rhs_new) ...
                    < obj.Options.sanitation_threshold;
                rhs_new(closeto0_list) = 0;

                % add the newly generated constraints to the end
                obj.Runtime.CurrentLPModel.A_ineq_cell{marg_id} = ...
                    [obj.Runtime.CurrentLPModel.A_ineq_cell{marg_id}; ...
                    -A_new];
                obj.Runtime.CurrentLPModel.rhs_ineq_cell{marg_id} = ...
                    [obj.Runtime.CurrentLPModel.rhs_ineq_cell{marg_id}; ...
                    -rhs_new];

                % add the added indices to the runtime environment
                obj.Runtime.CutInputs{marg_id} = ...
                    [obj.Runtime.CutInputs{marg_id}; ...
                    optimizers_i.inputs];
                obj.Runtime.CutQualities{marg_id} = ...
                    [obj.Runtime.CutQualities{marg_id}; ...
                    optimizers_i.qualities];
                obj.Runtime.CutNumList(marg_id) = ...
                    obj.Runtime.CutNumList(marg_id) + constr_num;
            end

            % the constraints in the model need to be assembled again
            obj.Runtime.CurrentLPModel.A = ...
                [obj.Runtime.CurrentLPModel.A_eq; ...
                blkdiag(obj.Runtime.CurrentLPModel.A_ineq_cell{:})];
            obj.Runtime.CurrentLPModel.rhs = ...
                [obj.Runtime.CurrentLPModel.rhs_eq; ...
                vertcat(obj.Runtime.CurrentLPModel.rhs_ineq_cell{:})];
            obj.Runtime.CurrentLPModel.sense = ...
                [repmat('=', length(obj.Runtime.CurrentLPModel.rhs_eq), ...
                1); repmat('>', length(obj.Runtime.CurrentLPModel.rhs) ...
                - length(obj.Runtime.CurrentLPModel.rhs_eq), 1)];

            if ~isempty(obj.Runtime.vbasis) && ~isempty(obj.Runtime.cbasis)
                obj.Runtime.CurrentLPModel.vbasis = obj.Runtime.vbasis;
                
                % set the warm-start basis of the new constraints to 0
                for marg_id = 1:marg_num
                    optimizers_i = optimizers{marg_id};
                    constr_num = size(optimizers_i.inputs, 1);
                    obj.Runtime.cbasis_ineq_cell{marg_id} = ...
                        [obj.Runtime.cbasis_ineq_cell{marg_id}; ...
                        zeros(constr_num, 1)];
                end

                obj.Runtime.CurrentLPModel.cbasis = ...
                    [obj.Runtime.cbasis_eq; ...
                    vertcat(obj.Runtime.cbasis_ineq_cell{:})];
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

            % if obj.Runtime.NumOfInitialConstraints is not set, it means
            % that this is the first call to obj.addConstraints which
            % generates the initial constraints; this number is stored in
            % the runtime environment
            if ~isfield(obj.Runtime, 'NumOfInitialConstraints') ...
                    || isempty(obj.Runtime.NumOfInitialConstraints)
                obj.Runtime.NumOfInitialConstraints = zeros(marg_num, 1);

                for marg_id = 1:marg_num
                    optimizers_i = optimizers{marg_id};
                    constr_num = length(optimizers_i.inputs);
                    obj.Runtime.NumOfInitialConstraints(marg_id) = ...
                        constr_num;
                end
            end
        end

        function reduceConstraints(obj, result)
            % Remove some of the constraints to speed up the LP solver
            % Input:
            %   result: the output from the gurobi LP solver

            if isinf(obj.Options.reduce.thres) ...
                    || obj.Runtime.iter <= 0 ...
                    || obj.Runtime.iter > obj.Options.reduce.max_iter ...
                    || mod(obj.Runtime.iter, obj.Options.reduce.freq) ~= 0
                return;
            end

            marg_num = length(obj.Marginals);
            quality_vert_num = ...
                size(obj.Quality.SimplicialTestFuncs.Vertices, 1);
            cut_num_list = obj.Runtime.CutNumList;

            cut_slack = result.slack(quality_vert_num:end);

            cut_counter = 0;

            for marg_id = 1:marg_num
                % the list of constraints to be kept (here since the
                % directions of the inequalities are all >=, the
                % slackness is non-positive; the threshold specifies
                % the maximum absolute value of slackness)
                slack = cut_slack(cut_counter ...
                    + (1:cut_num_list(marg_id)));
                keep_list = slack >= -obj.Options.reduce.thres;

                % always keep the initial constraints
                keep_list( ...
                    1:obj.Runtime.NumOfInitialConstraints(marg_id)) ...
                    = true;

                % update all variables
                obj.Runtime.CutInputs{marg_id} ...
                    = obj.Runtime.CutInputs{marg_id}(keep_list, :);
                obj.Runtime.CutQualities{marg_id} ...
                    = obj.Runtime.CutQualities{marg_id}(keep_list, :);
                obj.Runtime.CutNumList(marg_id) = sum(keep_list);
                obj.Runtime.CurrentLPModel.A_ineq_cell{marg_id} ...
                    = obj.Runtime.CurrentLPModel.A_ineq_cell{ ...
                    marg_id}(keep_list, :);
                obj.Runtime.CurrentLPModel.rhs_ineq_cell{marg_id} ...
                    = obj.Runtime.CurrentLPModel.rhs_ineq_cell{ ...
                    marg_id}(keep_list);

                if ~isempty(obj.Runtime.cbasis)
                    obj.Runtime.cbasis_ineq_cell{marg_id} ...
                        = obj.Runtime.cbasis_ineq_cell{marg_id}( ...
                        keep_list);
                end

                cut_counter = cut_counter + cut_num_list(marg_id);
            end

            if ~isempty(obj.Runtime.cbasis)
                obj.Runtime.CurrentLPModel.cbasis = ...
                    [obj.Runtime.cbasis_eq; ...
                    vertcat(obj.Runtime.cbasis_ineq_cell{:})];
            end

            % the constraints in the model need to be assembled again
            obj.Runtime.CurrentLPModel.A = ...
                [obj.Runtime.CurrentLPModel.A_eq; ...
                blkdiag(obj.Runtime.CurrentLPModel.A_ineq_cell{:})];
            obj.Runtime.CurrentLPModel.rhs = ...
                [obj.Runtime.CurrentLPModel.rhs_eq; ...
                vertcat(obj.Runtime.CurrentLPModel.rhs_ineq_cell{:})];
            obj.Runtime.CurrentLPModel.sense = ...
                [repmat('=', length(obj.Runtime.CurrentLPModel.rhs_eq), ...
                1); repmat('>', length(obj.Runtime.CurrentLPModel.rhs) ...
                - length(obj.Runtime.CurrentLPModel.rhs_eq), 1)];
        end

        function primal_sol = buildPrimalSolution(obj, result, ~)
            % Given the output from gurobi and a lower bound for the
            % optimal value of the global minimization oracle, build the
            % corresponding primal solution
            % Inputs:
            %   result: output of the gurobi LP solver
            %   violation: a lower bound for the global minimization
            %   problem
            % Output:
            %   primal_sol: the constructed potential functions on the
            %   support of the input measures as well as the transfer
            %   functions on the quality space

            marg_num = length(obj.Marginals);
            vec = result.x;
            primal_sol = cell(marg_num, 1);

            for marg_id = 1:marg_num
                primal_sol{marg_id} = struct;
                primal_sol{marg_id}.Constant = ...
                    vec(obj.Storage.DeciVarInterceptIndices(marg_id)) ...
                    - obj.Runtime.GlobalMin.MinVals(marg_id);

                % add back the first test function for the marginal
                primal_sol{marg_id}.Coefficients = [0; ...
                    vec(obj.Storage.DeciVarMargTestFuncIndices{marg_id})];

                % similarly, add back the first test function on the
                % quality space
                primal_sol{marg_id}.TransFuncCoefficients = [0; vec( ...
                    obj.Storage.DeciVarQualityTestFuncIndices{marg_id})];
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

            marg_num = length(obj.Marginals);
            quality_vert_num = ...
                    size(obj.Quality.SimplicialTestFuncs.Vertices, 1);
            cut_num_list = obj.Runtime.CutNumList;
            cut_dual = result.pi(quality_vert_num:end);

            dual_sol = cell(marg_num, 1);
            cut_counter = 0;

            for marg_id = 1:marg_num
                marg_cut_dual = cut_dual(cut_counter ...
                    + (1:cut_num_list(marg_id)));
                pos_list = marg_cut_dual > 0;
                dual_sol{marg_id} = struct;
                dual_sol{marg_id}.Probabilities = marg_cut_dual(pos_list);
                dual_sol{marg_id}.InputAtoms = ...
                    obj.Runtime.CutInputs{marg_id}(pos_list, :);
                dual_sol{marg_id}.QualityAtoms = ...
                    obj.Runtime.CutQualities{marg_id}(pos_list, :);

                % normalize the probabilities to resolve small numerical
                % inaccuracies
                dual_sol{marg_id}.Probabilities = ...
                    dual_sol{marg_id}.Probabilities / ...
                    sum(dual_sol{marg_id}.Probabilities);

                cut_counter = cut_counter + cut_num_list(marg_id);
            end
        end

        function [samps, disc_atom_index_samps] = ...
                doRandSampleFromPartialReassembly(obj, ...
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

            % make sure that the couplings between the candidate discrete 
            % measure and the discretized marginals are computed
            if ~isfield(obj.Runtime, 'DiscreteCouplingsComputed') ...
                    || ~obj.Runtime.DiscreteCouplingsComputed
                obj.prepareDiscreteOT();
            end

            % make sure that the couplings between the discretized
            % marginals and the true marginals are computed
            if ~isfield(obj.Runtime, 'OTComputed') ...
                    || ~obj.Runtime.OTComputed
                obj.prepareSemiDiscreteOT();
            end
            
            marg_num = length(obj.Marginals);
            disc_marg_atoms_cell = obj.Runtime.DiscretizedMarginalAtoms;
            disc_couplings_cell = obj.Runtime.DiscreteCouplings;

            cand = obj.Runtime.DiscreteQualityCandidate;
            cand_atom_num = length(cand.Probabilities);

            % generate random indices of the atoms in the candidate 
            % discrete measure according to the probabilities
            disc_atom_index_samps = randsample(rand_stream, ...
                cand_atom_num, samp_num, true, cand.Probabilities);

            samps = struct;
            samps.DiscreteInputs = cell(marg_num, 1);
            samps.ContinuousInputs = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                disc_marg_atoms = disc_marg_atoms_cell{marg_id};

                samp_marg_atom_indices = zeros(samp_num, 1);

                % sample the indices of atoms in the discretized marginal
                for cand_atom_id = 1:cand_atom_num
                    fill_list = disc_atom_index_samps == cand_atom_id;

                    [cond_atoms, ~, cond_probs] = find( ...
                        disc_couplings_cell{marg_id}(cand_atom_id, :)');
                    
                    if length(cond_probs) > 1
                        % normalize the conditional probabilities
                        cond_probs = cond_probs / sum(cond_probs);
    
                        cond_samp_num = sum(fill_list);
                        samp_marg_atom_indices(fill_list) = ...
                            randsample(rand_stream, cond_atoms, ...
                            cond_samp_num, true, cond_probs);
                    else
                        % set all samples to be equal to that atom
                        samp_marg_atom_indices(fill_list) = cond_atoms;
                    end
                end

                % store the coordinates of the samples with discrete
                % marginals
                samps.DiscreteInputs{marg_id} ...
                    = disc_marg_atoms(samp_marg_atom_indices, :);

                if ~ismember(marg_id, marg_to_reassemble)
                    % skip those continuous marginals that are not required
                    % to be sampled
                    continue;
                end

                samps.ContinuousInputs{marg_id} = zeros(samp_num, 2);

                % count the number of samples coupled with each of the
                % atoms in the discretized marginal
                atom_num_list = accumarray(samp_marg_atom_indices, 1, ...
                    [size(disc_marg_atoms, 1), 1]);

                % generate from the conditional distributions
                cont_samp_cell = marg.conditionalRandSample( ...
                    atom_num_list, rand_stream);

                % fill in the coupled samples from the continuous marginals
                for atom_id = 1:length(atom_num_list)
                    samps.ContinuousInputs{marg_id}( ...
                        samp_marg_atom_indices == atom_id, :) ...
                        = cont_samp_cell{atom_id};
                end
            end
        end
    end
end

