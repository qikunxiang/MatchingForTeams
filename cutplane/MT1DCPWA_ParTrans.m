classdef MT1DCPWA_ParTrans < MT1DCPWA
    % Class for matching for teams problems with one-dimensional marginals,
    % a quality space which is a two-dimensional polytope, and cost 
    % functions that have the form c_i(x_i, z) = l_i(x_i - s' * z), where 
    % l_i is a continuous piece-wise affine (CPWA) function. The problem is
    % solved via parametrizing the transfer functions on the quality space.

    methods(Access = public)
        function obj = MT1DCPWA_ParTrans(marginals, costfuncs, ...
                quality_vertices, quality_testfuncs, varargin)
            % Constructor method
            % Inputs:
            %   marginals: cell array containing objects of 
            %   ProbMeas1D_Interval
            %   costfuncs: cell array containing structs with fields
            %   weights, knots, and values
            %   quality_vertices: two-column matrix containing the vertices
            %   of the quality space
            %   quality_testfuncs: cell array containing the inputs to the
            %   function obj.setQualitySimplicialTestFuncs

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

            % build the representation of the quality space used by the
            % constructor of MT1DCPWA
            quality_repr = struct;
            quality_repr.dim = 2;
            quality_repr.aux_num = 0;
            quality_repr.ineq_A = hp_w;
            quality_repr.ineq_rhs = hp_b;

            obj@MT1DCPWA(marginals, costfuncs, quality_repr, varargin{:});

            obj.Quality.Vertices = q_vert_cc1;

            % if the test functions are specified, set the test functions
            if exist('quality_testfuncs', 'var') ...
                    && ~isempty(quality_testfuncs)
                obj.setQualitySimplicialTestFuncs(quality_testfuncs{:});
            end

            % set the default option for the formulation of global
            % minimization problem (options are 'LOG_DLOG' and 'INC_CC')
            if ~isfield(obj.Options, 'global_formulation') ...
                    || isempty(obj.Options.global_formulation)
                obj.Options.global_formulation = 'LOG_DLOG';
            end

            % set the default option for the index of the discrete measure
            % candidate on the quality space
            if ~isfield(obj.Options, 'discmeas_cand_index') ...
                    || isempty(obj.Options.discmeas_cand_index)
                obj.Options.discmeas_cand_index = 1;
            end

            % set the default option for sanitizing the entries in the
            % inequality constraints
            if ~isfield(obj.Options, 'sanitation_threshold') ...
                    || isempty(obj.Options.sanitation_threshold)
                obj.Options.sanitation_threshold = 0;
            end

            % store a reference point from the quality space
            obj.Quality.SamplePoint = obj.Quality.Vertices(1, :)';

            % this flag is used to track if the function
            % obj.initializeSimplicialTestFuncs has been called
            obj.Storage.SimplicialTestFuncsInitialized = false;

            marg_num = length(obj.Marginals);
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

            setSimplicialTestFuncs@MT1DCPWA(obj, args_cell);
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
                [old_coup_marg_atoms_cell{marg_id}, ~, uind_marg] ...
                    = unique(round(input_atoms, 6), 'stable');
                old_coup_marg_probs_cell{marg_id} = accumarray( ...
                    uind_marg, old_coup_probs_cell{marg_id});

                % store only the unique atoms in the quality space and
                % use their indices in the coupling with the marginal
                [old_coup_q_atoms_cell{marg_id}, ~, uind_quality] ...
                    = unique(round(quality_atoms, 6), 'rows', 'stable');
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
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Knots;
                new_probs = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;

                % the cost function is the absolute distance
                dist_mat = abs(old_coup_marg_atoms_cell{marg_id} ...
                    - new_atoms');

                % compute an optimal coupling between the original discrete
                % marginal and the new discrete marginal discrete optimal
                % transport
                [marg_coup_atom_indices, marg_coup_probs] ...
                    = discrete_OT(old_coup_marg_probs_cell{marg_id}, ...
                    new_probs, dist_mat);

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
            tri_num = size(quality_testfuncs.Triangles, 1);
            obj.Storage.QualityTriNum = tri_num;

            % store the number of knots in the test functions for each
            % marginal
            obj.Storage.MargKnotNumList = zeros(marg_num, 1);

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
                obj.Storage.MargKnotNumList(marg_id) = length( ...
                    marg.SimplicialTestFuncs.Knots);

                % note that the first test function for the marginal is
                % removed for identification purposes
                obj.Storage.DeciVarMargTestFuncIndices{marg_id} = ...
                    ind_counter ...
                    + (1:obj.Storage.MargKnotNumList(marg_id) - 1)';
                ind_counter = ind_counter ...
                    + obj.Storage.MargKnotNumList(marg_id) - 1;

                % note that the first test function on the quality space is
                % removed for identification purposes
                obj.Storage.DeciVarQualityTestFuncIndices{marg_id} = ...
                    ind_counter + (1:quality_vert_num - 1)';
                ind_counter = ind_counter + quality_vert_num - 1;
            end

            obj.Storage.DeciVarLength = ind_counter;
            obj.Storage.TotalMargKnotNum ...
                = sum(obj.Storage.MargKnotNumList);

            if strcmp(obj.Options.global_formulation, 'INC_CC')
                % build the gurobi models for the global minimization 
                % problems; the one-dimensional CPWA functions are modeled 
                % by the incremental (INC) formulation and the 
                % two-dimensional CPWA functions are modeled by the convex 
                % combination (CC) formulation
                obj.Storage.GlobalMinGurobiModels = cell(marg_num, 1);

                for marg_id = 1:marg_num
                    marg = obj.Marginals{marg_id};
                    cf = obj.CostFuncs{marg_id};
                    cf_knot_num = length(cf.Knots);
                    tf = marg.SimplicialTestFuncs;
                    tf_knot_num = length(tf.Knots);

                    gl_model = struct;

                    % first compute the indices of the variables
                    input_index = 1;
                    quality_indices = (2:3)';
                    var_counter = 3;

                    % variables between 0 and 1 that are used to 
                    % interpolate within the intervals for the cost 
                    % function
                    zeta_indices = var_counter + (1:cf_knot_num - 1)';
                    var_counter = var_counter + cf_knot_num - 1;

                    % binary-valued variables for indicating the interval 
                    % for the cost function
                    iota_indices = var_counter + (1:cf_knot_num - 2)';
                    var_counter = var_counter + cf_knot_num - 2;

                    % variables between 0 and 1 that are used to 
                    % interpolate within the intervals for the parametric 
                    % potential function
                    xi_indices = var_counter + (1:tf_knot_num - 1)';
                    var_counter = var_counter + tf_knot_num - 1;

                    % binary-valued variables for indicating the interval 
                    % for the parametric potential function
                    eta_indices = var_counter + (1:tf_knot_num - 2)';
                    var_counter = var_counter + tf_knot_num - 2;

                    % variables between 0 and 1 that are used to represent 
                    % a point in the quality space as a convex combination
                    beta_indices = var_counter + (1:quality_vert_num)';
                    var_counter = var_counter + quality_vert_num;

                    % binary-valued variables for indicating the triangle 
                    % for the parametric transfer function on the quality 
                    % space
                    chi_indices = var_counter + (1:tri_num)';
                    var_counter = var_counter + tri_num;

                    gl_model.modelsense = 'min';
                    gl_model.objcon = cf.Values(1);
                    gl_model.obj = ...
                        [0; 0; 0; cf.ValueDiff; ...
                        zeros(cf_knot_num - 2, 1); ...
                        zeros(tf_knot_num - 1, 1); ...
                        zeros(tf_knot_num - 2, 1); ...
                        zeros(quality_vert_num, 1); zeros(tri_num, 1)];
                    gl_model.lb = ...
                        [-inf; -inf; -inf; ...
                        zeros(cf_knot_num - 1, 1); ...
                        -inf(cf_knot_num - 2, 1); ...
                        zeros(tf_knot_num - 1, 1); ...
                        -inf(tf_knot_num - 2, 1); ...
                        zeros(quality_vert_num, 1); ...
                        -inf(tri_num, 1)];
                    gl_model.ub = ...
                        [inf; inf; inf; ...
                        1; inf(cf_knot_num - 2, 1); ...
                        inf(cf_knot_num - 2, 1); ...
                        1; inf(tf_knot_num - 2, 1); ...
                        inf(tf_knot_num - 2, 1); ...
                        inf(quality_vert_num, 1); ...
                        inf(tri_num, 1)];
                    gl_model.vtype = ...
                        [repmat('C', 3, 1); ...
                        repmat('C', cf_knot_num - 1, 1); ...
                        repmat('B', cf_knot_num - 2, 1); ...
                        repmat('C', tf_knot_num - 1, 1); ...
                        repmat('B', tf_knot_num - 2, 1); ...
                        repmat('C', quality_vert_num, 1); ...
                        repmat('B', tri_num, 1)];
                    gl_model.input_index = input_index;
                    gl_model.quality_indices = quality_indices;

                    % constraint: iota_{t} - zeta_{t} <= 0 and
                    % constraint: iota_{t} - zeta_{t+1} >= 0
                    A_cf_ord_r = repelem((1:(cf_knot_num - 2) * 2)', 2, 1);
                    A_cf_ord_c_mat = [iota_indices'; ...
                        zeta_indices(1:end - 1)'; ...
                        iota_indices'; ...
                        zeta_indices(2:end)'];
                    A_cf_ord_c = A_cf_ord_c_mat(:);
                    A_cf_ord_v = repmat([1; -1; 1; -1], ...
                        cf_knot_num - 2, 1);
                    A_cf_ord = sparse(A_cf_ord_r, A_cf_ord_c, ...
                        A_cf_ord_v, 2 * (cf_knot_num - 2), var_counter);
                    rhs_cf_ord = zeros(2 * (cf_knot_num - 2), 1);
                    sense_cf_ord = repmat(['<'; '>'], ...
                        cf_knot_num - 2, 1);

                    % constraint linking zeta and the inputs
                    A_cf_link = sparse(ones(cf_knot_num - 1 ...
                        + obj.Quality.Dim + 1, 1), ...
                        [zeta_indices; quality_indices; input_index], ...
                        [cf.KnotDiff; cf.Weights; -1], 1, var_counter);
                    rhs_cf_link = -cf.Knots(1);
                    sense_cf_link = '=';

                    % constraint: eta_{t} - xi_{t} <= 0 and
                    % constraint: eta_{t} - xi_{t+1} >= 0
                    A_tf_ord_r = repelem((1:(tf_knot_num - 2) * 2)', 2, 1);
                    A_tf_ord_c_mat = [eta_indices'; ...
                        xi_indices(1:end - 1)'; ...
                        eta_indices'; ...
                        xi_indices(2:end)'];
                    A_tf_ord_c = A_tf_ord_c_mat(:);
                    A_tf_ord_v = repmat([1; -1; 1; -1], ...
                        tf_knot_num - 2, 1);
                    A_tf_ord = sparse(A_tf_ord_r, ...
                        A_tf_ord_c, ...
                        A_tf_ord_v, ...
                        2 * (tf_knot_num - 2), var_counter);
                    rhs_tf_ord = zeros(2 * (tf_knot_num - 2), 1);
                    sense_tf_ord = repmat(['<'; '>'], tf_knot_num - 2, 1);

                    % constraint linking xi and the inputs
                    A_tf_link = sparse(ones(tf_knot_num - 1 + 1, 1), ...
                        [xi_indices; input_index], ...
                        [tf.KnotDiff; -1], ...
                        1, var_counter);
                    rhs_tf_link = -tf.Knots(1);
                    sense_tf_link = '=';

                    % constraint that the sum of all beta coefficients must
                    % be equal to 1
                    A_q_sum1 = sparse( ...
                        ones(quality_vert_num, 1), beta_indices, ...
                        ones(quality_vert_num, 1), 1, var_counter);
                    rhs_q_sum1 = 1;
                    sense_q_sum1 = '=';

                    % constraint that the sum of all chi variables must be
                    % equal to 1
                    A_q_sum2 = sparse( ...
                        ones(tri_num, 1), chi_indices, ...
                        ones(tri_num, 1), 1, var_counter);
                    rhs_q_sum2 = 1;
                    sense_q_sum2 = '=';

                    % constraint identifying non-zero beta's with chi
                    A_q_id_r_cell = cell(quality_vert_num, 1);
                    A_q_id_c_cell = cell(quality_vert_num, 1);
                    A_q_id_v_cell = cell(quality_vert_num, 1);

                    for bit_id = 1:quality_vert_num
                        % list of triangles that contain this vertex
                        tri_contain_list = find(any( ...
                            quality_testfuncs.Triangles == bit_id, 2));
                        tri_contain_num = length(tri_contain_list);

                        A_q_id_r_cell{bit_id} = bit_id ...
                            * ones(tri_contain_num + 1, 1);
                        A_q_id_c_cell{bit_id} = ...
                            [beta_indices(bit_id); ...
                            chi_indices(tri_contain_list)];
                        A_q_id_v_cell{bit_id} = ...
                            [1; -ones(tri_contain_num, 1)];
                    end

                    A_q_id = sparse(vertcat(A_q_id_r_cell{:}), ...
                        vertcat(A_q_id_c_cell{:}), ...
                        vertcat(A_q_id_v_cell{:}), ...
                        quality_vert_num, var_counter);
                    rhs_q_id = zeros(quality_vert_num, 1);
                    sense_q_id = repmat('<', quality_vert_num, 1);

                    % constraint linking beta and the point in the quality
                    % space
                    A_q_link = sparse( ...
                        repelem([1; 2], quality_vert_num + 1, 1), ...
                        [beta_indices; quality_indices(1); ...
                        beta_indices; quality_indices(2)], ...
                        [quality_testfuncs.Vertices(:, 1); -1; ...
                        quality_testfuncs.Vertices(:, 2); -1], ...
                        2, var_counter);
                    rhs_q_link = [0; 0];
                    sense_q_link = repmat('=', 2, 1);

                    gl_model.A = [A_cf_ord; A_cf_link; ...
                        A_tf_ord; A_tf_link; ...
                        A_q_sum1; A_q_sum2; A_q_id; A_q_link];
                    gl_model.rhs = [rhs_cf_ord; rhs_cf_link; ...
                        rhs_tf_ord; rhs_tf_link; ...
                        rhs_q_sum1; rhs_q_sum2; rhs_q_id; rhs_q_link];
                    gl_model.sense = [sense_cf_ord; sense_cf_link; ...
                        sense_tf_ord; sense_tf_link; ...
                        sense_q_sum1; sense_q_sum2; ...
                        sense_q_id; sense_q_link];

                    % indices in the decision vector to place the 
                    % transformed coefficients of the parametric potential 
                    % functions
                    obj_tfcoef_r = repelem((1:tf_knot_num - 1)', 2, 1);
                    obj_tfcoef_c_mat = [1:tf_knot_num - 1; ...
                        2:tf_knot_num];
                    obj_tfcoef_v = repmat([1; -1], tf_knot_num - 1, 1);
                    gl_model.obj_tfcoef_mat = sparse( ...
                        obj_tfcoef_r, ...
                        obj_tfcoef_c_mat(:), ...
                        obj_tfcoef_v, ...
                        tf_knot_num - 1, tf_knot_num);
                    gl_model.obj_tfcoef_indices = xi_indices;

                    % indices in the decision vector to place the 
                    % coefficients of the parametric transfer functions
                    gl_model.obj_qcoef_mat = -speye(quality_vert_num);
                    gl_model.obj_qcoef_indices = beta_indices;
                    
                    % matrix and intercept vector for directly extracting
                    % the values of the test functions for the marginal
                    % from the results of the solver
                    gl_model.tfcoef_extract_mat = sparse( ...
                        [(1:tf_knot_num - 1)'; (2:tf_knot_num)'], ...
                        repmat(xi_indices, 2, 1), ...
                        [-ones(tf_knot_num - 1, 1); ...
                        ones(tf_knot_num - 1, 1)], ...
                        tf_knot_num, var_counter);
                    gl_model.tfcoef_extract_intercept = ...
                        [1; zeros(tf_knot_num - 1, 1)];

                    % matrix for directly extracting the values of the test
                    % functions on the quality space from the result of the
                    % solver
                    gl_model.qcoef_extract_mat = sparse( ...
                        (1:quality_vert_num)', beta_indices, 1, ...
                        quality_vert_num, var_counter);

                    obj.Storage.GlobalMinGurobiModels{marg_id} = gl_model;
                end
            elseif strcmp(obj.Options.global_formulation, 'LOG_DLOG')
                % build the gurobi models for the global minimization 
                % problems; the one-dimensional CPWA functions are modeled 
                % by the logarithmic convex combination (LOG) formulation 
                % and the two-dimensional CPWA functions are modeled by the
                % logarithmic disaggregated convex combination (DLOG) 
                % formulation
                obj.Storage.GlobalMinGurobiModels = cell(marg_num, 1);

                tri_bit_length = ceil(log2(tri_num));
                tri_bisect_cell = cell(tri_bit_length, 2);

                for bit_id = 1:tri_bit_length
                    logi0 = bitget((0:tri_num - 1)', bit_id) == 0;
                    tri_bisect_cell{bit_id, 1} = find(logi0);
                    tri_bisect_cell{bit_id, 2} = find(~logi0);
                end

                % compute the list containing the vertex indices in the
                % triangles
                tri_vert_mat = quality_testfuncs.Triangles';
                tri_vert_list = tri_vert_mat(:);

                for marg_id = 1:marg_num
                    marg = obj.Marginals{marg_id};
                    cf = obj.CostFuncs{marg_id};
                    cf_knot_num = length(cf.Knots);
                    tf = marg.SimplicialTestFuncs;
                    tf_knot_num = length(tf.Knots);

                    % compute the bisections of the knots in the cost
                    % function
                    cf_bit_length = ceil(log2(cf_knot_num - 1));
                    cf_bisect_cell = ...
                        obj.compute1DBisection(cf_knot_num - 1);

                    % compute the bisections in the knots of the test
                    % functions
                    tf_bit_length = ceil(log2(tf_knot_num - 1));
                    tf_bisect_cell = ...
                        obj.compute1DBisection(tf_knot_num - 1);

                    gl_model = struct;

                    % first compute the indices of the variables
                    input_index = 1;
                    quality_indices = (2:3)';
                    var_counter = 3;

                    % variables between 0 and 1 that are used to 
                    % interpolate within the intervals for the cost 
                    % function
                    zeta_indices = var_counter + (1:cf_knot_num)';
                    var_counter = var_counter + cf_knot_num;

                    % binary-valued variables for indicating the interval 
                    % for the cost function
                    iota_indices = var_counter + (1:cf_bit_length)';
                    var_counter = var_counter + cf_bit_length;

                    % variables between 0 and 1 that are used to 
                    % interpolate within the intervals for the parametric 
                    % potential function
                    xi_indices = var_counter + (1:tf_knot_num)';
                    var_counter = var_counter + tf_knot_num;

                    % binary-valued variables for indicating the interval 
                    % for the parametric potential function
                    eta_indices = var_counter + (1:tf_bit_length)';
                    var_counter = var_counter + tf_bit_length;

                    % variables between 0 and 1 that are used to represent 
                    % a point in a triangle in the quality space as a 
                    % convex combination
                    beta_indices = cell(tri_num, 1);

                    for tri_id = 1:tri_num
                        beta_indices{tri_id} = var_counter + (1:3)';
                        var_counter = var_counter + 3;
                    end

                    % binary-valued variables for indicating the triangle 
                    % for the parametric transfer function on the quality 
                    % space
                    chi_indices = var_counter + (1:tri_bit_length)';
                    var_counter = var_counter + tri_bit_length;

                    gl_model.modelsense = 'min';
                    gl_model.objcon = 0;
                    gl_model.obj = ...
                        [0; 0; 0; cf.Values; ...
                        zeros(cf_bit_length, 1); ...
                        zeros(tf_knot_num, 1); ...
                        zeros(tf_bit_length, 1); ...
                        zeros(tri_num * 3, 1); ...
                        zeros(tri_bit_length, 1)];
                    gl_model.lb = ...
                        [-inf; -inf; -inf; ...
                        zeros(cf_knot_num, 1); ...
                        -inf(cf_bit_length, 1); ...
                        zeros(tf_knot_num, 1); ...
                        -inf(tf_bit_length, 1); ...
                        zeros(tri_num * 3, 1); ...
                        -inf(tri_bit_length, 1)];
                    gl_model.ub = ...
                        [inf; inf; inf; ...
                        inf(cf_knot_num, 1); ...
                        inf(cf_bit_length, 1); ...
                        inf(tf_knot_num, 1); ...
                        inf(tf_bit_length, 1); ...
                        inf(tri_num * 3, 1); ...
                        inf(tri_bit_length, 1)];
                    gl_model.vtype = ...
                        [repmat('C', 3, 1); ...
                        repmat('C', cf_knot_num, 1); ...
                        repmat('B', cf_bit_length, 1); ...
                        repmat('C', tf_knot_num, 1); ...
                        repmat('B', tf_bit_length, 1); ...
                        repmat('C', tri_num * 3, 1); ...
                        repmat('B', tri_bit_length, 1)];
                    gl_model.input_index = input_index;
                    gl_model.quality_indices = quality_indices;
                    gl_model.quality_convcombcoef_indices = ...
                        vertcat(beta_indices{:});

                    % constraint that the sum of all zeta coefficients must
                    % be equal to 1
                    A_cf_sum = sparse( ...
                        ones(cf_knot_num, 1), zeta_indices, ...
                        ones(cf_knot_num, 1), ...
                        1, var_counter);
                    rhs_cf_sum = 1;
                    sense_cf_sum = '=';

                    % constraints identifying non-zero zeta's with iota's
                    A_cf_id_r_cell = cell(cf_bit_length, 1);
                    A_cf_id_c_cell = cell(cf_bit_length, 1);
                    A_cf_id_v_cell = cell(cf_bit_length, 1);

                    for bit_id = 1:cf_bit_length
                        cf_bit0_list = cf_bisect_cell{bit_id, 1};
                        cf_bit0_num = length(cf_bit0_list);
                        cf_bit1_list = cf_bisect_cell{bit_id, 2};
                        cf_bit1_num = length(cf_bit1_list);

                        A_cf_id_r_cell{bit_id} = (bit_id - 1) * 2 ...
                            + [ones(cf_bit0_num + 1, 1); ...
                            2 * ones(cf_bit1_num + 1, 1)];
                        A_cf_id_c_cell{bit_id} = ...
                            [zeta_indices(cf_bit0_list); ...
                            iota_indices(bit_id); ...
                            zeta_indices(cf_bit1_list); ...
                            iota_indices(bit_id)];
                        A_cf_id_v_cell{bit_id} = ...
                            [ones(cf_bit0_num, 1); -1; ...
                            ones(cf_bit1_num, 1); 1];
                    end

                    A_cf_id = sparse(vertcat(A_cf_id_r_cell{:}), ...
                        vertcat(A_cf_id_c_cell{:}), ...
                        vertcat(A_cf_id_v_cell{:}), ...
                        2 * cf_bit_length, var_counter);
                    rhs_cf_id = repmat([0; 1], cf_bit_length, 1);
                    sense_cf_id = repmat('<', 2 * cf_bit_length, 1);

                    % constraint linking zeta and the inputs
                    A_cf_link = sparse(ones(cf_knot_num ...
                        + obj.Quality.Dim + 1, 1), ...
                        [zeta_indices; quality_indices; input_index], ...
                        [cf.Knots; cf.Weights; -1], 1, var_counter);
                    rhs_cf_link = 0;
                    sense_cf_link = '=';

                    % constraint that the sum of all xi coefficients must
                    % be equal to 1
                    A_tf_sum = sparse( ...
                        ones(tf_knot_num, 1), xi_indices, ...
                        ones(tf_knot_num, 1), ...
                        1, var_counter);
                    rhs_tf_sum = 1;
                    sense_tf_sum = '=';

                    % constraints identifying non-zero xi's with eta's
                    A_tf_id_r_cell = cell(tf_bit_length, 1);
                    A_tf_id_c_cell = cell(tf_bit_length, 1);
                    A_tf_id_v_cell = cell(tf_bit_length, 1);

                    for bit_id = 1:tf_bit_length
                        tf_bit0_list = tf_bisect_cell{bit_id, 1};
                        tf_bit0_num = length(tf_bit0_list);
                        tf_bit1_list = tf_bisect_cell{bit_id, 2};
                        tf_bit1_num = length(tf_bit1_list);

                        A_tf_id_r_cell{bit_id} = (bit_id - 1) * 2 ...
                            + [ones(tf_bit0_num + 1, 1); ...
                            2 * ones(tf_bit1_num + 1, 1)];
                        A_tf_id_c_cell{bit_id} = ...
                            [xi_indices(tf_bit0_list); ...
                            eta_indices(bit_id); ...
                            xi_indices(tf_bit1_list); ...
                            eta_indices(bit_id)];
                        A_tf_id_v_cell{bit_id} = ...
                            [ones(tf_bit0_num, 1); -1; ...
                            ones(tf_bit1_num, 1); 1];
                    end

                    A_tf_id = sparse(vertcat(A_tf_id_r_cell{:}), ...
                        vertcat(A_tf_id_c_cell{:}), ...
                        vertcat(A_tf_id_v_cell{:}), ...
                        2 * tf_bit_length, var_counter);
                    rhs_tf_id = repmat([0; 1], tf_bit_length, 1);
                    sense_tf_id = repmat('<', 2 * tf_bit_length, 1);

                    % constraint linking xi and the inputs
                    A_tf_link = sparse(ones(tf_knot_num + 1, 1), ...
                        [xi_indices; input_index], ...
                        [tf.Knots; -1], ...
                        1, var_counter);
                    rhs_tf_link = 0;
                    sense_tf_link = '=';

                    % constraint that the sum of all beta coefficients must 
                    % be equal to 1
                    A_q_sum = sparse( ...
                        ones(tri_num * 3, 1), vertcat(beta_indices{:}), ...
                        ones(tri_num * 3, 1), 1, var_counter);
                    rhs_q_sum = 1;
                    sense_q_sum = '=';

                    % constraints identifying non-zero beta's with chi's
                    A_q_id_r_cell = cell(tri_bit_length, 1);
                    A_q_id_c_cell = cell(tri_bit_length, 1);
                    A_q_id_v_cell = cell(tri_bit_length, 1);

                    for bit_id = 1:tri_bit_length
                        % list of beta variables that are representing
                        % triangles with that specific bit being 0 or 1
                        tri_bit0_beta_list = ...
                            vertcat(beta_indices{ ...
                            tri_bisect_cell{bit_id, 1}});
                        tri_bit0_beta_num = length(tri_bit0_beta_list);
                        tri_bit1_beta_list = ...
                            vertcat(beta_indices{ ...
                            tri_bisect_cell{bit_id, 2}});
                        tri_bit1_beta_num = length(tri_bit1_beta_list);

                        A_q_id_r_cell{bit_id} = (bit_id - 1) * 2 ...
                            + [ones(tri_bit1_beta_num + 1, 1); ...
                            2 * ones(tri_bit0_beta_num + 1, 1)];
                        A_q_id_c_cell{bit_id} = ...
                            [tri_bit1_beta_list; chi_indices(bit_id); ...
                            tri_bit0_beta_list; chi_indices(bit_id)];
                        A_q_id_v_cell{bit_id} = ...
                            [ones(tri_bit1_beta_num, 1); -1; ...
                            ones(tri_bit0_beta_num, 1); 1];
                    end

                    A_q_id = sparse(vertcat(A_q_id_r_cell{:}), ...
                        vertcat(A_q_id_c_cell{:}), ...
                        vertcat(A_q_id_v_cell{:}), ...
                        2 * tri_bit_length, var_counter);
                    rhs_q_id = repmat([0; 1], tri_bit_length, 1);
                    sense_q_id = repmat('<', tri_bit_length * 2, 1);

                    % constraint linking beta and the point in the quality
                    % space
                    A_q_link = sparse( ...
                        repelem([1; 2], tri_num * 3 + 1, 1), ...
                        [vertcat(beta_indices{:}); quality_indices(1); ...
                        vertcat(beta_indices{:}); quality_indices(2)], ...
                        [quality_testfuncs.Vertices(tri_vert_list, 1); -1; ...
                        quality_testfuncs.Vertices(tri_vert_list, 2); -1], ...
                        2, var_counter);
                    rhs_q_link = [0; 0];
                    sense_q_link = repmat('=', 2, 1);

                    gl_model.A = [A_cf_sum; A_cf_id; A_cf_link; ...
                        A_tf_sum; A_tf_id; A_tf_link; ...
                        A_q_sum; A_q_id; A_q_link];
                    gl_model.rhs = [rhs_cf_sum; rhs_cf_id; rhs_cf_link; ...
                        rhs_tf_sum; rhs_tf_id; rhs_tf_link; ...
                        rhs_q_sum; rhs_q_id; rhs_q_link];
                    gl_model.sense = ...
                        [sense_cf_sum; sense_cf_id; sense_cf_link; ...
                        sense_tf_sum; sense_tf_id; sense_tf_link; ...
                        sense_q_sum; sense_q_id; sense_q_link];

                    % indices in the decision vector to place the 
                    % transformed coefficients of the parametric potential 
                    % functions
                    gl_model.obj_tfcoef_mat = -speye(tf_knot_num);
                    gl_model.obj_tfcoef_indices = xi_indices;

                    % indices in the decision vector to place the 
                    % coefficients of the parametric transfer functions
                    gl_model.obj_qcoef_mat = ...
                        sparse((1:tri_num * 3)', tri_vert_list, -1, ...
                        tri_num * 3, quality_vert_num);
                    gl_model.obj_qcoef_indices = vertcat(beta_indices{:});

                    % matrix and intercept vector for directly extracting
                    % the values of the test functions for the marginal
                    % from the results of the solver
                    gl_model.tfcoef_extract_mat = sparse( ...
                        (1:tf_knot_num)', xi_indices, 1, ...
                        tf_knot_num, var_counter);
                    gl_model.tfcoef_extract_intercept = ...
                        zeros(tf_knot_num, 1);

                    % matrix for directly extracting the values of the test
                    % functions on the quality space from the result of the
                    % solver
                    gl_model.qcoef_extract_mat = sparse( ...
                        tri_vert_list, ...
                        vertcat(beta_indices{:}), 1, ...
                        quality_vert_num, var_counter);

                    obj.Storage.GlobalMinGurobiModels{marg_id} = gl_model;
                end
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
            LP_options.FeasibilityTol = 1e-9;

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
                    coup_atom_num, marg_knot_num_list(marg_id));
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

            marg = obj.Marginals{marg_id};
            primal_sol = obj.Runtime.PrimalSolution{marg_id};
            vals = marg.evaluateWeightedSumOfSimplicialTestFuncs( ...
                pts, primal_sol.Coefficients, batch_size) ...
                + primal_sol.Constant;
        end

        function vals_mat = evaluateOptTransferFuncs(obj, pts, ref_pt, ...
                batch_size)
            % Evaluate the transfer functions resulted from optimized
            % parametric functions. These transfer functions is part of
            % approximate matching equilibria. The computation is done
            % in batches if necessary.
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

            assert(all(obj.checkIfInsideQualitySpace(pts)), ...
                'some points are outside the quality space');

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
                    primal_sol = obj.Runtime.PrimalSolution{marg_id};
                    coef = primal_sol.Coefficients;

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

            candidate = obj.Runtime.DualSolution{ ...
                obj.Options.discmeas_cand_index};

            distri = struct;
            distri.Probabilities = candidate.Probabilities;
            distri.Atoms = candidate.QualityAtoms;
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
                marg = obj.Marginals{marg_id};
                dual_sol = obj.Runtime.DualSolution{marg_id};
                costfunc = obj.CostFuncs{marg_id};

                marg_atoms = dual_sol.InputAtoms;
                marg_probs = dual_sol.Probabilities;

                % combine atoms that are very close to each other
                [~, uind, mapping] = unique(round(marg_atoms, 6), ...
                    'stable');
                marg_atoms = marg_atoms(uind);
                marg_probs = accumarray(mapping, marg_probs, ...
                    [length(marg_atoms), 1]);

                marg.setCoupledDiscreteMeasure(marg_atoms, marg_probs);
                EB = EB + costfunc.LipschitzConstQuality ...
                    * obj.Runtime.DiscreteOTCosts(marg_id) ...
                    + costfunc.LipschitzConstMarg ...
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
                costfunc = obj.CostFuncs{marg_id};
                
                EB = EB + costfunc.LipschitzConstQuality ...
                    * 2 * obj.Quality.SimplicialTestFuncs.MeshSize ...
                    + costfunc.LipschitzConstMarg ...
                    * 2 * marg.SimplicialTestFuncs.MeshSize;
            end
        end
    end

    methods(Access = protected)

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
                obj.Runtime.CutInputs{marg_id} = zeros(0, 1);
                obj.Runtime.CutQualities{marg_id} = zeros(0, ...
                    obj.Quality.Dim);
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
            marg_knot_num_list = obj.Storage.MargKnotNumList;
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
                    1 + marg_knot_num_list(marg_id) - 1 ...
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
                    'inputs', pool_inputs(pool_neg_list), ...
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
            %   marginals at the knots
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

            model = obj.Storage.GlobalMinGurobiModels{marg_id};

            % this vector has the coefficients of all test functions set to
            % zero; thus it can be used to evaluate the sum of the cost
            % functions at the approximate optimizers on the quality space
            costfunc_obj = model.obj;
            costfunc_obj_const = model.objcon;

            % fill in the coefficients in the parametric potential function
            model.obj(model.obj_tfcoef_indices) = model.obj_tfcoef_mat ...
                * marg_vals;
            model.objcon = model.objcon - objective_const;

            % fill in the coefficients in the parametric transfer function
            model.obj(model.obj_qcoef_indices) = model.obj_qcoef_mat ...
                * quality_vert_vals;

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

            pool_inputs = pool_cuts(:, model.input_index);
            pool_qualities = pool_cuts(:, model.quality_indices);
            pool_marg_testfuncs = pool_cuts * model.tfcoef_extract_mat' ...
                + model.tfcoef_extract_intercept';
            pool_quality_testfuncs = pool_cuts * model.qcoef_extract_mat';

            % make sure that the inputs are within their respective bounds
            % (they may be slightly out of bound due to small numerical
            % inaccuracies)
            pool_inputs = min(max(pool_inputs, ...
                obj.Storage.MarginalLowerBounds(marg_id)), ...
                obj.Storage.MarginalUpperBounds(marg_id));

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
                constr_num = length(optimizers_i.inputs);

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
                    constr_num = length(optimizers_i.inputs);
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
                    = obj.Runtime.CutInputs{marg_id}(keep_list);
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
                    obj.Runtime.CutInputs{marg_id}(pos_list);
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
            [cand_atoms, ~, cand_uind] = unique( ...
                round(cand.QualityAtoms, 6), 'rows', 'stable');
            cand_atom_num = size(cand_atoms, 1);
            cand_probs = accumarray(cand_uind, cand.Probabilities, ...
                [cand_atom_num, 1]);

            marg_num = length(obj.Marginals);

            marg_atoms_cell = cell(marg_num, 1);
            coup_cell = cell(marg_num, 1);
            discreteOT_cost_list = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                dual_sol = obj.Runtime.DualSolution{marg_id};
                [marg_atoms, ~, m_uind] = unique( ...
                    round(dual_sol.InputAtoms, 6), 'rows', 'stable');
                marg_atoms_cell{marg_id} = marg_atoms;
                marg_atom_num = length(marg_atoms);

                if marg_id == cand_id
                    % simply store the coupling between the candidate
                    % measure and the discretized marginal
                    coup_cell{marg_id} = sparse(cand_uind, ...
                        m_uind, cand.Probabilities, ...
                        cand_atom_num, marg_atom_num);
                    discreteOT_cost_list(marg_id) = 0;

                    continue;
                end

                [quality_atoms, ~, q_uind] = ...
                    unique(dual_sol.QualityAtoms, 'rows', 'stable');
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

            marg_num = length(obj.Marginals);
            disc_marg_atoms_cell = obj.Runtime.DiscretizedMarginalAtoms;
            disc_couplings_cell = obj.Runtime.DiscreteCouplings;

            if ~isfield(obj.Runtime, 'CouplingDone') ...
                || ~obj.Runtime.CouplingDone
                % set the discrete measures that the marginals are coupled 
                % with
                for marg_id = 1:marg_num
                    disc_marg_atoms = disc_marg_atoms_cell{marg_id};
                    disc_marg_probs = ...
                        sum(disc_couplings_cell{marg_id}, 1)';
                        
                    obj.Marginals{marg_id}.setCoupledDiscreteMeasure( ...
                        disc_marg_atoms, disc_marg_probs);
                end

                obj.Runtime.CouplingDone = true;
            end

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
                    = disc_marg_atoms(samp_marg_atom_indices);

                if ~ismember(marg_id, marg_to_reassemble)
                    % skip those continuous marginals that are not required
                    % to be sampled
                    continue;
                end

                samps.ContinuousInputs{marg_id} = zeros(samp_num, 1);

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

