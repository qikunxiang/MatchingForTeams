classdef MT2DQuad_ParTrans < MT2DQuad
    % Class for matching for teams problems with two-dimensional marginals,
    % two-dimensional quality space, and quadratic cost functions. The
    % problem is solved via parametrizing the transfer functions on the
    % quality space.

    methods(Access = public)
        function obj = MT2DQuad_ParTrans(marginals, weights, ...
                quality_cell, quality_testfuncs, varargin)
            % Constructor method
            % Inputs:
            %   marginals: cell array containing objects of
            %   ProbMeas2D_ConvexPolytope
            %   weights: vector containing the weights corresponding to the
            %   marginals, the weights will sum up to 1
            %   quality_cell: cell array where each cell contains a
            %   two-column matrix indicating the vertices of a polytope
            %   quality_testfuncs: cell array containing the inputs to the
            %   function obj.setQualitySimplicialTestFuncs

            obj@MT2DQuad(marginals, weights, quality_cell, varargin{:});

            % if the test functions are specified, set the test functions
            if exist('quality_testfuncs', 'var') ...
                    && ~isempty(quality_testfuncs)
                obj.setQualitySimplicialTestFuncs(quality_testfuncs{:});
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

            % set the default options for the global minimization oracle
            if ~isfield(obj.GlobalOptions, 'pool_size') ...
                    || isempty(obj.GlobalOptions.pool_size)
                obj.GlobalOptions.pool_size = 100;
            end

            if ~isfield(obj.GlobalOptions, 'display') ...
                    || isempty(obj.GlobalOptions.display)
                obj.GlobalOptions.display = true;
            end

            if ~isfield(obj.GlobalOptions, 'log_file') ...
                    || isempty(obj.GlobalOptions.log_file)
                obj.GlobalOptions.log_file = '';
            end

            % this flag is used to track if the function
            % obj.initializeSimplicialTestFuncs has been called
            obj.Storage.SimplicialTestFuncsInitialized = false;

            marg_num = length(obj.MarginalWeights);
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

            marg_num = length(obj.MarginalWeights);

            % retrieve the dual solution resulted from the cutting-plane
            % algorithm
            dual_sol = obj.Runtime.DualSolution;
            old_coup_probs_cell = cell(marg_num, 1);
            old_coup_indices_cell = cell(marg_num, 1);
            old_coup_q_atoms_cell = cell(marg_num, 1);
            old_coup_q_probs_cell = cell(marg_num, 1);

            % retrieve the old discretized marginals
            old_atoms_cell = cell(marg_num, 1);
            old_probs_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                quality_atoms = dual_sol{marg_id}.QualityAtoms;
                old_coup_probs_cell{marg_id} = ...
                    dual_sol{marg_id}.Probabilities;

                % store only the unique atoms in the quality space and
                % use their indices in the coupling with the marginal
                [old_coup_q_atoms_cell{marg_id}, ~, uind] ...
                    = unique(quality_atoms, 'rows', 'stable');
                old_coup_q_probs_cell{marg_id} = accumarray(uind, ...
                    old_coup_probs_cell{marg_id});
                old_coup_indices_cell{marg_id} = ...
                    [dual_sol{marg_id}.VertexIndices, uind];

                if marg_id == fixed_index
                    keep_list = old_coup_q_probs_cell{marg_id} >= 1e-12;
                    fixed_cand_q_probs = ...
                        old_coup_q_probs_cell{marg_id}(keep_list);
                    fixed_cand_q_probs = fixed_cand_q_probs ...
                        / sum(fixed_cand_q_probs);
                    fixed_cand_q_atoms = ...
                        old_coup_q_atoms_cell{marg_id}(keep_list, :);
                end

                old_atoms_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices;
                old_probs_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;
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

                % the cost function is the Euclidean distance
                dist_mat = pdist2(old_atoms_cell{marg_id}, ...
                    new_atoms, 'euclidean');

                % compute an optimal coupling between the original discrete
                % marginal and the new discrete marginal discrete optimal
                % transport
                [marg_coup_atom_indices, marg_coup_probs] ...
                    = discrete_OT(old_probs_cell{marg_id}, ...
                    new_probs, dist_mat);

                % compute an optimal coupling between the measure on
                % the quality space and the fixed measure on the
                % quality space; this is to ensure that all measures on
                % the quality space satisfy the constraints with
                % respect to the test functions
                q_dist_mat = pdist2(old_coup_q_atoms_cell{marg_id}, ...
                    fixed_meas_q_atoms, 'euclidean');
                [q_coup_atom_indices, q_coup_probs] = discrete_OT( ...
                    old_coup_q_probs_cell{marg_id}, ...
                    fixed_meas_q_probs, q_dist_mat);

                % perform discrete reassembly to get the new coupling;
                % in this case, we need to reassemble both marginals in
                % the coupling
                [coup_indices, coup_probs] = discrete_reassembly( ...
                    old_coup_indices_cell{marg_id}, ...
                    old_coup_probs_cell{marg_id}, ...
                    {marg_coup_atom_indices; q_coup_atom_indices}, ...
                    {marg_coup_probs; q_coup_probs}, 1:2);
                coup_points = fixed_meas_q_atoms( ...
                    coup_indices(:, 2), :);
                coup_marg_vertices = new_atoms(coup_indices(:, 1), :);
                coup_atom_num = length(coup_probs);
                coup_atom_testfunc_vals = sparse(1:coup_atom_num, ...
                    coup_indices(:, 2), 1, coup_atom_num, ...
                    fixed_meas_q_atom_num);

                % compute the corresponding cost function values
                coup_costfunc_vals = obj.MarginalWeights(marg_id) ...
                    * sum(coup_points .* (coup_points ...
                    - 2 * coup_marg_vertices), 2);

                coup_cell{marg_id} = struct( ...
                    'probabilities', coup_probs, ...
                    'vertex_indices', coup_indices(:, 1), ...
                    'points', coup_points, ...
                    'testfunc_vals', coup_atom_testfunc_vals, ...
                    'costfunc_vals', coup_costfunc_vals);
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

            setSimplicialTestFuncs@MT2DQuad(obj, args_cell);
        end

        function initializeSimplicialTestFuncs(obj)
            % Initialize some quantities related to the simplicial test
            % functions of the marginals and the simplicial test functions
            % on the quality space

            marg_num = length(obj.MarginalWeights);

            quality_testfuncs = obj.Quality.SimplicialTestFuncs;
            quality_vert_num = size(quality_testfuncs.Vertices, 1);
            tri_num = size(quality_testfuncs.Triangles, 1);
            obj.Storage.QualityTriNum = tri_num;

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
                obj.Storage.MargVertNumList(marg_id) = size( ...
                    marg.SimplicialTestFuncs.Vertices, 1);

                % note that the first test function for the marginal is
                % removed for identification purposes
                obj.Storage.DeciVarMargTestFuncIndices{marg_id} = ...
                    ind_counter ...
                    + (1:obj.Storage.MargVertNumList(marg_id) - 1)';
                ind_counter = ind_counter ...
                    + obj.Storage.MargVertNumList(marg_id) - 1;

                % note that the first test function on the quality sapce is
                % removed for identification purposes
                obj.Storage.DeciVarQualityTestFuncIndices{marg_id} = ...
                    ind_counter + (1:quality_vert_num - 1)';
                ind_counter = ind_counter + quality_vert_num - 1;
            end

            obj.Storage.DeciVarLength = ind_counter;
            obj.Storage.TotalMargVertNum ...
                = sum(obj.Storage.MargVertNumList);

            % struct to store information for the global minimization
            % oracle
            obj.Storage.GlobalMin = struct;

            % compute the quantities related to the KKT conditions for the
            % minimization problem within each triangle in the quality
            % space
            obj.Storage.GlobalMin.KKTMat = cell(marg_num, 1);
            obj.Storage.GlobalMin.KKTFixed = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg_testfuncs = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs;
                marg_vert = marg_testfuncs.Vertices;

                KKT_mat_cell = cell(tri_num, 1);
                KKT_rhs_cell = cell(tri_num, 1);

                for tri_id = 1:tri_num
                    tri_vert = quality_testfuncs.Vertices( ...
                        quality_testfuncs.Triangles(tri_id, :), :);

                    KKT_mat_full = [(2 * obj.MarginalWeights(marg_id)) ...
                        * (tri_vert * tri_vert'), -eye(3), ones(3, 1); ...
                        ones(1, 3), zeros(1, 4)];

                    KKT_case1 = KKT_mat_full(:, [1, 2, 3, 7]) \ eye(4);
                    KKT_case2 = KKT_mat_full(:, [1, 2, 6, 7]) \ eye(4);
                    KKT_case3 = KKT_mat_full(:, [1, 3, 5, 7]) \ eye(4);
                    KKT_case4 = KKT_mat_full(:, [2, 3, 4, 7]) \ eye(4);
                    KKT_case5 = KKT_mat_full(:, [1, 5, 6, 7]) \ eye(4);
                    KKT_case6 = KKT_mat_full(:, [2, 4, 6, 7]) \ eye(4);
                    KKT_case7 = KKT_mat_full(:, [3, 4, 5, 7]) \ eye(4);

                    KKT_mat_cell{tri_id} = sparse( ...
                        [KKT_case1(1:3, :); ...
                        KKT_case2(1:3, :); ...
                        KKT_case3(1:3, :); ...
                        KKT_case4(1:3, :); ...
                        zeros(1, 4); ...
                        KKT_case5(2:3, :); ...
                        zeros(1, 4); ...
                        KKT_case6(2:3, :); ...
                        zeros(1, 4); ...
                        KKT_case7(2:3, :)]);
                    KKT_rhs_cell{tri_id} = ...
                        [(2 * obj.MarginalWeights(marg_id)) ...
                        * tri_vert * marg_vert'; ...
                        ones(1, size(marg_vert, 1))];
                end

                obj.Storage.GlobalMin.KKTMat{marg_id} ...
                    = blkdiag(KKT_mat_cell{:});
                obj.Storage.GlobalMin.KKTFixed{marg_id} = ...
                    obj.Storage.GlobalMin.KKTMat{marg_id} ...
                    * vertcat(KKT_rhs_cell{:}) ...
                    + repmat([zeros(12, 1); 1; 0; 0; 1; 0; 0; 1; 0; 0], ...
                    tri_num, 1);
            end

            quality_tri_transp = quality_testfuncs.Triangles';
            obj.Storage.GlobalMin.QualityTestFuncIndices = ...
                quality_tri_transp(:);
            obj.Storage.GlobalMin.KKTRHSIndices = find( ...
                repmat([1; 1; 1; 0], tri_num, 1));

            comb_r = repelem((1:7 * tri_num)', 3, 1);
            comb_c = (1:21 * tri_num)';

            horz_coord_cell = cell(tri_num, 1);
            vert_coord_cell = cell(tri_num, 1);

            row_counter = 0;
            col_counter = 0;

            transfunc_row_cell = cell(tri_num, 1);
            transfunc_col_cell = cell(tri_num, 1);
            transfunc_index_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert_indices = quality_testfuncs.Triangles(tri_id, :)';
                tri_vert = [quality_testfuncs.Vertices( ...
                    tri_vert_indices, :); zeros(1, 2)];
                horz_coord_cell{tri_id} = tri_vert( ...
                    [1, 2, 3, ...
                    1, 2, 4, 1, 3, 4, 2, 3, 4, ...
                    1, 4, 4, 2, 4, 4, 3, 4, 4], 1);
                vert_coord_cell{tri_id} = tri_vert( ...
                    [1, 2, 3, ...
                    1, 2, 4, 1, 3, 4, 2, 3, 4, ...
                    1, 4, 4, 2, 4, 4, 3, 4, 4], 2);

                transfunc_row_cell{tri_id} = row_counter ...
                    + [1; 1; 1; 2; 2; 3; 3; 4; 4; 5; 6; 7];
                transfunc_col_cell{tri_id} = col_counter ...
                    + [1; 2; 3; 4; 5; 7; 8; 10; 11; 13; 16; 19];
                row_counter = row_counter + 7;
                col_counter = col_counter + 21;

                transfunc_index_cell{tri_id} = tri_vert_indices( ...
                    [1; 2; 3; 1; 2; 1; 3; 2; 3; 1; 2; 3]);
            end

            obj.Storage.GlobalMin.KKTCoefSumMat = ...
                sparse(comb_r, comb_c, 1, 7 * tri_num, 21 * tri_num);
            obj.Storage.GlobalMin.ConvCombHorzMat = ...
                sparse(comb_r, comb_c, vertcat(horz_coord_cell{:}), ...
                7 * tri_num, 21 * tri_num);
            obj.Storage.GlobalMin.ConvCombVertMat = ...
                sparse(comb_r, comb_c, vertcat(vert_coord_cell{:}), ...
                7 * tri_num, 21 * tri_num);
            obj.Storage.GlobalMin.ConvCombCoefMat = cell(3, 1);
            obj.Storage.GlobalMin.ConvCombCoefMat{1} = sparse( ...
                comb_r, comb_c, repmat( ...
                [1; 0; 0; 1; 0; 0; 1; 0; 0; 0; 0; 0; ...
                1; 0; 0; 0; 0; 0; 0; 0; 0], tri_num, 1));
            obj.Storage.GlobalMin.ConvCombCoefMat{2} = sparse( ...
                comb_r, comb_c, repmat( ...
                [0; 1; 0; 0; 1; 0; 0; 0; 0; 1; 0; 0; ...
                0; 0; 0; 1; 0; 0; 0; 0; 0], tri_num, 1));
            obj.Storage.GlobalMin.ConvCombCoefMat{3} = sparse( ...
                comb_r, comb_c, repmat( ...
                [0; 0; 1; 0; 0; 0; 0; 1; 0; 0; 1; 0; ...
                0; 0; 0; 0; 0; 0; 1; 0; 0], tri_num, 1));

            obj.Storage.GlobalMin.TransFuncRows = ...
                vertcat(transfunc_row_cell{:});
            obj.Storage.GlobalMin.TransFuncCols = ...
                vertcat(transfunc_col_cell{:});
            obj.Storage.GlobalMin.TransFuncIndices = ...
                vertcat(transfunc_index_cell{:});

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

            marg_num = length(obj.MarginalWeights);
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
                rhs_ineq_cell{marg_id} = obj.MarginalWeights(marg_id) ...
                    * sum(marg_grid_pts .* (marg_grid_pts ...
                    - 2 * quality_grid_pts), 2);

                constr_num_list(marg_id) = marg_vert_num ...
                    * quality_vert_num;
                marg_indices_cell{marg_id} = marg_grid_indices;
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
                coup_marg_indices = ...
                    marg_indices_cell{marg_id}(pos_list, :);
                coup_quality_indices = ...
                    quality_indices_cell{marg_id}(pos_list);
                coup_quality_points = ...
                    quality_points_cell{marg_id}(pos_list, :);
                coup_testfunc_vals = sparse(1:coup_atom_num, ...
                    coup_quality_indices, 1, coup_atom_num, ...
                    quality_vert_num);
                coup_cost = rhs_ineq_cell{marg_id}(pos_list, :);

                coup_cell{marg_id} = struct( ...
                    'probabilities', coup_probs, ...
                    'vertex_indices', coup_marg_indices, ...
                    'points', coup_quality_points, ...
                    'testfunc_vals', coup_testfunc_vals, ...
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
            %   pts: two-column matrix containing the input points
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
            % approximate matching equilibria. The computation is done in
            % batches if necessary.
            % Inputs:
            %   pts: two-column matrix containing the input points
            %   ref_pt: two-element vector indicating a reference point
            %   where the transfer function will evaluate to 0 (default is
            %   the first vertex characterizing the quality space)
            %   batch_size: the maximum number of input points to be
            %   handled at the same time in the vectorized procedure
            %   (default is 10001)
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
                ref_pt = obj.Quality.Vertices(1, :)';
            else
                assert(obj.checkIfInsideQualitySpace(ref_pt'), ...
                    'the reference point is not in the quality space');
            end

            marg_num = length(obj.MarginalWeights);

            % add the reference point
            pts = [ref_pt'; pts];
            pt_num = size(pts, 1);

            vals_mat = zeros(pt_num - 1, marg_num);
            batch_num = ceil(pt_num / batch_size);

            % only evaluate the first (marg_num - 1) transfer functions, 
            % the last one is determined by the balance condition
            for marg_id = 1:marg_num - 1
                marg = obj.Marginals{marg_id};
                marg_weight = obj.MarginalWeights(marg_id);
                marg_vertices = marg.SimplicialTestFuncs.Vertices;
                primal_sol = obj.Runtime.PrimalSolution{marg_id};

                % read the constant and the coefficients of test functions
                % for the marginals
                marg_constant = primal_sol.Constant;
                marg_coefficients = primal_sol.Coefficients;

                vals_cell = cell(batch_num, 1);

                for batch_id = 1:batch_num
                    cur_pts = pts(((batch_id - 1) * batch_size ...
                        + 1):min(batch_id * batch_size, pt_num), :);

                    CPWA_val_mat = (2 * marg_weight) * cur_pts ...
                        * marg_vertices' + marg_coefficients';
                    vals_cell{batch_id} = sum(cur_pts .^ 2, 2) ...
                        .* marg_weight - marg_constant ...
                        - max(CPWA_val_mat, [], 2);

                    if batch_id == 1
                        % store the values of the transfer functions before
                        % shifting evaluated at the reference point
                        ref_vals = vals_cell{batch_id}(1);

                        % shift the transfer functions to make the values
                        % vanish at the reference point, then remove the
                        % reference point
                        vals_cell{batch_id} = ...
                            vals_cell{batch_id}(2:end) - ref_vals;
                    else
                        % shift the transfer functions to make the values
                        % vanish at the reference point
                        vals_cell{batch_id} = vals_cell{batch_id} ...
                            - ref_vals;
                    end
                end

                vals_mat(:, marg_id) = vertcat(vals_cell{:});
            end

            % the last transfer function is obtained from the balance
            % condition, the resulting transfer functions are guaranteed to
            % add up to 0
            vals_mat(:, end) = -sum(vals_mat(:, 1:end - 1), 2);
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
                samp_num, rand_stream, batch_size)
            % Generate independent random sample from the joint
            % distribution consisting of the candidate discrete measure on 
            % the quality space, the corresponding coupling of the 
            % discretized marginals resulted from binding, the coupled 
            % continuous marginals, and the corresponding continuous 
            % measure in the quality space
            % Inputs:
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the
            %   continuous distribution on the quality space (default is
            %   1e4)
            % Output:
            %   samps: struct containing the following fields:
            %       DiscreteIndices: matrix containing the indices of the
            %       atoms in the discretimized marginals; each row
            %       represents a sample and each column represents a
            %       marginal
            %       DiscretePoints: cell array containing the coordinates
            %       of samples from the coupling of the discretized
            %       marginals
            %       DiscreteQualities: two-column matrix containing the
            %       coordinates of samples from the discrete measure on the
            %       quality space
            %       ContinuousPoints: cell array containing the coordinates
            %       of samples from the continuous marginals
            %       ContinuousQualities: two-column matrix containing the
            %       coordinates of samples from the continuous measure on
            %       the quality space

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

            marg_num = length(obj.MarginalWeights);
            samps = obj.doRandSampleFromPartialReassembly( ...
                (1:marg_num)', samp_num, rand_stream);

            % first, compute the weights in the quadratic minimization
            % problem
            weights = zeros(samp_num, 2);

            for marg_id = 1:marg_num
                weights = weights + obj.MarginalWeights(marg_id) ...
                    * samps.ContinuousPoints{marg_id};
            end

            % then, the samples from the continuous measure on the quality
            % space are obtained by solving the quadratic minimization
            % problem over the quality space
            samps.ContinuousQualities = obj.computeQuadMin(weights, ...
                batch_size);
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
            %       of samples from the discrete measure on the
            %       quality space
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

            samps = obj.doRandSampleFromPartialReassembly( ...
                marg_id, samp_num, rand_stream);

            coup = struct;
            coup.AgentTypes = samps.ContinuousPoints{marg_id};

            % get the corresponding samples from the discrete measure on
            % the quality space
            coup.Qualities = samps.DiscreteQualities;
        end

        function ql_samps = randSampleFromOptContinuousQualityDistr( ...
                obj, samp_num, rand_stream, batch_size)
            % Generate independent random sample from the continuous
            % measure on the quality space
            % Inputs:
            %   marg_id: index of the marginal
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the
            %   continuous distribution on the quality space (default is
            %   1e4)
            % Output:
            %   ql_samps: two-column matrix containing the coordinates
            %   of samples from the discrete measure on the quality space

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
                rand_stream, batch_size);
            ql_samps = samps.ContinuousQualities;
        end

        function coup = randSampleFromOptContinuousCoupling(obj, ...
                marg_id, samp_num, rand_stream, batch_size)
            % Generate independent random sample from the coupling between
            % the continuous measure on the quality space and a continuous
            % marginal
            % Inputs:
            %   marg_id: index of the marginal
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the
            %   continuous distribution on the quality space (default is
            %   1e4)
            % Output:
            %   coup: struct containing the following fields:
            %       Qualities: two-column matrix containing the coordinates
            %       of samples from the continuous measure on the
            %       quality space
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

            if ~exist('batch_size', 'var') ...
                    || isempty(batch_size)
                batch_size = 1e4;
            end

            samps = obj.randSampleFromOptJointDistr(samp_num, ...
                rand_stream, batch_size);
            coup = struct;
            coup.Qualities = samps.ContinuousQualities;
            coup.AgentTypes = samps.ContinuousPoints{marg_id};
        end

        function [UB_disc, UB_cont, samps] = getMTUpperBounds(obj, ...
                samp_num, rand_stream, batch_size)
            % Compute two upper bounds for the matching for teams problem
            % based on the discrete measure on the quality space and based
            % on the continuous measure on the quality space. The bounds
            % are approximated by Monte Carlo integration.
            % Inputs:
            %   samp_num: number of samples for Monte Carlo integration
            %   rand_stream: RandStream object (default is
            %   RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the
            %   continuous distribution on the quality space (default is
            %   1e4)
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

            marg_num = length(obj.MarginalWeights);

            samps = obj.randSampleFromOptJointDistr(samp_num, ...
                rand_stream, batch_size);

            UB_disc_mat = zeros(samp_num, marg_num);
            UB_cont_mat = zeros(samp_num, marg_num);

            for marg_id = 1:marg_num
                scaled_marg_samps = 2 * samps.ContinuousPoints{marg_id};
                ql_disc = samps.DiscreteQualities;
                ql_cont = samps.ContinuousQualities;

                % compute each term that contributes to the upper bounds
                % via Monte Carlo integration
                UB_disc_mat(:, marg_id) = sum((ql_disc ...
                    - scaled_marg_samps) .* ql_disc, 2);
                UB_cont_mat(:, marg_id) = sum((ql_cont ...
                    - scaled_marg_samps) .* ql_cont, 2);
            end

            UB_disc_list = UB_disc_mat * obj.MarginalWeights;
            UB_cont_list = UB_cont_mat * obj.MarginalWeights;

            UB_disc = mean(UB_disc_list);
            UB_cont = mean(UB_cont_list);

            samps.UBDiscrete = UB_disc_list;
            samps.UBContinuous = UB_cont_list;
        end

        function EB = getMTErrorBoundBasedOnOT(obj, LSIP_tolerance)
            % Compute the error bound for the objective value of the 
            % matching for teams problem based on the Lipschitz constant of 
            % the cost functions and the optimal transport distances 
            % between the marginals and their discretizations as well as
            % the supremum optimal transport distances between discrete
            % measures on the quality space and the candidate
            % Input:
            %   LSIP_tolerance: tolerance value used in the computation of
            %   the LSIP problem
            % Output:
            %   EB: the computed error bound

            % make sure that the semi-discrete OT problems are solved
            if ~isfield(obj.Storage, 'OTComputed') ...
                    || ~obj.Storage.OTComputed
                obj.performReassembly();
            end

            marg_num = length(obj.MarginalWeights);
            EB = LSIP_tolerance;

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                EB = EB + obj.MarginalWeights(marg_id) ...
                    * (2 * obj.Quality.MaxNorm * marg.OT.Cost);

                if marg_id ~= obj.Options.discmeas_cand_index
                    EB = EB + 2 * obj.MarginalWeights(marg_id) ...
                        * (obj.Quality.MaxNorm + marg.Supp.MaxNorm) ...
                        * 2 * obj.Quality.SimplicialTestFuncs.MeshSize;
                end
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

            marg_num = length(obj.MarginalWeights);

            EB = LSIP_tolerance;

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                EB = EB + 2 * obj.MarginalWeights(marg_id) ...
                    * obj.Quality.MaxNorm ...
                    * 2 * marg.SimplicialTestFuncs.MeshSize;

                if marg_id ~= obj.Options.discmeas_cand_index
                    EB = EB + 2 * obj.MarginalWeights(marg_id) ...
                        * (obj.Quality.MaxNorm + marg.Supp.MaxNorm) ...
                        * 2 * obj.Quality.SimplicialTestFuncs.MeshSize;
                end
            end
        end
    end

    methods(Access = protected)

        function prepareRuntime(obj)
            % Prepare the runtime environment by initializing some
            % variables
            obj.Runtime = struct;

            marg_num = length(obj.MarginalWeights);

            % the warm-start basis for the constraints are stored in a cell
            % array
            obj.Runtime.cbasis_eq = [];
            obj.Runtime.cbasis_ineq_cell = cell(marg_num, 1);
            obj.Runtime.cbasis = [];
            obj.Runtime.vbasis = [];

            obj.Runtime.CutIndices = cell(marg_num, 1);
            obj.Runtime.CutPoints = cell(marg_num, 1);
            obj.Runtime.CutNumList = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                obj.Runtime.cbasis_ineq_cell{marg_id} = zeros(0, 1);

                % initialize the cuts to be empty
                obj.Runtime.CutIndices{marg_id} = zeros(0, 1);
                obj.Runtime.CutPoints{marg_id} = zeros(0, 2);
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

                marg_num = length(obj.MarginalWeights);
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
            marg_num = length(obj.MarginalWeights);
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
                    marg_vert_num_list(marg_id) - 1);
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
            %   vertex_indices and points that correspond to the
            %   approximate optimizers of the global minimization problems

            marg_num = length(obj.MarginalWeights);
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
                [pool_inds, pool_pts, pool_testfuncs, ...
                    pool_vals, pool_costs, min_vals_list(marg_id)] ...
                    = obj.computeGlobalMin(marg_id, quality_vert_vals, ...
                    marg_vert_vals, objective_const);
                pool_neg_list = pool_vals < 0;
                optimizers{marg_id} = struct( ...
                    'vertex_indices', pool_inds(pool_neg_list), ...
                    'points', pool_pts(pool_neg_list, :), ...
                    'testfuncs_vals', pool_testfuncs(pool_neg_list, :), ...
                    'costfunc_vals', pool_costs(pool_neg_list), ...
                    'min_val', min_vals_list(marg_id));
            end

            min_lb = sum(min_vals_list);
        end

        function [pool_inds, pool_pts, pool_testfuncs, ...
                pool_vals, pool_costs, min_val] ...
                = computeGlobalMin(obj, marg_id, quality_vert_vals, ...
                marg_vert_vals, objective_const)
            % Solve the global minimization problem via the KKT conditions
            % Input:
            %   marg_id: the index of the marginal with respect to which
            %   the global minimization problem needs to be solved
            %   quality_vert_vals: values of the simplicial test functions
            %   on the quality spaces at the vertices
            %   marg_vert_vals: values of the simplicial test functions
            %   for the marginals at the vertices
            %   objective_const: constant part of the objective function
            % Outputs:
            %   pool_inds: matrix containing approximate optimizers of the
            %   global minimization problem in the form of vertex indices
            %   in the triangulations of the input measures; the algorithm
            %   will retain up to obj.GlobalOptions.pool_size approximate
            %   optimizers
            %   pool_pts: two-column matrix containing the approximate
            %   optimizers in the quality space
            %   pool_vals: the corresponding objective values of the
            %   approximate optimizers
            %   pool_costs: the corresponding values of the cost function
            %   evaluated at the approximate optimizers
            %   min_val: the computed minimal value

            % open the log file
            if ~isempty(obj.GlobalOptions.log_file)
                log_file = fopen(obj.GlobalOptions.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- Function call starts ---\n');
            end

            marg_weight = obj.MarginalWeights(marg_id);
            marg_vertices = ...
                obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices;
            marg_vert_num = size(marg_vertices, 1);

            qualtity_triangles = obj.Quality.SimplicialTestFuncs.Triangles;
            quality_vert_num = ...
                size(obj.Quality.SimplicialTestFuncs.Vertices, 1);
            tri_num = size(qualtity_triangles, 1);

            % variable part of the right-hand side of the KKT equations
            KKT_rhs = zeros(tri_num * 4, 1);
            KKT_rhs(obj.Storage.GlobalMin.KKTRHSIndices) = ...
                quality_vert_vals( ...
                obj.Storage.GlobalMin.QualityTestFuncIndices);

            % matrix containing the computed coefficients from the KKT
            % equations
            KKT_coef_mat = obj.Storage.GlobalMin.KKTFixed{marg_id} ...
                + obj.Storage.GlobalMin.KKTMat{marg_id} * KKT_rhs;

            % logical matrix indicating those case that are valid
            % (positivity constraints on the variables are satisfied)
            KKT_valid_mat = (obj.Storage.GlobalMin.KKTCoefSumMat ...
                * (KKT_coef_mat < -1e-14)) == 0;

            % horizontal and vertical coordinates of the candidate
            % optimizers (regardless of whether they are valid solutions)
            opt_horz = obj.Storage.GlobalMin.ConvCombHorzMat ...
                * KKT_coef_mat;
            opt_vert = obj.Storage.GlobalMin.ConvCombVertMat ...
                * KKT_coef_mat;

            % coefficients in the convex combinations of the candidate
            % optimizers
            coef1 = obj.Storage.GlobalMin.ConvCombCoefMat{1} ...
                * KKT_coef_mat;
            coef2 = obj.Storage.GlobalMin.ConvCombCoefMat{2} ...
                * KKT_coef_mat;
            coef3 = obj.Storage.GlobalMin.ConvCombCoefMat{3} ...
                * KKT_coef_mat;

            % prepare a sparse matrix containing the values of the transfer
            % functions on the quality space
            transfunc_mat = sparse( ...
                obj.Storage.GlobalMin.TransFuncRows, ...
                obj.Storage.GlobalMin.TransFuncCols, ...
                -quality_vert_vals( ...
                obj.Storage.GlobalMin.TransFuncIndices), ...
                7 * tri_num, 21 * tri_num);

            % compute the objectives at the candidate optimizers (note that
            % some are invalid optimizers)
            costfunc_mat = marg_weight * (opt_horz .^ 2 + opt_vert .^ 2 ...
                - 2 * (opt_horz .* marg_vertices(:, 1)' ...
                + opt_vert .* marg_vertices(:, 2)'));
            objective_mat = transfunc_mat * KKT_coef_mat ...
                + costfunc_mat - marg_vert_vals' - objective_const;

            % check if only a single case is valid for each set of KKT
            % conditions
            if sum(sum(KKT_valid_mat)) == (tri_num * marg_vert_num)
                objective_tri_min = objective_mat(KKT_valid_mat);
                min_costs = costfunc_mat(KKT_valid_mat);

                % compute the minimizers
                min_inds = repelem((1:marg_vert_num)', tri_num, 1);
                min_pts = [opt_horz(KKT_valid_mat), ...
                    opt_vert(KKT_valid_mat)];
                min_coefs = [coef1(KKT_valid_mat), ...
                    coef2(KKT_valid_mat), coef3(KKT_valid_mat)];
            else
                % set the objectives of the invalid optimizers to infinity
                objective_mat(~KKT_valid_mat) = inf;

                % compute the minimum inside each triangle over the 7 cases
                objective_mat_re = reshape(objective_mat(:), ...
                    7, tri_num * marg_vert_num);
                [objective_tri_min, tri_min_case_id] = ...
                    min(objective_mat_re, [], 1);
                objective_tri_min = objective_tri_min';

                % compute the minimizers
                min_inds = repelem((1:marg_vert_num)', tri_num, 1);
                tri_min_lin_id = sub2ind([7 * tri_num, marg_vert_num], ...
                    tri_min_case_id' ...
                    + repmat((0:tri_num - 1)' * 7, marg_vert_num, 1), ...
                    min_inds);
                min_pts = [opt_horz(tri_min_lin_id), ...
                    opt_vert(tri_min_lin_id)];
                min_coefs = [coef1(tri_min_lin_id), ...
                    coef2(tri_min_lin_id), coef3(tri_min_lin_id)];
                min_costs = costfunc_mat(tri_min_lin_id);
            end

            % evaluate the simplicial test functions at all these points
            quantity_triangles_dup = repmat(qualtity_triangles, ...
                marg_vert_num, 1);
            min_testfuncs = sparse(repmat((1:tri_num * marg_vert_num)', ...
                3, 1), quantity_triangles_dup(:), min_coefs(:), ...
                tri_num * marg_vert_num, quality_vert_num);

            % remove (approximately) duplicate optimizers
            [~, unique_inds] = unique(round([min_pts, min_inds], 6), ...
                'rows', 'stable');
            objective_tri_min = objective_tri_min(unique_inds);
            min_costs = min_costs(unique_inds);
            min_pts = min_pts(unique_inds, :);
            min_inds = min_inds(unique_inds);
            min_testfuncs = min_testfuncs(unique_inds, :);

            % sort the optimizers
            [min_vals, sorted_order] = sort(objective_tri_min, 1, ...
                'ascend');
            min_costs = min_costs(sorted_order);
            min_inds = min_inds(sorted_order);
            min_pts = min_pts(sorted_order, :);
            min_testfuncs = min_testfuncs(sorted_order, :);

            % retain up to obj.GlobalOptions.pool_size optimizers
            pool_size = min(obj.GlobalOptions.pool_size, length(min_vals));
            pool_inds = min_inds(1:pool_size);
            pool_pts = min_pts(1:pool_size, :);
            pool_testfuncs = min_testfuncs(1:pool_size, :);
            pool_vals = min_vals(1:pool_size);
            pool_costs = min_costs(1:pool_size);
            min_val = pool_vals(1);

            % some points may be slightly outside the quality space due to 
            % numerical errors; project them back to the quality space
            pool_pts = obj.computeQuadMin(pool_pts);

            if obj.GlobalOptions.display
                fprintf(['Found %3d candidates, ' ...
                    'minimum value = %10.4f\n'], pool_size, min_val);
            end

            if ~isempty(obj.GlobalOptions.log_file)
                fprintf(log_file, ...
                    'Found %3d candidates, minimum value = %10.4f\n', ...
                    pool_size, min_val);
            end

            % close the log file
            if ~isempty(obj.GlobalOptions.log_file)
                fprintf(log_file, '--- Function call ends ---\n\n');
                fclose(log_file);
            end
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

            marg_num = length(obj.MarginalWeights);
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

            marg_num = length(obj.MarginalWeights);

            for marg_id = 1:marg_num
                optimizers_i = optimizers{marg_id};
                constr_num = size(optimizers_i.vertex_indices, 1);

                if isfield(optimizers_i, 'testfuncs_vals')
                    quality_testfunc_vals = optimizers_i.testfuncs_vals;
                else
                    quality_testfunc_vals = ...
                        obj.evaluateQualitySimplicialTestFuncs( ...
                        optimizers_i.points);
                end

                % sanitize the output to avoid potential numerical issues;
                % first, check for entries that are very close to 1, set
                % them to 1 and set the remaining entries in those rows to 
                % 0; then, set entries that are close to 0 to 0
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

                % first generate a matrix containing all test functions 
                % (each row corresponds to an approximate optimizer, each 
                % column corresponds to a test function)
                A_marg_full = sparse((1:constr_num)', ...
                    optimizers_i.vertex_indices, 1, ...
                    constr_num, obj.Storage.MargVertNumList(marg_id));

                % remove the first test function for identification, then 
                % prepend a column of 1 and add the values of the
                % simplicial test functions on the quality space on the
                % right
                A_new = [sparse(ones(constr_num, 1)), ...
                    A_marg_full(:, 2:end), ...
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
                obj.Runtime.CutIndices{marg_id} = ...
                    [obj.Runtime.CutIndices{marg_id}; ...
                    optimizers_i.vertex_indices];
                obj.Runtime.CutPoints{marg_id} = ...
                    [obj.Runtime.CutPoints{marg_id}; ...
                    optimizers_i.points];
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
                    constr_num = size(optimizers_i.vertex_indices, 1);
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
                    constr_num = size(optimizers_i.vertex_indices, 1);
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

            marg_num = length(obj.MarginalWeights);
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
                obj.Runtime.CutIndices{marg_id} ...
                    = obj.Runtime.CutIndices{marg_id}(keep_list);
                obj.Runtime.CutPoints{marg_id} ...
                    = obj.Runtime.CutPoints{marg_id}(keep_list, :);
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

            marg_num = length(obj.MarginalWeights);
            vec = result.x;
            primal_sol = cell(marg_num, 1);

            for marg_id = 1:marg_num
                primal_sol{marg_id} = struct;
                primal_sol{marg_id}.Constant = ...
                    vec(obj.Storage.DeciVarInterceptIndices(marg_id)) ...
                    + obj.Runtime.GlobalMin.MinVals(marg_id);

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

            marg_num = length(obj.MarginalWeights);
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
                dual_sol{marg_id}.VertexIndices = ...
                    obj.Runtime.CutIndices{marg_id}(pos_list);
                dual_sol{marg_id}.QualityAtoms = ...
                    obj.Runtime.CutPoints{marg_id}(pos_list, :);

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
            [cand_atoms, ~, cand_uind] = unique(cand.QualityAtoms, ...
                'rows', 'stable');
            cand_atom_num = size(cand_atoms, 1);
            cand_probs = accumarray(cand_uind, cand.Probabilities, ...
                [cand_atom_num, 1]);

            marg_num = length(obj.MarginalWeights);

            coup_cell = cell(marg_num, 1);
            discreteOT_cost_list = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                marg_atom_num = obj.Storage.MargVertNumList(marg_id);

                if marg_id == cand_id
                    % simply store the coupling between the candidate
                    % measure and the discretized marginal
                    coup_cell{marg_id} = sparse(cand_uind, ...
                        cand.VertexIndices, cand.Probabilities, ...
                        cand_atom_num, marg_atom_num);
                    discreteOT_cost_list(marg_id) = 0;

                    continue;
                end

                dual_sol = obj.Runtime.DualSolution{marg_id};
                [quality_atoms, ~, uind] = ...
                    unique(dual_sol.QualityAtoms, 'rows', 'stable');
                quality_atom_num = size(quality_atoms, 1);

                % construct the sparse matrix representing the coupling
                % between the discrete measure on the quality space and the
                % discretized marginal
                q2m_coup_mat = sparse(uind, dual_sol.VertexIndices, ...
                    dual_sol.Probabilities, quality_atom_num, ...
                    marg_atom_num);

                % compute the optimal coupling between the candidate
                % discrete measure and this discrete measure on the quality
                % space
                dist_mat = pdist2(cand_atoms, quality_atoms, 'euclidean');
                quality_probs = accumarray(uind, ...
                    dual_sol.Probabilities, [quality_atom_num, 1]);
                [coup_atoms, coup_probs, discreteOT_cost_list(marg_id)] ...
                    = discrete_OT(cand_probs, quality_probs, dist_mat);
                
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
            obj.Runtime.DiscreteCouplings = coup_cell;
            obj.Runtime.DiscreteCouplingsComputed = true;
        end

        function samps = doRandSampleFromPartialReassembly(obj, ...
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
            %       DiscreteIndices: matrix containing the indices of the
            %       atoms in the discretimized marginals; each row
            %       represents a sample and each column represents a
            %       marginal
            %       DiscretePoints: cell array containing the coordinates
            %       of samples from the coupling of the discretized
            %       marginals
            %       ContinuousPoints: cell array containing the coordinates
            %       of samples from the continuous marginals (only those
            %       included in marg_to_reassmeble will be sampled)

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            % first make sure that the semi-discrete OT problems are solved
            if ~isfield(obj.Storage, 'OTComputed') ...
                    || ~obj.Storage.OTComputed
                obj.performReassembly();
            end

            % then make sure that the couplings between the candidate
            % discrete measure and the discretized marginals are computed
            if ~isfield(obj.Runtime, 'DiscreteCouplingsComputed') ...
                    || ~obj.Runtime.DiscreteCouplingsComputed
                obj.prepareDiscreteOT();
            end

            marg_num = length(obj.MarginalWeights);
            cand = obj.Runtime.DiscreteQualityCandidate;
            cand_atom_num = length(cand.Probabilities);

            % generate random indices of the atoms in the candidate 
            % discrete measure according to the probabilities
            disc_atom_index_samps = randsample(rand_stream, ...
                cand_atom_num, samp_num, true, cand.Probabilities);

            samps = struct;
            samps.DiscreteQualities ...
                = cand.Atoms(disc_atom_index_samps, :);
            samps.DiscreteIndices = zeros(samp_num, marg_num);
            samps.DiscretePoints = cell(marg_num, 1);
            samps.ContinuousPoints = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_atoms = marg.SimplicialTestFuncs.Vertices;

                % sample the indices of atoms in the discretized marginal
                for cand_atom_id = 1:cand_atom_num
                    fill_list = disc_atom_index_samps == cand_atom_id;

                    [cond_atoms, ~, cond_probs] = find( ...
                        obj.Runtime.DiscreteCouplings{marg_id}( ...
                        cand_atom_id, :)');
                    
                    if length(cond_probs) > 1
                        % normalize the conditional probabilities
                        cond_probs = cond_probs / sum(cond_probs);
    
                        cond_samp_num = sum(fill_list);
                        samps.DiscreteIndices(fill_list, marg_id) = ...
                            randsample(rand_stream, cond_atoms, ...
                            cond_samp_num, true, cond_probs);
                    else
                        % set all samples to be equal to that atom
                        samps.DiscreteIndices(fill_list, marg_id) = ...
                            cond_atoms;
                    end
                end

                % store the coordinates of the samples with discrete
                % marginals
                samps.DiscretePoints{marg_id} ...
                    = marg_atoms(samps.DiscreteIndices(:, marg_id), :);

                if ~ismember(marg_id, marg_to_reassemble)
                    % skip those continuous marginals that are not required
                    % to be sampled
                    continue;
                end

                samps.ContinuousPoints{marg_id} = zeros(samp_num, 2);

                % count the number of samples coupled with each of the
                % atoms in the discretized marginal
                atom_num_list = accumarray( ...
                    samps.DiscreteIndices(:, marg_id), 1, ...
                    [size(marg_atoms, 1), 1]);

                % generate from the conditional distributions
                cont_samp_cell = marg.conditionalRandSample( ...
                    atom_num_list, rand_stream);

                % fill in the coupled samples from the continuous marginals
                for atom_id = 1:length(atom_num_list)
                    samps.ContinuousPoints{marg_id}( ...
                        samps.DiscreteIndices(:, marg_id) ...
                        == atom_id, :) = cont_samp_cell{atom_id};
                end
            end
        end
    end
end

