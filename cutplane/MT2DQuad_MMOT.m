classdef MT2DQuad_MMOT < MT2DQuad
    % Class for matching for teams problems with two-dimensional marginals,
    % two-dimensional quality space, and quadratic cost functions. The
    % problem is solved via the multi-marginal optimal transport (MMOT)
    % formulation.

    methods(Access = public)
        function obj = MT2DQuad_MMOT(marginals, weights, quality_cell, ...
                varargin)
            % Constructor method
            % Inputs:
            %   marginals: cell array containing objects of
            %   ProbMeas2D_ConvexPolytope
            %   weights: vector containing the weights corresponding to the
            %   marginals, the weights will sum up to 1
            %   quality_cell: cell array where each cell contains a
            %   two-column matrix indicating the vertices of a polytope

            obj@MT2DQuad(marginals, weights, quality_cell, varargin{:});

            % set the default options for the global minimization oracle
            if ~isfield(obj.GlobalOptions, 'tolerance') ...
                    || isempty(obj.GlobalOptions.tolerance)
                obj.GlobalOptions.tolerance = 1e-4;
            end

            if ~isfield(obj.GlobalOptions, 'batch_size') ...
                    || isempty(obj.GlobalOptions.batch_size)
                obj.GlobalOptions.batch_size = 10300;
            end

            if ~isfield(obj.GlobalOptions, 'pool_size') ...
                    || isempty(obj.GlobalOptions.pool_size)
                obj.GlobalOptions.pool_size = 100;
            end

            if ~isfield(obj.GlobalOptions, 'branching') ...
                    || isempty(obj.GlobalOptions.branching)
                obj.GlobalOptions.branching = [100; 100];
            end

            if ~isfield(obj.GlobalOptions, 'termination_thres') ...
                    || isempty(obj.GlobalOptions.termination_thres)
                obj.GlobalOptions.termination_thres = -0.1;
            end

            if ~isfield(obj.GlobalOptions, 'display') ...
                    || isempty(obj.GlobalOptions.display)
                obj.GlobalOptions.display = true;
            end

            if ~isfield(obj.GlobalOptions, 'log_file') ...
                    || isempty(obj.GlobalOptions.log_file)
                obj.GlobalOptions.log_file = '';
            end

            % information for the global minimization problem
            obj.Storage.GlobalMin = struct;

            % information for the grid searching problem in the global
            % minimization problem
            obj.Storage.GlobalMin.Grid = struct;

            poly_num = length(quality_cell);

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

                % rearrange the vertices in counterclockwise order
                vertices_cc = vertices(ccorder(1:end - 1), :);

                first_vert(poly_id, :) = vertices_cc(1, :);
            end

            % store one vertex from each polytope
            obj.Storage.GlobalMin.Grid.AdditionalPoints = first_vert;

            % set the initial rectangle
            obj.Storage.GlobalMin.Grid.InitRectangle = struct;
            obj.Storage.GlobalMin.Grid.InitRectangle.lb ...
                = min(obj.Quality.Vertices, [], 1)';
            obj.Storage.GlobalMin.Grid.InitRectangle.ub ...
                = max(obj.Quality.Vertices, [], 1)';

            % compute all hyperplanes characterizing the polytopes in the
            % form of inequalities w' * z <= b
            hp_w_cell = cell(poly_num, 1);
            hp_b_cell = cell(poly_num, 1);

            % store the number of hyperplanes that impose upper/lower
            % bounds on vertical/horizontal lines
            grid_v_ub_num_list = zeros(poly_num, 1);
            grid_v_lb_num_list = zeros(poly_num, 1);
            grid_h_ub_num_list = zeros(poly_num, 1);
            grid_h_lb_num_list = zeros(poly_num, 1);

            % store information about the hyperplanes that impose
            % upper/lower bounds on vertical/horizontal lines
            grid_v_ub_cell = cell(poly_num, 2);
            grid_v_lb_cell = cell(poly_num, 2);
            grid_h_ub_cell = cell(poly_num, 2);
            grid_h_lb_cell = cell(poly_num, 2);

            % store lower/upper limit restrictions for the
            % vertical/horizontal lines that are induced by
            % horizontal/vertical edges in the polytopes
            grid_v_lim = zeros(poly_num, 2);
            grid_h_lim = zeros(poly_num, 2);

            for poly_id = 1:poly_num
                % retrieve the vertices in the polytopes and their next
                % counterclockwise neighbor
                vertices_cc1 = obj.Quality.Polytopes{poly_id};
                vertices_cc2 = vertices_cc1([2:end, 1], :);

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

                % hyperplanes that impose upper bounds on vertical lines
                v_ub_list = hp_w_cell{poly_id}(:, 2) > 1e-12;
                v_ub_w1 = hp_w_cell{poly_id}(v_ub_list, 1);
                v_ub_w2 = hp_w_cell{poly_id}(v_ub_list, 2);
                v_ub_b = hp_b_cell{poly_id}(v_ub_list);
                grid_v_ub_num_list(poly_id) = length(v_ub_b);
                grid_v_ub_cell{poly_id, 1} = -v_ub_w1 ./ v_ub_w2;
                grid_v_ub_cell{poly_id, 2} = v_ub_b ./ v_ub_w2;

                % hyperplanes that impose lower bounds on vertical lines
                v_lb_list = hp_w_cell{poly_id}(:, 2) < -1e-12;
                v_lb_w1 = hp_w_cell{poly_id}(v_lb_list, 1);
                v_lb_w2 = hp_w_cell{poly_id}(v_lb_list, 2);
                v_lb_b = hp_b_cell{poly_id}(v_lb_list);
                grid_v_lb_num_list(poly_id) = length(v_lb_b);
                grid_v_lb_cell{poly_id, 1} = -v_lb_w1 ./ v_lb_w2;
                grid_v_lb_cell{poly_id, 2} = v_lb_b ./ v_lb_w2;

                % hyperplanes that are vertical and thus impose limits in
                % the horizontal direction
                v_vv_list = abs(hp_w_cell{poly_id}(:, 2)) <= 1e-12;
                v_vv_w1 = hp_w_cell{poly_id}(v_vv_list, 1);
                v_vv_b = hp_b_cell{poly_id}(v_vv_list);

                % hyperplanes that are facing right and imposing lower
                % limits 
                v_vv_left_list = v_vv_w1 < 0;

                if sum(v_vv_left_list) > 0
                    grid_h_lim(poly_id, 1) = max(v_vv_b(v_vv_left_list) ...
                        ./ v_vv_w1(v_vv_left_list)) - 1e-12;
                else
                    grid_h_lim(poly_id, 1) = -inf;
                end

                % hyperplanes that are facing left and imposing upper
                % limits
                if sum(~v_vv_left_list) > 0
                    grid_h_lim(poly_id, 2) = min(v_vv_b(~v_vv_left_list) ...
                        ./ v_vv_w1(~v_vv_left_list)) + 1e-12;
                else
                    grid_h_lim(poly_id, 2) = inf;
                end

                % hyperplanes that impose upper bounds on horizontal lines
                h_ub_list = hp_w_cell{poly_id}(:, 1) > 1e-12;
                h_ub_w1 = hp_w_cell{poly_id}(h_ub_list, 1);
                h_ub_w2 = hp_w_cell{poly_id}(h_ub_list, 2);
                h_ub_b = hp_b_cell{poly_id}(h_ub_list);
                grid_h_ub_num_list(poly_id) = length(h_ub_b);
                grid_h_ub_cell{poly_id, 1} = -h_ub_w2 ./ h_ub_w1;
                grid_h_ub_cell{poly_id, 2} = h_ub_b ./ h_ub_w1;

                % hyperplanes that impose lower bounds on horizontal lines
                h_lb_list = hp_w_cell{poly_id}(:, 1) < -1e-12;
                h_lb_w1 = hp_w_cell{poly_id}(h_lb_list, 1);
                h_lb_w2 = hp_w_cell{poly_id}(h_lb_list, 2);
                h_lb_b = hp_b_cell{poly_id}(h_lb_list);
                grid_h_lb_num_list(poly_id) = length(h_lb_b);
                grid_h_lb_cell{poly_id, 1} = -h_lb_w2 ./ h_lb_w1;
                grid_h_lb_cell{poly_id, 2} = h_lb_b ./ h_lb_w1;

                % hyperplanes that are horizontal and thus impose limits in
                % the vertical direction
                h_hh_list = abs(hp_w_cell{poly_id}(:, 1)) <= 1e-12;
                h_hh_w2 = hp_w_cell{poly_id}(h_hh_list, 2);
                h_hh_b = hp_b_cell{poly_id}(h_hh_list);

                % hyperplanes that are facing up and imposing lower limits 
                h_hh_bot_list = h_hh_w2 < 0;

                if sum(h_hh_bot_list) > 0
                    grid_v_lim(poly_id, 1) = max(h_hh_b(h_hh_bot_list) ...
                        ./ h_hh_w2(h_hh_bot_list)) - 1e-12;
                else
                    grid_v_lim(poly_id, 1) = -inf;
                end

                % hyperplanes that are facing down and imposing upper
                % limits
                if sum(~h_hh_bot_list) > 0
                    grid_v_lim(poly_id, 2) = min(h_hh_b(~h_hh_bot_list) ...
                        ./ h_hh_w2(~h_hh_bot_list)) + 1e-12;
                else
                    grid_v_lim(poly_id, 2) = inf;
                end
            end

            % calculate the maximum number of upper/lower bounding
            % constraints for each vertical/horizontal line
            grid_v_ub_num_max = max(grid_v_ub_num_list);
            grid_v_lb_num_max = max(grid_v_lb_num_list);
            grid_h_ub_num_max = max(grid_h_ub_num_list);
            grid_h_lb_num_max = max(grid_h_lb_num_list);

            for poly_id = 1:poly_num
                % add redundant constraints
                grid_v_ub_cell{poly_id, 1} = [grid_v_ub_cell{poly_id, 1};
                    zeros(grid_v_ub_num_max ...
                    - grid_v_ub_num_list(poly_id), 1)];
                grid_v_ub_cell{poly_id, 2} = [grid_v_ub_cell{poly_id, 2};
                    inf(grid_v_ub_num_max ...
                    - grid_v_ub_num_list(poly_id), 1)];

                grid_v_lb_cell{poly_id, 1} = [grid_v_lb_cell{poly_id, 1};
                    zeros(grid_v_lb_num_max ...
                    - grid_v_lb_num_list(poly_id), 1)];
                grid_v_lb_cell{poly_id, 2} = [grid_v_lb_cell{poly_id, 2};
                    -inf(grid_v_lb_num_max ...
                    - grid_v_lb_num_list(poly_id), 1)];

                grid_h_ub_cell{poly_id, 1} = [grid_h_ub_cell{poly_id, 1};
                    zeros(grid_h_ub_num_max ...
                    - grid_h_ub_num_list(poly_id), 1)];
                grid_h_ub_cell{poly_id, 2} = [grid_h_ub_cell{poly_id, 2};
                    inf(grid_h_ub_num_max ...
                    - grid_h_ub_num_list(poly_id), 1)];

                grid_h_lb_cell{poly_id, 1} = [grid_h_lb_cell{poly_id, 1};
                    zeros(grid_h_lb_num_max ...
                    - grid_h_lb_num_list(poly_id), 1)];
                grid_h_lb_cell{poly_id, 2} = [grid_h_lb_cell{poly_id, 2};
                    -inf(grid_h_lb_num_max ...
                    - grid_h_lb_num_list(poly_id), 1)];
            end

            % store the coefficients, intercepts, and number of constraints
            % of the four types
            obj.Storage.GlobalMin.Grid.VertUB = struct( ...
                'coef', vertcat(grid_v_ub_cell{:, 1}), ...
                'intercept', vertcat(grid_v_ub_cell{:, 2}), ...
                'num', grid_v_ub_num_max);
            obj.Storage.GlobalMin.Grid.VertLB = struct( ...
                'coef', vertcat(grid_v_lb_cell{:, 1}), ...
                'intercept', vertcat(grid_v_lb_cell{:, 2}), ...
                'num', grid_v_lb_num_max);
            obj.Storage.GlobalMin.Grid.HorzUB = struct( ...
                'coef', vertcat(grid_h_ub_cell{:, 1}), ...
                'intercept', vertcat(grid_h_ub_cell{:, 2}), ...
                'num', grid_h_ub_num_max);
            obj.Storage.GlobalMin.Grid.HorzLB = struct( ...
                'coef', vertcat(grid_h_lb_cell{:, 1}), ...
                'intercept', vertcat(grid_h_lb_cell{:, 2}), ...
                'num', grid_h_lb_num_max);

            % store the lower/upper limits induced by the
            % vertical/horizontal edges
            obj.Storage.GlobalMin.Grid.VertLim = grid_v_lim;
            obj.Storage.GlobalMin.Grid.HorzLim = grid_h_lim;

            % prepare some indices to speed up the branch-and-bound
            % algorithm for the global minimization problem

            % branching factor in the horizontal/vertical direction; each
            % rectangle will be divided into (n-1) * (m-1) sub-rectangles
            % with n * m resulting grid points
            n = obj.GlobalOptions.branching(1) + 1;
            m = obj.GlobalOptions.branching(2) + 1;

            grid_ind = reshape(1:(n * m), n, m);
            grid_bl = grid_ind(1:(n-1), 1:(m-1));
            grid_br = grid_ind(2:n, 1:(m-1));
            grid_tl = grid_ind(1:(n-1), 2:m);
            grid_tr = grid_ind(2:n, 2:m);

            % these indices will select points that are
            % bottom-left/bottom-right/top-left/top-right corners of the
            % rectangles in the grid
            obj.Storage.GlobalMin.Grid.Indices = struct;
            obj.Storage.GlobalMin.Grid.Indices.BottomLeft = grid_bl(:);
            obj.Storage.GlobalMin.Grid.Indices.BottomRight = grid_br(:);
            obj.Storage.GlobalMin.Grid.Indices.TopLeft = grid_tl(:);
            obj.Storage.GlobalMin.Grid.Indices.TopRight = grid_tr(:);

            % pre-compute the scan lines to save time
            obj.Storage.GlobalMin.Grid.PointHorzScanLine ...
                = repmat((1:n)', 1, m * poly_num);
            obj.Storage.GlobalMin.Grid.RectVertScanLine ...
                = repmat((1:m - 1)', 1, n * poly_num)';
            obj.Storage.GlobalMin.Grid.RectHorzScanLine ...
                = repmat((1:n - 1)', 1, m * poly_num);

            % sparse matrices for logical operations where
            % vertical/horizontal scan lines are combined to select grid
            % points and rectangles
            obj.Storage.GlobalMin.Grid.PointHorzScanLineCombMat ...
                = repmat(speye(m), poly_num, 1);
            obj.Storage.GlobalMin.Grid.VertScanLineCombMat = (kron( ...
                speye(poly_num), sparse( ...
                [(1:n-1)'; (2:n)'], [(1:n-1)'; (1:n-1)'], 1, n, n - 1)) ...
                * repmat(speye(n - 1), poly_num, 1))' > 0;

            obj.Storage.GlobalMin.Grid.HorzScanLineCombMat = (kron( ...
                speye(poly_num), sparse( ...
                [(1:m-1)'; (2:m)'], [(1:m-1)'; (1:m-1)'], 1, m, m - 1)) ...
                * repmat(speye(m - 1), poly_num, 1)) > 0;

            % sparse matrices that are used for selecting the corners of
            % selected rectangles
            obj.Storage.GlobalMin.Grid.Rect2GridMats = struct;
            obj.Storage.GlobalMin.Grid.Rect2GridMats.BottomLeft ...
                = sparse(obj.Storage.GlobalMin.Grid.Indices.BottomLeft, ...
                (1:(n - 1) * (m - 1))', 1, ...
                n * m, (n - 1) * (m - 1));
            obj.Storage.GlobalMin.Grid.Rect2GridMats.BottomRight ...
                = sparse( ...
                obj.Storage.GlobalMin.Grid.Indices.BottomRight, ...
                (1:(n - 1) * (m - 1))', 1, ...
                n * m, (n - 1) * (m - 1));
            obj.Storage.GlobalMin.Grid.Rect2GridMats.TopLeft ...
                = sparse(obj.Storage.GlobalMin.Grid.Indices.TopLeft, ...
                (1:(n - 1) * (m - 1))', 1, ...
                n * m, (n - 1) * (m - 1));
            obj.Storage.GlobalMin.Grid.Rect2GridMats.TopRight ...
                = sparse(obj.Storage.GlobalMin.Grid.Indices.TopRight, ...
                (1:(n - 1) * (m - 1))', 1, ...
                n * m, (n - 1) * (m - 1));

            % this flag is used to track if the function
            % obj.initializeSimplicialTestFuncs has been called
            obj.Storage.SimplicialTestFuncsInitialized = false;

            marg_num = length(obj.MarginalWeights);
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

        function [coup_indices, coup_probs, opt_qualities] ...
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
            %   coup_indices: the indices of atoms in the coupled discrete
            %   measure
            %   coup_probs: the probabilities of atoms in the coupled
            %   discrete measure
            %   opt_qualities: the optimal point in the quality space
            %   corresponding to each atom in the coupling

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            marg_num = length(obj.MarginalWeights);

            % retrieve the dual solution resulted from the cutting-plane
            % algorithm
            old_coup_indices = obj.Runtime.DualSolution.VertexIndices;
            old_coup_probs = obj.Runtime.DualSolution.Probabilities;

            % retrieve the old discretized marginals
            old_atoms_cell = cell(marg_num, 1);
            old_probs_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                old_atoms_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices;
                old_probs_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;
            end

            % set the new simplicial test functions
            obj.setSimplicialTestFuncs(args_cell);

            % the optimal couplings between the original discrete marginals
            % and the new discrete marginals where the cost function is
            % given by the Euclidean distance
            marg_coup_atoms_cell = cell(marg_num, 1);
            marg_coup_probs_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                new_atoms = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices;
                new_probs = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;

                % the cost function is the Euclidean distance
                dist_mat = pdist2(old_atoms_cell{marg_id}, ...
                    new_atoms, 'squaredeuclidean');

                % compute an optimal coupling between the original discrete
                % marginal and the new discrete marginal discrete optimal 
                % transport 
                [marg_coup_atoms_cell{marg_id}, ...
                    marg_coup_probs_cell{marg_id}] ...
                    = discrete_OT(old_probs_cell{marg_id}, ...
                    new_probs, dist_mat);
            end

            % perform discrete reassembly to get the new coupling
            [coup_indices, coup_probs] ...
                = discrete_reassembly(old_coup_indices, old_coup_probs, ...
                marg_coup_atoms_cell, marg_coup_probs_cell);

            % solve for the optimizers in the quality space
            coup_weights = obj.computeWeightedSumOfVertices(coup_indices);
            opt_qualities = obj.computeQuadMin(coup_weights);
        end

        function initializeSimplicialTestFuncs(obj)
            % Initialize some quantities related to the simplicial test 
            % functions of the marginals

            marg_num = length(obj.MarginalWeights);

            obj.Storage.MargVertNumList = zeros(marg_num, 1);

            weighted_vertices_cell = cell(marg_num, 1);
            deci_logical_cell = cell(marg_num, 1);
            
            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                obj.Storage.MargVertNumList(marg_id) = size( ...
                    marg.SimplicialTestFuncs.Vertices, 1);

                weighted_vertices_cell{marg_id} ...
                    = obj.MarginalWeights(marg_id) ...
                    * marg.SimplicialTestFuncs.Vertices;

                % the coefficient corresponding to the first test function
                % will not be included in the decision variable for
                % identification purposes
                deci_logical_cell{marg_id} = [0; ...
                    ones(obj.Storage.MargVertNumList(marg_id) - 1, 1)];
            end

            obj.Storage.TotalVertNum = sum(obj.Storage.MargVertNumList);
            obj.Storage.WeightedMargVertices ...
                = vertcat(weighted_vertices_cell{:});

            % compute the offset of the vertices in the marginals in the
            % matrix containing all vertices
            vert_num_cumsum = cumsum(obj.Storage.MargVertNumList);
            obj.Storage.MargVertNumOffsets = [0; ...
                vert_num_cumsum(1:end - 1)];

            % store the indices to place the decision variables in a vector
            % containing the coefficients of the test functions
            obj.Storage.DeciVarIndicesInTestFuncs = find( ...
                vertcat(deci_logical_cell{:}));

            % prepare information for facilitating the computation of the
            % concave part of the objective function in the global
            % minimization problem
            gm_vert_cell = cell(marg_num, 1);
            gm_ind_cell = cell(marg_num, 1);

            % compute the maximum number of vertices in the triangulations
            % of the marginals 
            gm_vert_num_max = max(obj.Storage.MargVertNumList);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                vertices = marg.SimplicialTestFuncs.Vertices;
                vert_num = size(vertices, 1);

                % add redundant vertices to the end for ease of computation
                gm_vert_cell{marg_id} = [(2 ...
                    * obj.MarginalWeights(marg_id)) ...
                    * vertices; zeros(gm_vert_num_max - vert_num, 2)];

                gm_ind_cell{marg_id} = [true(vert_num, 1); 
                    false(gm_vert_num_max - vert_num, 1)];
            end

            obj.Storage.GlobalMin.Concave = struct;
            obj.Storage.GlobalMin.Concave.lin = vertcat(gm_vert_cell{:});

            % indices for the non-redundant entries in the intercept
            % vector; this allows the values to be quickly filled into the
            % intercept vector
            obj.Storage.GlobalMin.Concave.intercept_indices ...
                = find(vertcat(gm_ind_cell{:}));

            obj.Storage.GlobalMin.Concave.vertex_num = gm_vert_num_max;

            obj.Storage.SimplicialTestFuncsInitialized = true;

            % updating the simplicial test functions will invalidate all
            % quantities in the runtime environment, thus all variables in
            % the runtime environment need to be flushed
            obj.Runtime = [];
        end

        function [coup_indices, coup_probs, opt_qualities] ...
                = generateHeuristicCoupling(obj, projection_dir)
            % Heuristically couple the marginals by first projecting all
            % two-dimensional vertices onto a fixed direction and then
            % apply comonotone coupling for the resulting one-dimensional
            % measures
            % Input:
            %   proj: vector representing the projection direction (default
            %   is [-1; 1], which means projecting onto the diagonal line)
            % Outputs:
            %   coup_indices: the indices of atoms in the coupled discrete
            %   measure
            %   coup_probs: the probabilities of atoms in the coupled
            %   discrete measure
            %   opt_qualities: the optimal point in the quality space
            %   corresponding to each atom in the coupling

            if ~exist('projection_dir', 'var') || isempty(projection_dir)
                projection_dir = [-1; 1];
            end

            marg_num = length(obj.MarginalWeights);

            % retrieve some information from the marginals
            atoms_cell = cell(marg_num, 1);
            probs_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                atoms_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices ...
                    * projection_dir;
                probs_cell{marg_id} = ...
                    obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;
            end

            [coup_indices, coup_probs] = comonotone_coupling( ...
                atoms_cell, probs_cell);

            coup_weights = obj.computeWeightedSumOfVertices(coup_indices);
            opt_qualities = obj.computeQuadMin(coup_weights);
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

            marg_num = length(obj.MarginalWeights);
            marg = obj.Marginals{marg_id};

            % divide the constant term by marg_num and add that to each of
            % the parametric functions
            vals = marg.evaluateWeightedSumOfSimplicialTestFuncs(pts, ...
                obj.Runtime.PrimalSolution.Coefficients{marg_id}, ...
                batch_size) + obj.Runtime.PrimalSolution.Constant ...
                / marg_num;
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

            total_vert_num = obj.Storage.GlobalMin.Concave.vertex_num ...
                * (marg_num - 1);
            intercept_indices = ...
                obj.Storage.GlobalMin.Concave.intercept_indices;
            intercept_indices = intercept_indices(intercept_indices ...
                <= total_vert_num);
            concave_lin = obj.Storage.GlobalMin.Concave.lin( ...
                1:total_vert_num, :);
            concave_intercept = -inf(total_vert_num, 1);
            concave_intercept(intercept_indices) = vertcat( ...
                obj.Runtime.PrimalSolution.Coefficients{1:end - 1});

            batch_num = ceil(pt_num / batch_size);
            vals_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                cur_pts = pts(((batch_id - 1) * batch_size ...
                    + 1):min(batch_id * batch_size, pt_num), :);
                cur_batch_size = size(cur_pts, 1);
                
                % only evaluate the first (marg_num - 1) transfer
                % functions, the last one is determined by the balance
                % condition
                vals_cell{batch_id} = sum(cur_pts .^ 2, 2) ...
                    .* obj.MarginalWeights(1:end - 1)';

                % matrix where each column represents one constituent
                % continuous piece-wise affine function corresponding to a
                % marginal evaluated with respect to each vertex in that
                % marginal (before the maximum is taken)
                concave_vals = reshape(concave_lin * cur_pts' ...
                    + concave_intercept, ...
                    obj.Storage.GlobalMin.Concave.vertex_num, ...
                    (marg_num - 1) * cur_batch_size);

                vals_cell{batch_id} = vals_cell{batch_id} ...
                    + reshape(-max(concave_vals, [], 1), marg_num - 1, ...
                    cur_batch_size)';

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

        function [coup_indices, coup_probs, coup_qualities, coup_costs] ...
                = generateDiscreteCoupling(obj, quality_vertices)
            % Generate a feasible dual solution by solving the discretized
            % version of the problem with fixed atoms in the quality space
            % Input:
            %   quality_vertices: fixed atoms in the quality space used for
            %   formulating a discrete version of the problem
            % Outputs:
            %   coup_indices: the atoms in the input space in the coupled 
            %   discrete measure represented as indices of vertices in the
            %   simplicial test functions
            %   coup_qualities: the atoms in the quality space in the 
            %   coupled discrete measure
            %   coup_costs: the corresponding values of the cost functions
            %   evaluated at the combinations of inputs and qualities
            %   coup_probs: the probabilities of atoms in the coupled
            %   discrete measure

            marg_num = length(obj.MarginalWeights);
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
                rhs_ineq_cell{marg_id} = obj.MarginalWeights(marg_id) ...
                    * sum(marg_grid_pts .* (marg_grid_pts ...
                    - 2 * quality_grid_pts), 2);

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

            % compute the minimization over the quality space
            weighted_sums = obj.computeWeightedSumOfVertices(coup_indices);
            [coup_qualities, coup_costs] = ...
                obj.computeQuadMin(weighted_sums);
        end

        function samps = randSampleFromOptJointDistr(obj, ...
                samp_num, rand_stream, batch_size)
            % Generate independent random sample from the joint
            % distribution consisting of the coupling of the discretized
            % marginals, the corresponding minimizers in the quality space,
            % the coupled continuous marginals, and the corresponding
            % continuous measure in the quality space
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
            [samps, atom_indices] ...
                = obj.doRandSampleFromPartialReassembly( ...
                (1:marg_num)', samp_num, rand_stream);
            
            % get the corresponding samples from the discrete measure on
            % the quality space by the atom indices
            samps.DiscreteQualities = ...
                obj.Runtime.DualSolution.QualityAtoms(atom_indices, :);

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

            [samps, atom_indices] ...
                = obj.doRandSampleFromPartialReassembly( ...
                marg_id, samp_num, rand_stream);

            coup = struct;
            coup.AgentTypes = samps.ContinuousPoints{marg_id};

            % get the corresponding samples from the discrete measure on
            % the quality space by the atom indices
            coup.Qualities = ...
                obj.Runtime.DualSolution.QualityAtoms(atom_indices, :);
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

            UB_disc_list = UB_disc_mat * obj.MarginalWeights ...
                + obj.QuadraticConstant;
            UB_cont_list = UB_cont_mat * obj.MarginalWeights ...
                + obj.QuadraticConstant;

            UB_disc = mean(UB_disc_list);
            UB_cont = mean(UB_cont_list);

            samps.UBDiscrete = UB_disc_list;
            samps.UBContinuous = UB_cont_list;
        end

        function EB = getMTErrorBoundBasedOnOT(obj, LSIP_tolerance)
            % Compute the error bound for the objective value of the 
            % matching for teams problem based on the Lipschitz constant of 
            % the cost functions and the optimal transport distances 
            % between the marginals and their discretizations
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
                EB = EB + 2 * obj.MarginalWeights(marg_id) ...
                    * obj.Quality.MaxNorm * marg.OT.Cost;
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
            end
        end
    end

    methods(Access = protected)

        function prepareRuntime(obj)
            % Prepare the runtime environment by initializing some
            % variables
            prepareRuntime@LSIPMinCuttingPlaneAlgo(obj);

            % initialize the cuts to be empty
            obj.Runtime.CutIndices = zeros(0, length(obj.MarginalWeights));
        end

        function [rect_lb, rect_ub, grid_pts, inner_grid_indices, ...
                outer_grid_indices, outer_rect_indices, ...
                outer_rect_corners] ...
                = computeGlobalMinGrids(obj, bound_rect_lb, bound_rect_ub)
            % Compute the relevant grids for the global minimization
            % problem. Given the bottom-left and top-right coordinates of a
            % rectangle, the rectangle is divided according to the
            % branching factor in obj.GlobalOptions. Subsequently, two
            % subsets of grid points are computed. The first subset of grid
            % points (the inner grid) are those that are contained in the
            % quality space. The second subset of grid points (the outer
            % grid) are those who are vertices of sub-rectangles that
            % overlap with the quality space. 
            % Inputs: 
            %   bound_rect_lb: the coordinate of the bottom-left corner of
            %   the bounding rectangle that will be subdivided
            %   bount_rect_ub: the coordinate of the top-right corner of
            %   the bounding rectangle that will be subdivided
            % Outputs: 
            %   rect_lb: two-column matrix containing the coordinates of
            %   the bottom-left corners of the sub-rectangles
            %   rect_ub: two-column matrix containing the coordinates of
            %   the top-right corners of the sub-rectangles
            %   grid_pts: two-column matrix containing the coordinates of
            %   all grid points
            %   inner_grid_indices: vector containing indices of the inner
            %   grid in output grid_pts
            %   outer_grid_indices: vector containing indices of the outer
            %   grid in output grid_pts
            %   outer_rect_indices: vector containing indices of the
            %   rectangles enclosed by the outer grid
            %   outer_rect_corners: four-column logical matrix for
            %   selecting the values on the corners of the selected
            %   rectangles from the outer grid

            n = obj.GlobalOptions.branching(1) + 1;
            m = obj.GlobalOptions.branching(2) + 1;
            poly_num = length(obj.Quality.Polytopes);

            % subdivide the bounding rectangle in the horizontal and
            % vertical directions; after subdivision the grid is n by m and
            % the resulting sub-rectangles are (n-1) by (m-1)
            grid_knot1 = linspace(bound_rect_lb(1), bound_rect_ub(1), n)';
            grid_knot2 = linspace(bound_rect_lb(2), bound_rect_ub(2), m)';
            [grid2, grid1] = meshgrid(grid_knot2, grid_knot1);
            grid_pts = [grid1(:), grid2(:)];

            % compute the coordinates of the bottom-left and top-right
            % corners of the sub-rectangles
            grid1_bl = grid1(1:n-1, 1:m-1);
            grid2_bl = grid2(1:n-1, 1:m-1);
            grid1_tr = grid1(2:n, 2:m);
            grid2_tr = grid2(2:n, 2:m);
            rect_lb = [grid1_bl(:), grid2_bl(:)];
            rect_ub = [grid1_tr(:), grid2_tr(:)];
            
            Grid = obj.Storage.GlobalMin.Grid;

            % for each vertical line, compute the upper end point that is
            % in each polytope; this takes into account the upper limits
            % induced by horizontal edges as well
            vert_ub = min(reshape(min(reshape(grid_knot1' ...
                .* Grid.VertUB.coef + Grid.VertUB.intercept, ...
                Grid.VertUB.num, poly_num * n), [], 1), poly_num, n)', ...
                Grid.VertLim(:, 2)');

            % for each vertical line, compute the lower end point that is
            % in each polytope; this takes into account the lower limits
            % induced by horizontal edges as well
            vert_lb = max(reshape(max(reshape(grid_knot1' ...
                .* Grid.VertLB.coef + Grid.VertLB.intercept, ...
                Grid.VertLB.num, poly_num * n), [], 1), poly_num, n)', ...
                Grid.VertLim(:, 1)');

            % exclude vertical lines lie outside the horizontal limits
            % induced by vertical edges of each polytope
            vert_line_oob = grid_knot1 < Grid.HorzLim(:, 1)' ...
                | grid_knot1 > Grid.HorzLim(:, 2)';
            vert_lb(vert_line_oob) = inf;

            % for each horizontal line, compute the upper end point that is
            % in each polytope; this takes into account the upper limits
            % induced by vertical edges as well
            horz_ub = min(reshape(min(reshape(grid_knot2' ...
                .* Grid.HorzUB.coef + Grid.HorzUB.intercept, ...
                Grid.HorzUB.num, poly_num * m), [], 1), poly_num, m)', ...
                Grid.HorzLim(:, 2)');

            % for each horizontal line, compute the lower end point that is
            % in each polytope; this takes into account the lower limits
            % induced by vertical edges as well
            horz_lb = max(reshape(max(reshape(grid_knot2' ...
                .* Grid.HorzLB.coef + Grid.HorzLB.intercept, ...
                Grid.HorzLB.num, poly_num * m), [], 1), poly_num, m)', ...
                Grid.HorzLim(:, 1)');

            % exclude horizontal lines lie outside the vertical limits
            % induced by horizontal edges of each polytope
            horz_line_oob = grid_knot2 < Grid.VertLim(:, 1)' ...
                | grid_knot2 > Grid.VertLim(:, 2)';
            horz_lb(horz_line_oob) = inf;

            % normalize each interval according to the lower end point of
            % the bounding rectangle and the mesh size; the resulting
            % intervals are represented by (possibly fractional) indices
            % with respect to the grid
            vert_min = bound_rect_lb(2);
            vert_mesh = grid_knot2(2) - grid_knot2(1);
            horz_min = bound_rect_lb(1);
            horz_mesh = grid_knot1(2) - grid_knot1(1);
            vert_ub_norm = (vert_ub(:) - vert_min) / vert_mesh + 1;
            vert_lb_norm = (vert_lb(:) - vert_min) / vert_mesh + 1;
            horz_ub_norm = (horz_ub(:) - horz_min) / horz_mesh + 1;
            horz_lb_norm = (horz_lb(:) - horz_min) / horz_mesh + 1;

            % taking the floor with result in indices for the outer
            % sub-rectangles
            vert_lb_outer = floor(vert_lb_norm);
            horz_lb_outer = floor(horz_lb_norm);

            % lower limits of empty intervals are set to inf to completely
            % remove them; else something like [1.5, 1.4] will create a
            % non-empty line segment
            vert_lb_outer(vert_lb_norm > vert_ub_norm) = inf;
            horz_lb_outer(horz_lb_norm > horz_ub_norm) = inf;

            % scan lines for computing the inner grid; for the inner grid,
            % a horizontal scan is sufficient
            grid_scanline = Grid.PointHorzScanLine;
            inner_grid_indices = (((grid_scanline <= horz_ub_norm') ...
                & (grid_scanline >= horz_lb_norm')) ...
                * Grid.PointHorzScanLineCombMat) > 0;
            inner_grid_indices = inner_grid_indices(:);

            % vertical and horizontal scan lines for computing the outer
            % grid (these will check whether each edge of the sub-rectangle
            % intersects with each polytope)
            rect_vert_scanline = Grid.RectVertScanLine;
            rect_horz_scanline = Grid.RectHorzScanLine;

            % vertical edges of the sub-rectangles that intersect with
            % each polytope
            rect_vert_overlap = (rect_vert_scanline <= vert_ub_norm) ...
                & (rect_vert_scanline >= vert_lb_outer);

            % horizontal edges of the sub-rectangles that intersect with
            % each polytope
            rect_horz_overlap = (rect_horz_scanline <= horz_ub_norm') ...
                & (rect_horz_scanline >= horz_lb_outer');

            outer_rect_indices = ((Grid.VertScanLineCombMat ...
                * rect_vert_overlap) ...
                + (rect_horz_overlap * Grid.HorzScanLineCombMat)) > 0;
            outer_rect_indices = outer_rect_indices(:);

            % final step: check the chosen vertices from the polytopes for
            % a special case where a polytope is entirely contained in a
            % sub-rectangle

            % compute the row and column indices of these vertices
            poly_vert_rect = floor((Grid.AdditionalPoints ...
                - [horz_min, vert_min]) ./ [horz_mesh, vert_mesh] + 1);
            poly_vert_rect = poly_vert_rect( ...
                poly_vert_rect(:, 1) >= 1 ...
                & poly_vert_rect(:, 1) <= n - 1 ...
                & poly_vert_rect(:, 2) >= 1 ...
                & poly_vert_rect(:, 2) <= m - 1, :);

            % make sure that these rectangles are chosen in the outer grid
            outer_rect_indices((poly_vert_rect(:, 2) - 1) * (n - 1) ...
                + poly_vert_rect(:, 1)) = 1;

            % select grid points that are the corners of the selected
            % rectangles
            outer_rect_corners = ...
                [Grid.Rect2GridMats.BottomLeft * outer_rect_indices, ...
                Grid.Rect2GridMats.BottomRight * outer_rect_indices, ...
                Grid.Rect2GridMats.TopLeft * outer_rect_indices, ...
                Grid.Rect2GridMats.TopRight * outer_rect_indices];

            % the outer grid is the union of those corners
            outer_grid_indices = sum(outer_rect_corners, 2) > 0;

            % compute the logical indices for selecting the corners of the
            % selected rectangles from the outer grid
            outer_rect_corners = logical( ...
                outer_rect_corners(outer_grid_indices, :));
        end

        function model = generateInitialMinModel(obj)
            % Generate the initial linear programming model for gurobi
            % Output:
            %   model: struct containing the linear programming model in
            %   gurobi

            model = struct;
            marg_num = length(obj.MarginalWeights);
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

            % the constant in the objective is the sum of the quadratic
            % constants that do not affect the matching
            model.objcon = -obj.QuadraticConstant;

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
            %   global minimization problem in the form of vertex indices
            %   in the triangulations for each of the input measures as
            %   well as the the corresponding points in the two-dimensional
            %   space

            % disassemble the vector resulted from solving the LP problem

            % the first component corresponds to the constant term; the
            % second component onwards will be filled into the
            % corresponding entries of the vector storing the coefficients
            % corresponding to the test functions
            objective_const = vec(1);
            vert_vals = zeros(obj.Storage.TotalVertNum, 1);
            vert_vals(obj.Storage.DeciVarIndicesInTestFuncs) = vec(2:end);

            % call the branch-and-bound algorithm
            [pool_inds, pool_pts, pool_vals, LB] ...
                = obj.computeGlobalMinBnB(vert_vals, objective_const);

            min_lb = LB;

            % ignore those approximate minimizers whose objective values
            % are non-negative since they do not generate cuts
            pool_neg_list = pool_vals < 0;

            optimizers = struct;
            optimizers.vertex_indices = pool_inds(pool_neg_list, :);
            optimizers.points = pool_pts(pool_neg_list, :);
        end

        function [pool_inds, pool_pts, pool_vals, LB] ...
                = computeGlobalMinBnB(obj, vert_vals, objective_const)
            % Approximately solve the global minimization problem via the
            % branch-and-bound algorithm
            % Input:
            %   vert_vals: values of the simplicial test functions at the
            %   vertices concatenated into a vector
            %   objective_const: constant part of the objective function
            % Outputs:
            %   pool_inds: matrix containing approximate optimizers of the
            %   global minimization problem in the form of vertex indices
            %   in the triangulations for each of the input measures; the
            %   algorithm will retain up to obj.GlobalOptions.pool_size
            %   approximate optimizers
            %   pool_pts: two-column matrix containing the approximate
            %   optimizers in the two-dimensional space; this will be used
            %   to generate linear inequality constraints
            %   pool_vals: the corresponding objective values of the
            %   approximate optimizers
            %   LB: the best lower bound found throughout the algorithm

            % open the log file
            if ~isempty(obj.GlobalOptions.log_file)
                log_file = fopen(obj.GlobalOptions.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- Function call starts ---\n');
            end

            % set up variables in the runtime environment
            obj.Runtime.GlobalMin = struct;
            obj.Runtime.GlobalMin.Concave = struct;

            % fill in the intercept values for the concave part of the
            % objective
            total_vert_num = size(obj.Storage.GlobalMin.Concave.lin, 1);
            concave_intercept = -inf(total_vert_num, 1);
            concave_intercept( ...
                obj.Storage.GlobalMin.Concave.intercept_indices) ...
                = vert_vals;
            obj.Runtime.GlobalMin.Concave.intercept = concave_intercept;

            % retrieve options for the global optimization problem
            tolerance = obj.GlobalOptions.tolerance;
            batch_size = obj.GlobalOptions.batch_size;
            pool_size = obj.GlobalOptions.pool_size;
            termination_thres = obj.GlobalOptions.termination_thres;

            % lower/upper bounds of the rectangles
            rect_lb = obj.Storage.GlobalMin.Grid.InitRectangle.lb';
            rect_ub = obj.Storage.GlobalMin.Grid.InitRectangle.ub';

            % lower/upper bounds of the current rectangle to branch
            cur_rect_lb = rect_lb';
            cur_rect_ub = rect_ub';

            % points where the objective is evaluated in the initial step
            pts = [cur_rect_lb'; [cur_rect_lb(1), cur_rect_ub(2)]; ...
                [cur_rect_ub(1), cur_rect_lb(2)]; cur_rect_ub'; ...
                obj.Quality.Vertices];

            [vals, min_inds] = obj.evaluateGlobalMinConcave(pts, ...
                batch_size);
            vals = vals + sum(pts .^ 2, 2) - objective_const;

            % vertices of the polytopes are the initial candidates
            pool_vals = vals(5:end);
            pool_inds = min_inds(5:end, :);
            pool_pts = pts(5:end, :);

            % sort the candidates and make sure that the pool size does not
            % exceed the desired pool size
            [pool_vals, pool_vals_order] = sort(pool_vals, 'ascend');
            pool_inds = pool_inds(pool_vals_order, :);
            pool_pts = pool_pts(pool_vals_order, :);

            % remove duplicate candidates that have the same vertex indices
            [pool_inds, pool_unique_indices] = unique(pool_inds, ...
                    'rows', 'stable');
            pool_vals = pool_vals(pool_unique_indices);

            pool_vals_retain_num = min(length(pool_vals), pool_size);

            % pool containing the objective values of the best candidates
            % found so far (sorted into ascending order)
            pool_vals = pool_vals(1:pool_vals_retain_num);

            % pool containing the best candidates found so far (sorted into
            % ascending order in the corresponding objective values)
            pool_inds = pool_inds(1:pool_vals_retain_num, :);
            pool_pts = pool_pts(1:pool_vals_retain_num, :);

            % lower bounds on the objective value in each outer rectangle
            rect_val_lb = min(vals(1:4)) ...
                - sum((cur_rect_ub - cur_rect_lb) .^ 2) / 4;

            % the global upper bound which is given by the objective value
            % of the best candidate found so far
            UB = min(pool_vals);

            % the global lower bound is given by the minimum of lower
            % bounds within each outer rectangle
            LB = rect_val_lb;

            % iteration counter
            iter = 1;

            % loop until either the difference between the upper and lower
            % bounds is below the desired tolerance or a candidate with
            % objective value less than the termination threshold has been
            % found
            while true
                if obj.GlobalOptions.display
                    % display the current upper and lower bounds
                    fprintf(['%s: ' ...
                        'iteration %4d: LB = %10.4f, UB = %10.4f, ' ...
                        'UB - LB = %10.6f\n'], ...
                        class(obj), iter, LB, UB, UB - LB);
                end

                % logging
                if ~isempty(obj.GlobalOptions.log_file)
                    fprintf(log_file, ['%s: iteration %4d: ' ...
                        'LB = %10.4f, UB = %10.4f, ' ...
                        'UB - LB = %10.6f\n'], ...
                        class(obj), iter, LB, UB, UB - LB);
                end

                if (UB - LB <= tolerance) || (UB < termination_thres)
                    break;
                end
                
                % compute the grid and the associated quantities
                [new_rect_lb, new_rect_ub, new_pts, inner_grid_indices, ...
                    outer_grid_indices, outer_square_indices, ...
                    outer_rect_corners] ...
                    = obj.computeGlobalMinGrids(cur_rect_lb, cur_rect_ub);

                % compute the objective function only at the outer grid
                % points
                pts_to_evaluate = new_pts(outer_grid_indices, :);
                [vals, min_inds] = obj.evaluateGlobalMinConcave( ...
                    pts_to_evaluate, batch_size);

                % add in the quadratic part of the objective function and
                % the constant part
                vals = vals + sum(pts_to_evaluate .^ 2, 2) ...
                    - objective_const;

                % add the new candidates into the pool
                new_cand_indices = inner_grid_indices(outer_grid_indices);
                pool_inds = [pool_inds; ...
                    min_inds(new_cand_indices, :)]; %#ok<AGROW>
                pool_pts = [pool_pts; ...
                    pts_to_evaluate(new_cand_indices, :)]; %#ok<AGROW> 
                pool_vals = [pool_vals; ...
                    vals(new_cand_indices)]; %#ok<AGROW> 

                % sort the candidates and make sure that the pool size does
                % not exceed the desired pool size
                [pool_vals, pool_vals_order] = sort(pool_vals, 'ascend');
                pool_inds = pool_inds(pool_vals_order, :);
                pool_pts = pool_pts(pool_vals_order, :);

                % remove duplicate candidates that have the same vertex
                % indices
                [pool_inds, pool_unique_indices] = unique(pool_inds, ...
                    'rows', 'stable');
                pool_vals = pool_vals(pool_unique_indices);

                pool_vals_retain_num = min(length(pool_vals), pool_size);
                pool_vals = pool_vals(1:pool_vals_retain_num);
                pool_inds = pool_inds(1:pool_vals_retain_num, :);
                pool_pts = pool_pts(1:pool_vals_retain_num, :);

                % update the global upper bound
                UB = pool_vals(1);

                % compute the lower bound within each outer rectangle by
                % computing the minimum objective value on the four corners
                % minus the relaxation error
                new_rect_val_lb = min([vals(outer_rect_corners(:, 1)), ...
                    vals(outer_rect_corners(:, 2)), ...
                    vals(outer_rect_corners(:, 3)), ...
                    vals(outer_rect_corners(:, 4))], [], 2) ...
                    - sum((new_rect_ub(1, :) - new_rect_lb(1, :)) .^ 2) ...
                    / 4;

                % add the new rectangles into the list
                rect_lb = [rect_lb(2:end, :); ...
                    new_rect_lb(outer_square_indices, :)];
                rect_ub = [rect_ub(2:end, :); ...
                    new_rect_ub(outer_square_indices, :)];

                % sort the lower bounds within outer rectangle into
                % ascending order
                [rect_val_lb, rect_val_lb_order] ...
                    = sort([rect_val_lb(2:end); new_rect_val_lb], ...
                    'ascend');

                % the number of rectangles to retain; the rest have lower
                % bounds no less than the current upper bound and thus can
                % be discarded
                rect_val_lb_retain_num = find(rect_val_lb < UB, 1, 'last');
                rect_val_lb_retain_list = rect_val_lb_order( ...
                    1:rect_val_lb_retain_num);

                % update the outer rectangles
                rect_lb = rect_lb(rect_val_lb_retain_list, :);
                rect_ub = rect_ub(rect_val_lb_retain_list, :);
                rect_val_lb = rect_val_lb(1:rect_val_lb_retain_num);

                % choose the outer rectangle with the smallest lower bound
                % to branch in the next iteration
                cur_rect_lb = rect_lb(1, :)';
                cur_rect_ub = rect_ub(1, :)';
                LB = rect_val_lb(1);

                % update the iteration counter
                iter = iter + 1;
            end

            % close the log file
            if ~isempty(obj.GlobalOptions.log_file)
                fprintf(log_file, '--- Function call ends ---\n\n');
                fclose(log_file);
            end
        end

        function [vals, min_inds] ...
                = evaluateGlobalMinConcave(obj, pts, batch_size)
            % Evaluate the concave part of the objective function in the
            % global minimization problem. This function is only invoked
            % after setting obj.Runtime.GlobalMin.Concave.intercept.
            % Evaluation is done in batches if necessary to avoid excessive
            % memory usage. 
            % Inputs:
            %   pts: two-column matrix containing the input points
            %   batch_size: maximum number of inputs to be evaluated
            %   together via a vectorized routine (default is 10300)
            % Output:
            %   vals: the computed function values
            %   min_inds: the selected minimum indices in each continuous
            %   piece-wise affine part of the function
            
            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 10300;
            end

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);
            min_inds_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                [vals_cell{batch_id}, min_inds_cell{batch_id}] ...
                    = obj.doEvaluateGlobalMinConcave( ...
                    pts(((batch_id - 1) * batch_size + 1):min(batch_id ...
                    * batch_size, input_num), :));
            end

            vals = vertcat(vals_cell{:});
            min_inds = vertcat(min_inds_cell{:});
        end

        function [vals, min_inds] = doEvaluateGlobalMinConcave(obj, pts)
            % Evaluate the concave part of the objective function in the
            % global minimization problem. This function is only invoked
            % after setting obj.Runtime.GlobalMin.Concave.intercept.
            % Input:
            %   pts: two-column matrix containing the input points
            % Output:
            %   vals: the computed function values
            %   min_inds: the selected minimum indices in each continuous
            %   piece-wise affine part of the function
            
            pt_num = size(pts, 1);
            marg_num = length(obj.MarginalWeights);

            % matrix where each column represents one constituent
            % continuous piece-wise affine function corresponding to a
            % marginal evaluated with respect to each vertex in that
            % marginal (before the maximum is taken)
            vals_per_index = reshape(obj.Storage.GlobalMin.Concave.lin ...
                * pts' + obj.Runtime.GlobalMin.Concave.intercept, ...
                obj.Storage.GlobalMin.Concave.vertex_num, ...
                marg_num * pt_num);

            [vals_min, min_inds_unwrapped] = max(vals_per_index, [], 1);

            vals = -sum(reshape(vals_min, marg_num, pt_num), 1)';
            min_inds = reshape(min_inds_unwrapped, marg_num, pt_num)';
        end

        function addConstraints(obj, optimizers)
            % Given a collection of approximate optimizers from the global
            % minimization oracle, generate and add the corresponding
            % linear constraints
            % Inputs:
            %   optimizers: output of the method callGlobalOracle

            constr_num = size(optimizers.vertex_indices, 1);
            marg_num = length(obj.MarginalWeights);
            
            col_indices = optimizers.vertex_indices' ...
                + obj.Storage.MargVertNumOffsets;

            % first generate a matrix containing all test functions (each
            % row corresponds to an approximate optimizer, each column
            % corresponds to a test function)
            A_full = sparse(repelem((1:constr_num)', marg_num, 1), ...
                col_indices(:), 1, constr_num, obj.Storage.TotalVertNum);
            
            % filter out those test functions whose coefficients are not
            % included in the decision variable, then prepend a column of
            % 1
            A_new = [sparse(ones(constr_num, 1)), ...
                A_full(:, obj.Storage.DeciVarIndicesInTestFuncs)];

            gm_concave_indices = optimizers.vertex_indices' ...
                + (0:marg_num - 1)' ...
                * obj.Storage.GlobalMin.Concave.vertex_num;
            sum_over_marg_mat = sparse(repelem((1:constr_num)', ...
                marg_num, 1), 1:constr_num * marg_num, 1);
            rhs_new = sum(optimizers.points .^ 2, 2) ...
                - sum(optimizers.points .* ...
                (sum_over_marg_mat * obj.Storage.GlobalMin.Concave.lin( ...
                gm_concave_indices, :)), 2);
            
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

            % add the added indices to the runtime environment
            obj.Runtime.CutIndices = [obj.Runtime.CutIndices; ...
                optimizers.vertex_indices];

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
                obj.Runtime.CutIndices ...
                    = obj.Runtime.CutIndices(keep_list, :);
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
            primal_sol.Constant = vec(1) - violation;

            vert_coefs = zeros(obj.Storage.TotalVertNum, 1);
            vert_coefs(obj.Storage.DeciVarIndicesInTestFuncs) = vec(2:end);
            marg_num = length(obj.MarginalWeights);
            primal_sol.Coefficients = cell(marg_num, 1);

            for marg_id = 1:marg_num
                primal_sol.Coefficients{marg_id} = vert_coefs( ...
                    obj.Storage.MargVertNumOffsets(marg_id) ...
                    + (1:obj.Storage.MargVertNumList(marg_id)));
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
            dual_sol.VertexIndices = obj.Runtime.CutIndices(pos_list, :);
            weighted_sum = obj.computeWeightedSumOfVertices( ...
                dual_sol.VertexIndices);
            dual_sol.QualityAtoms = obj.computeQuadMin(weighted_sum);
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
            %   disc_atom_index_samps: vector containing the indices of
            %   atoms in the discrete coupling

            % first make sure that the semi-discrete OT problems are solved
            if ~isfield(obj.Storage, 'OTComputed') ...
                    || ~obj.Storage.OTComputed
                obj.performReassembly();
            end

            marg_num = length(obj.MarginalWeights);
            n = length(obj.Runtime.DualSolution.Probabilities);
            
            % generate random indices of the atoms according to the
            % probabilities
            disc_atom_index_samps = randsample(rand_stream, n, ...
                samp_num, true, obj.Runtime.DualSolution.Probabilities);

            samps = struct;
            samps.DiscreteIndices ...
                = obj.Runtime.DualSolution.VertexIndices( ...
                disc_atom_index_samps, :);
            samps.DiscretePoints = cell(marg_num, 1);
            samps.ContinuousPoints = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_atoms = marg.SimplicialTestFuncs.Vertices;

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

