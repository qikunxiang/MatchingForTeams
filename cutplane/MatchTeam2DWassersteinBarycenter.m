classdef MatchTeam2DWassersteinBarycenter < LSIPMinCuttingPlaneAlgo
    % Class for matching for teams problems with two-dimensional marginals, two-dimensional quality space, and quadratic cost
    % functions. The problem is solved via parametrizing the transfer functions on the quality space.

    properties(Constant)
        % numerical tolerance for deciding whether a point is inside a polytope
        INSIDE_TOLERANCE = 1e-12;
    end

    properties(GetAccess = public, SetAccess = protected)
        % cell array containing the marginals
        Marginals;

        % vector containing the weights
        MarginalWeights;

        % struct containing information about the quality space
        Quality;

        % constant part of the quadratic cost functions that does not affect the matching
        QuadraticConstant = 0;
    end
    
    methods(Access = public)
        function obj = MatchTeam2DWassersteinBarycenter( ...
                marginals, ...
                weights, ...
                quality_cell, ...
                quality_testfuncs, ...
                varargin)
            % Constructor method
            % Inputs:
            %   marginals: cell array containing objects of ProbMeas2DConvexPolytope
            %   weights: vector containing the weights corresponding to the marginals, the weights will sum up to 1
            %   quality_cell: cell array where each cell contains a two-column matrix indicating the vertices of a polytope
            %   quality_testfuncs: cell array containing the inputs to the function obj.setQualitySimplicialTestFuncs

            obj@LSIPMinCuttingPlaneAlgo(varargin{:});

            % set the default options for semi-discrete optimal transport
            if ~isfield(obj.Options, 'OT') || isempty(obj.Options.OT)
                obj.Options.OT = struct;
            end

            if ~isfield(obj.Options.OT, 'angle_num') || isempty(obj.Options.OT.angle_num)
                obj.Options.OT.angle_num = repmat({[]}, length(weights), 1);
            elseif ~iscell(obj.Options.OT.angle_num)
                obj.Options.OT.angle_num = repmat({obj.Options.OT.angle_num}, length(weights), 1);
            end

            if ~isfield(obj.Options.OT, 'pp_angle_num') || isempty(obj.Options.OT.pp_angle_num)
                obj.Options.OT.pp_angle_num = repmat({[]}, length(weights), 1);
            elseif ~iscell(obj.Options.OT.pp_angle_num)
                obj.Options.OT.pp_angle_num = repmat({obj.Options.OT.pp_angle_num}, length(weights), 1);
            end

            if ~isfield(obj.Options.OT, 'normalize_gradient') || isempty(obj.Options.OT.normalize_gradient)
                obj.Options.OT.normalize_gradient = repmat({[]}, length(weights), 1);
            elseif ~iscell(obj.Options.OT.normalize_gradient)
                obj.Options.OT.normalize_gradient = repmat({obj.Options.OT.normalize_gradient}, length(weights), 1);
            end

            if ~isfield(obj.Options.OT, 'optimization_options')
                obj.Options.OT.optimization_options = cell(length(weights), 1);
            elseif ~iscell(obj.Options.OT.optimization_options)
                optim_options = obj.Options.OT.optimization_options;
                obj.Options.OT.optimization_options = cell(length(weights), 1);
                obj.Options.OT.optimization_options(:, :) = {optim_options};
            end

            if ~isfield(obj.Options.OT, 'prob_thres') || isempty(obj.Options.OT.prob_thres)
                % threshold to filter out atoms with probabilities close to 0 to prevent numerical issues in the computation of
                % semi-discrete OT
                obj.Options.OT.prob_thres = repmat({0}, length(weights), 1);
            elseif ~iscell(obj.Options.OT.prob_thres)
                obj.Options.OT.prob_thres = repmat({obj.Options.OT.prob_thres}, length(weights), 1);
            end

            % set the default options for semi-discrete Wasserstein-2 optimal transport
            if ~isfield(obj.Options, 'W2OT') || isempty(obj.Options.W2OT)
                obj.Options.W2OT = struct;
            end

            if ~isfield(obj.Options.W2OT, 'optimization_options')
                obj.Options.W2OT.optimization_options = cell(length(weights), 1);
            elseif ~iscell(obj.Options.W2OT.optimization_options)
                optim_options = obj.Options.W2OT.optimization_options;
                obj.Options.W2OT.optimization_options = cell(length(weights), 1);
                obj.Options.W2OT.optimization_options(:, :) = {optim_options};
            end

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

            % set the default options for reducing constraints; two thresholds are used for determining which constraints to remove 
            % based on their slackness values: whenever the (negative) slackness is above obj.Options.reduce.thres or its quantile
            % among all the non-tight constraints is below obj.Options.reduce.thres_quantile, the constraint is removed
            if ~isfield(obj.Options.reduce, 'thres_quantile') || isempty(obj.Options.reduce.thres_quantile)
                obj.Options.reduce.thres_quantile = 1;
            end

            % minimum value of slackness (in absolute value) for a constraint to be considered not tight where a constraint is only
            % flagged as removable only if |slack| > obj.Options.reduce.min_slack
            if ~isfield(obj.Options.reduce, 'min_slack') || isempty(obj.Options.reduce.min_slack)
                obj.Options.reduce.min_slack = 0;
            end

            % boolean indicating whether the initial constraints should be preserved throughout the cutting-plane algorithm
            if ~isfield(obj.Options.reduce, 'preserve_init_constr') || isempty(obj.Options.reduce.preserve_init_constr)
                obj.Options.reduce.preserve_init_constr = true;
            end

            % set the default option for the index of the discrete measure candidate on the quality space
            if ~isfield(obj.Options, 'discmeas_cand_index') || isempty(obj.Options.discmeas_cand_index)
                obj.Options.discmeas_cand_index = 1;
            end

            % set the default option for the sparse parametrization; if sparse_parametrization is true, rather than using a constant
            % intercept for the test functions on the support of a marginal, use the full set of the simplicial test functions; this
            % increases the sparsity of the constraint matrix in the LP and improves numerical stability
            if ~isfield(obj.Options, 'sparse_parametrization') || isempty(obj.Options.sparse_parametrization)
                obj.Options.sparse_parametrization = false;
            end

            % set the default options for the minimum probability value for an atom in the dual LSIP solution to be kept in the sense
            % that an atom is kept in the dual LSIP solution only if probability > obj.Options.dual_prob_thres
            if ~isfield(obj.Options, 'dual_prob_thres') || isempty(obj.Options.dual_prob_thres)
                obj.Options.dual_prob_thres = 0;
            end

            % set the default options for the global minimization oracle
            if ~isfield(obj.GlobalOptions, 'pool_size') || isempty(obj.GlobalOptions.pool_size)
                obj.GlobalOptions.pool_size = 100;
            end

            if ~isfield(obj.GlobalOptions, 'boundary_buffer') || isempty(obj.GlobalOptions.boundary_buffer)
                % this value controls how "close" approximate global minimizers can get; some approximate global minimizers are
                % represented as convex combinations of two or three vertices in the triangulation on the quality space, and when any
                % coefficient in a convex combination is below this threshold it is regarded as invalid; effectively, an approximate
                % global minimizer in the interior of a triangle that is too close to an edge will be discarded, and similarly, an
                % approximate global minimizer in the relative interior of an edge of a triangle that is too close to a vertex will be
                % discarded; this guarantees that there will not be approximate global minimizers that are almost identical to prevent
                % numerical issues
                obj.GlobalOptions.boundary_buffer = 1e-4;
            end

            if ~isfield(obj.GlobalOptions, 'objective_threshold') || isempty(obj.GlobalOptions.objective_threshold)
                % this value controls the threshold below which an approximate optimizer will be used to generate new cuts; using 0 as
                % the threshold might result in a large number of new cuts that are close to being non-violations
                obj.GlobalOptions.objective_threshold = -1e-6;
            end

            if ~isfield(obj.GlobalOptions, 'display') || isempty(obj.GlobalOptions.display)
                obj.GlobalOptions.display = true;
            end

            if ~isfield(obj.GlobalOptions, 'log_file') || isempty(obj.GlobalOptions.log_file)
                obj.GlobalOptions.log_file = '';
            end

            marg_num = length(weights);
            assert(length(marginals) == marg_num, 'input mis-specified');
            assert(abs(sum(weights) - 1) < 1e-12, 'weights do not sum up to 1');
            
            % normalize the weights to remove numerical inaccuracies
            weights = weights / sum(weights);

            obj.Marginals = marginals;
            obj.MarginalWeights = weights;

            % compute the constant terms in the cost functions that are related to the quadratic expectation with respect to the
            % marginals
            quad_consts = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                quad_consts(marg_id) = sum(diag(obj.Marginals{marg_id}.SecondMomentMat));
            end

            obj.QuadraticConstant = quad_consts' * obj.MarginalWeights;

            obj.Quality = struct;
            poly_num = length(quality_cell);
            obj.Quality.Polytopes = cell(poly_num, 1);

            % number of vertices in each polytope
            vert_num_list = zeros(poly_num, 1);

            % store all vertices of the polytopes which will be used in the first iteration of the global minimization algorithm
            vert_cell = cell(poly_num, 1);

            % store one vertex from each polytope in order to compute the rectangles that overlap with the polytopes; this is to
            % resolve an edge case where the entire polytope is contained in a rectangle
            first_vert = zeros(poly_num, 2);

            for poly_id = 1:poly_num
                vertices = quality_cell{poly_id};

                % compute the convex hull of the vertices
                ccorder = convhull(vertices(:, 1), vertices(:, 2), 'Simplify', true);

                vert_num_list(poly_id) = length(ccorder) - 1;

                % rearrange the vertices in counterclockwise order
                vertices_cc = vertices(ccorder(1:end - 1), :);

                obj.Quality.Polytopes{poly_id} = vertices_cc;
                vert_cell{poly_id} = vertices_cc;
                first_vert(poly_id, :) = vertices_cc(1, :);
            end

            % store all vertices of the polytopes
            obj.Quality.Vertices = unique(vertcat(vert_cell{:}), 'rows');

            % store the bounding box of the quality space
            obj.Quality.BoundingBox = [min(obj.Quality.Vertices, [], 1); max(obj.Quality.Vertices, [], 1)];

            % store the maximum norm of points in the quality space
            obj.Quality.MaxNorm = sqrt(max(sum(obj.Quality.Vertices .^ 2, 2)));

            % compute all hyperplanes characterizing the polytopes in the form of inequalities w' * z <= b
            hp_w_cell = cell(poly_num, 1);
            hp_b_cell = cell(poly_num, 1);

            % hyperplanes with redundant inequality as padding such that each polytope has the same number of inequalities
            hp_w_cell_rd = cell(poly_num, 1);
            hp_b_cell_rd = cell(poly_num, 1);

            % compute the maximum number of vertices among the polytopes; this will be used to add redundant inequalities such that all
            % polytopes are characterized by the same number of inequalities
            vert_num_max = max(vert_num_list);

            for poly_id = 1:poly_num
                % retrieve the vertices in the polytopes and their next counterclockwise neighbor
                vertices_cc1 = obj.Quality.Polytopes{poly_id};
                vertices_cc2 = vertices_cc1([2:end, 1], :);
                vert_num = vert_num_list(poly_id);

                % compute the weights and the intercepts of the hyperplanes characterizing the polytope
                hp_w_cell{poly_id} = [vertices_cc2(:, 2) - vertices_cc1(:, 2), vertices_cc1(:, 1) - vertices_cc2(:, 1)];
                hp_b_cell{poly_id} = vertices_cc2(:, 2) .* vertices_cc1(:, 1) - vertices_cc1(:, 2) .* vertices_cc2(:, 1);

                % correct numerical inaccuracies and make sure that the vertices themselves are contained in all of the half-spaces
                hp_b_cell{poly_id} = max(hp_b_cell{poly_id}, max(hp_w_cell{poly_id} * vertices_cc1', [], 2));

                % add redundant inequalities 0' * z <= 1 below such that the number of inequalities for each polytope is the same
                hp_w_cell_rd{poly_id} = [hp_w_cell{poly_id}; zeros(vert_num_max - vert_num, 2)];
                hp_b_cell_rd{poly_id} = [hp_b_cell{poly_id}; ones(vert_num_max - vert_num, 1)];
            end

            % store the information about the hyperplanes characterizing the polytopes
            obj.Quality.Hyperplane = struct;
            obj.Quality.Hyperplane.num = vert_num_max;
            obj.Quality.Hyperplane.w = vertcat(hp_w_cell_rd{:});
            obj.Quality.Hyperplane.b = vertcat(hp_b_cell_rd{:});

            % store the information about the hyperplanes characterizing the polytopes without the redundant constraints
            obj.Quality.HyperplaneConcise = struct;
            obj.Quality.HyperplaneConcise.w = hp_w_cell;
            obj.Quality.HyperplaneConcise.b = hp_b_cell;

            % prepare quantities for the minimization of quadratic function over the quality space
            obj.Storage.QuadMin = struct;

            convcomb_coef_cell = cell(poly_num, 1);
            convcomb_intercept_cell = cell(poly_num, 1);
            candidate_cell = cell(poly_num, 1);

            for poly_id = 1:poly_num
                % retrieve the vertices in the polytopes and their next counterclockwise neighbor
                vertices_cc1 = obj.Quality.Polytopes{poly_id};
                vertices_cc2 = vertices_cc1([2:end, 1], :);

                % difference of the end points of each edge in x and y coordinates
                vertices_diff_x = vertices_cc2(:, 1) - vertices_cc1(:, 1);
                vertices_diff_y = vertices_cc2(:, 2) - vertices_cc1(:, 2);

                % retrieve the inequalities characterizing the polytope
                hp_w = hp_w_cell{poly_id};
                hp_w_ss = sum(hp_w .^ 2, 2);
                hp_b = hp_b_cell{poly_id};

                % compute the coefficient and the intercept when computing the weight of a point as the convex combination of two end 
                % points of an edge of the polytope
                convcomb_coef_cell{poly_id} = [hp_w(:, 2) .^ 2, -hp_w(:, 1) .* hp_w(:, 2)] ./ hp_w_ss ./ vertices_diff_x;
                convcomb_intercept_cell{poly_id} = (hp_w(:, 1) .* hp_b ./ hp_w_ss - vertices_cc1(:, 1)) ./ vertices_diff_x;

                % if the edge is vertical, the method above will fail, and we thus determine the weight using the y-coordinates
                convcomb_coef_y = [-hp_w(:, 1) .* hp_w(:, 2), hp_w(:, 1) .^ 2] ./ hp_w_ss ./ vertices_diff_y;
                convcomb_intercept_y = (hp_w(:, 2) .* hp_b ./ hp_w_ss - vertices_cc1(:, 2)) ./ vertices_diff_y;

                vertical_list = abs(vertices_diff_x) < 1e-12;
                convcomb_coef_cell{poly_id}(vertical_list, :) = convcomb_coef_y(vertical_list, :);
                convcomb_intercept_cell{poly_id}(vertical_list) = convcomb_intercept_y(vertical_list);

                % the candidate optimizer before scaling
                candidate_cell{poly_id} = hp_w ./ hp_w_ss;
            end

            obj.Storage.QuadMin.IneqWeights = vertcat(hp_w_cell{:});
            obj.Storage.QuadMin.IneqRHS = vertcat(hp_b_cell{:});
            obj.Storage.QuadMin.ConvCombCoefs = vertcat(convcomb_coef_cell{:});
            obj.Storage.QuadMin.ConvCombIntercepts = vertcat(convcomb_intercept_cell{:});
            obj.Storage.QuadMin.UnscaledEdgeCandidates = vertcat(candidate_cell{:});
            obj.Storage.QuadMin.VerticesSS = sum(obj.Quality.Vertices .^ 2, 2);

            % if the test functions are specified, set the test functions
            if exist('quality_testfuncs', 'var') && ~isempty(quality_testfuncs)
                obj.setQualitySimplicialTestFuncs(quality_testfuncs{:});
            end

            % this flag is used to track if the function obj.initializeSimplicialTestFuncs has been called
            obj.Storage.SimplicialTestFuncsInitialized = false;

            marg_num = length(obj.MarginalWeights);
            marg_testfunc_set = false(marg_num, 1);

            for marg_id = 1:marg_num
                marg_testfunc_set(marg_id) = ~isempty(obj.Marginals{marg_id}.SimplicialTestFuncs);
            end

            if isfield(obj.Quality, 'SimplicialTestFuncs') ...
                    && ~isempty(obj.Quality.SimplicialTestFuncs) ...
                    && all(marg_testfunc_set)
                % initialize the simplicial test functions at the end of the constructor if they have already been set
                obj.initializeSimplicialTestFuncs();
            end
        end

        function inside = checkIfInsideQualitySpace(obj, ...
                pts, ...
                batch_size)
            % Check if the points are inside the quality space, which is the union of polytopes. Compute in batches if necessary to
            % avoid excessive memory use
            % Inputs:
            %   pts: two-column matrix containing the input points
            %   batch_size: maximum number of inputs to be evaluated together via a vectorized routine (default is 10300)
            % Output:
            %   inside: boolean vector indicating whether each of the points is inside the quality space

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 10300;
            end

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            inside_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                [inside_cell{batch_id}] = obj.doCheckIfInsideQualitySpace(pts(((batch_id - 1) * batch_size + 1):min(batch_id ...
                    * batch_size, input_num), :));
            end

            inside = vertcat(inside_cell{:});
        end

        function setSimplicialTestFuncs(obj, quality_args, args_cell)
            % Set the simplicial test functions on the quality space as well as the simplicial test functions for all marginals at the 
            % same time
            % Input:
            %   quality_arts: cell array containing all inputs to the method setQualitySimplicialTestFuncs
            %   args_cell: cell array where each cell is a cell array containing all inputs to the method setSimplicialTestFuncs of 
            % each marginal

            obj.setQualitySimplicialTestFuncs(quality_args{:});

            for marg_id = 1:length(obj.MarginalWeights)
                obj.Marginals{marg_id}.setSimplicialTestFuncs(args_cell{marg_id}{:});
            end

            % after setting the simplicial functions, initialize the quantities for the cutting-plane algorithm
            obj.initializeSimplicialTestFuncs();
        end

        function setQualitySimplicialTestFuncs(obj, ...
                vertices, ...
                triangles, ...
                quadratic_testfunc, ...
                check_simpcover)
            % Initialize the test functions with respect to a simplicial cover and compute the integrals of the test functions with
            % respect to the probability measure.
            % Warning: the function does not check if the union of the triangles in the simplicial cover equals the quality space given
            % by the union of polytopes.
            % Inputs:
            %   vertices: two-column matrix containing the vertices used in the triangulation
            %   quadratic_testfunc: boolean indicating the use of a quadratic function on the quality space (default is true)
            %   triangles: three-column matrix containing the triangulation
            %   check_simpcover: check whether the given triangulation forms a simplicial cover (default is true)

            if ~exist('quadratic_testfunc', 'var') || isempty(quadratic_testfunc)
                quadratic_testfunc = true;
            end

            if ~exist('check_simpcover', 'var') || isempty(check_simpcover)
                check_simpcover = true;
            end

            % check if all vertices are contained in the quality space
            assert(all(obj.checkIfInsideQualitySpace(vertices)), 'there are vertices outside the quality space');

            % check if there are duplicate vertices
            assert(size(unique(round(vertices, 6), 'rows'), 1) == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not appear in any triangle)
            assert(length(unique(triangles(:))) == size(vertices, 1), 'there are redundant vertices');

            if check_simpcover
                % check if the triangulation forms a simplicial cover
                [empty, is_simpcover] = check_triangulation_simplicial_cover(vertices, triangles);

                if any(empty)
                    error('some triangles are empty');
                end

                if ~is_simpcover
                    error('the triangulation does not form a simplicial cover');
                end
            end

            obj.Quality.SimplicialTestFuncs = struct;
            obj.Quality.SimplicialTestFuncs.Vertices = vertices;
            obj.Quality.SimplicialTestFuncs.Triangles = triangles;
            obj.Quality.SimplicialTestFuncs.Quadratic = quadratic_testfunc;

            tri_num = size(triangles, 1);

            % compute the inverse transformation matrix which transform a coordinate into weights in the affine combination
            invtrans_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert = vertices(triangles(tri_id, :), :);

                % compute the inverse transformation matrix
                invtrans_cell{tri_id} = [tri_vert'; ones(1, 3)] \ eye(3);
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

            obj.Quality.SimplicialTestFuncs.MeshSize = sqrt(max(max(max(edge1sq, edge2sq), edge3sq)));

            % if a quadratic test function is used, include the centroid of the first triangle as an extra point
            if quadratic_testfunc
                extra_pt = mean(vertices(triangles(1, :)', :), 1)';
                obj.Quality.SimplicialTestFuncs.QuadraticExtraPoint = extra_pt;

                % the quadratic test function evaluates to the quadratic function minus the continuous piece-wise affine interpolation 
                % of the quadratic function; hence, we pre-compute the values of the quadratic function at the vertices
                quad_vals = sum(vertices.^2, 2);
                obj.Quality.SimplicialTestFuncs.QuadraticValues = quad_vals;

                % the value of this test function evaluated at the extra point
                vertex_num = size(vertices, 1);
                eval_val = sum(extra_pt.^2) - mean(quad_vals(triangles(1, :)'));
                obj.Quality.SimplicialTestFuncs.QuadraticExtraPointEval = sparse( ...
                    ones(4, 1), ...
                    [triangles(1, :)'; vertex_num + 1], ...
                    [ones(3, 1) / 3; eval_val], ...
                    1, vertex_num + 1);
            end
        end

        function [vals, inside] = evaluateQualitySimplicialTestFuncs(obj, ...
                pts, ...
                batch_size)
            % Evaluate the test functions on the quality space at given locations; the input locations must be inside the quality space
            % Input:
            %   pts: two-column matrix containing the input points
            %   batch_size: maximum number of inputs to be evaluated together via a vectorized routine (default is 1e4)
            % Output:
            %   vals: sparse matrix containing the computed function values where each row corresponds to an input and each column 
            %   corresponds to a test function
            %   inside: sparse boolean matrix indicating whether an input point is inside each of the triangles; each row corresponds
            %   to an input and each column corresponds to a triangle

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            input_num = size(pts, 1);
            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);
            inside_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                [vals_cell{batch_id}, inside_cell{batch_id}] = obj.doEvaluateQualitySimplicialTestFuncs( ...
                    pts(((batch_id - 1) * batch_size + 1):min(batch_id * batch_size, input_num), :));
            end

            vals = vertcat(vals_cell{:});
            inside = vertcat(inside_cell{:});
        end

        function coup_cell = updateSimplicialTestFuncs(obj, ...
                quality_args, ...
                args_cell, ...
                fixed_index)
            % Update the simplicial test functions after an execution of the cutting-plane algorithm. Besides setting the new 
            % simplicial test functions, a new set of couplings of discretized marginals are generated via reassembly of the dual 
            % solution from the cutting-plane algorithm with the new discretized marginals. These couplings can be used to generate 
            % initial constraints for the new LSIP problem with the updated test functions.
            % Input:
            %   quality_args: cell array containing all inputs to the method setQualitySimplicialTestFuncs
            %   args_cell: cell array where each cell is a cell array containing all inputs to the methods setSimplicialTestFuncs of 
            %   each marginal
            %   fixed_index: index of the marginal whose corresponding measure on the quality space is used in the updated couplings 
            %   (default is 1)
            % Outputs:
            %   coup_cell: cell array containing structs with fields encoding the couplings between a marginal and the measure on the 
            %   quality space

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('fixed_index', 'var') || isempty(fixed_index)
                fixed_index = 1;
            end

            marg_num = length(obj.MarginalWeights);

            % retrieve the dual solution resulted from the cutting-plane algorithm
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
                old_coup_probs_cell{marg_id} = dual_sol{marg_id}.Probabilities;

                % store only the unique atoms in the quality space and use their indices in the coupling with the marginal
                [~, uind, umap] = unique(round(quality_atoms, 6), 'rows', 'stable');
                old_coup_q_atoms_cell{marg_id} = quality_atoms(uind, :);
                old_coup_q_probs_cell{marg_id} = accumarray(umap, old_coup_probs_cell{marg_id});
                old_coup_indices_cell{marg_id} = [dual_sol{marg_id}.VertexIndices, umap];

                if marg_id == fixed_index
                    keep_list = old_coup_q_probs_cell{marg_id} >= 1e-12;
                    fixed_cand_q_probs = old_coup_q_probs_cell{marg_id}(keep_list);
                    fixed_cand_q_probs = fixed_cand_q_probs / sum(fixed_cand_q_probs);
                    fixed_cand_q_atoms = old_coup_q_atoms_cell{marg_id}(keep_list, :);
                end

                old_atoms_cell{marg_id} = obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices;
                old_probs_cell{marg_id} = obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;
            end

            % set the new simplicial test functions
            obj.setSimplicialTestFuncs(quality_args, args_cell);

            % move the atoms in the fixed candidate on the quality space to the closest vertex in the updated test functions; then,
            % take a discrete measure on the quality space that is the mixture of the fixed candidate and a discrete measure that is
            % uniform on all vertices in the new triangulation
            fixed_meas_q_atoms = obj.Quality.SimplicialTestFuncs.Vertices;

            if obj.Quality.SimplicialTestFuncs.Quadratic
                fixed_meas_q_atoms = [fixed_meas_q_atoms; obj.Quality.SimplicialTestFuncs.QuadraticExtraPoint'];
            end

            fixed_meas_interp_const = 0.05;
            fixed_dist_mat = pdist2(fixed_cand_q_atoms, fixed_meas_q_atoms, 'euclidean');
            [~, min_dist_ind] = min(fixed_dist_mat, [], 2);
            fixed_meas_q_atom_num = size(fixed_meas_q_atoms, 1);
            fixed_meas_q_probs = ones(fixed_meas_q_atom_num, 1) / fixed_meas_q_atom_num * fixed_meas_interp_const;
            fixed_meas_q_probs_additional = accumarray(min_dist_ind, fixed_cand_q_probs, [fixed_meas_q_atom_num, 1]);
            fixed_meas_q_probs = fixed_meas_q_probs + fixed_meas_q_probs_additional * (1 - fixed_meas_interp_const);
            fixed_meas_q_probs = fixed_meas_q_probs / sum(fixed_meas_q_probs);

            coup_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                new_atoms = obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices;
                new_probs = obj.Marginals{marg_id}.SimplicialTestFuncs.Integrals;

                % compute an optimal coupling between the original discrete marginal and the new discrete marginal discrete optimal
                % transport

                % the cost function is the Euclidean distance
                dist_mat = pdist2(old_atoms_cell{marg_id}, new_atoms, 'euclidean');
                if size(dist_mat, 1) * size(dist_mat, 2) > 1e6
                    % if there are too many atoms in the discrete measures, a direct computation of discrete OT may cause memory
                    % throttling; thus, we resort to a constraint generation scheme
                    cp_options = struct('display', false);
                    OT = OTDiscrete( ...
                        old_probs_cell{marg_id}, ...
                        new_probs, ...
                        dist_mat, ...
                        cp_options);
                    [hcoup_indices, ~] = OT.generateHeuristicCoupling();
                    OT.run(hcoup_indices, 1e-6);
                    coup = OT.Runtime.DualSolution;
                    marg_coup_atom_indices = coup.CoupIndices;
                    marg_coup_probs = coup.Probabilities;
                else
                    [marg_coup_atom_indices, marg_coup_probs] = discrete_OT( ...
                        old_probs_cell{marg_id}, ...
                        new_probs, ...
                        dist_mat);
                end

                % compute an optimal coupling between the measure on the quality space and the fixed measure on the quality space; this
                % is to ensure that all measures on the quality space satisfy the constraints with respect to the test functions
                q_dist_mat = pdist2(old_coup_q_atoms_cell{marg_id}, fixed_meas_q_atoms, 'euclidean');

                if size(q_dist_mat, 1) * size(q_dist_mat, 2) > 1e6
                    % if there are too many atoms in the discrete measures, a direct computation of discrete OT may cause memory
                    % throttling; thus, we resort to a constraint generation scheme
                    cp_options = struct('display', false);
                    OT = OTDiscrete( ...
                        old_coup_q_probs_cell{marg_id}, ...
                        fixed_meas_q_probs, ...
                        q_dist_mat, ...
                        cp_options);
                    [hcoup_indices, ~] = OT.generateHeuristicCoupling();
                    OT.run(hcoup_indices, 1e-6);
                    coup = OT.Runtime.DualSolution;
                    q_coup_atom_indices = coup.CoupIndices;
                    q_coup_probs = coup.Probabilities;
                else
                    [q_coup_atom_indices, q_coup_probs] = discrete_OT( ...
                        old_coup_q_probs_cell{marg_id}, ...
                        fixed_meas_q_probs, ...
                        q_dist_mat);
                end

                % perform discrete reassembly to get the new coupling; in this case, we need to reassemble both marginals in the 
                % coupling
                [coup_indices, coup_probs] = discrete_reassembly( ...
                    old_coup_indices_cell{marg_id}, ...
                    old_coup_probs_cell{marg_id}, ...
                    {marg_coup_atom_indices; q_coup_atom_indices}, ...
                    {marg_coup_probs; q_coup_probs}, ...
                    1:2);
                coup_points = fixed_meas_q_atoms(coup_indices(:, 2), :);
                coup_marg_vertices = new_atoms(coup_indices(:, 1), :);
                coup_atom_num = length(coup_probs);

                if ~obj.Quality.SimplicialTestFuncs.Quadratic
                    coup_atom_testfunc_vals = sparse( ...
                        1:coup_atom_num, ...
                        coup_indices(:, 2), ...
                        1, ...
                        coup_atom_num, ...
                        fixed_meas_q_atom_num);
                else
                    coup_atom_testfunc_vals = sparse( ...
                        1:coup_atom_num, ...
                        coup_indices(:, 2), ...
                        1, ...
                        coup_atom_num, ...
                        fixed_meas_q_atom_num);

                    % modify the test function values evaluated at the extra point
                    extra_point_list = coup_indices(:, 2) == fixed_meas_q_atom_num;

                    if any(extra_point_list)
                        coup_atom_testfunc_vals(extra_point_list, :) = ...
                            repmat(obj.Quality.SimplicialTestFuncs.QuadraticExtraPointEval, sum(extra_point_list), 1); %#ok<SPRIX>
                    end
                end

                % compute the corresponding cost function values
                coup_costfunc_vals = obj.MarginalWeights(marg_id) * sum(coup_points .* (coup_points - 2 * coup_marg_vertices), 2);

                coup_cell{marg_id} = struct( ...
                    'probabilities', coup_probs, ...
                    'vertex_indices', coup_indices(:, 1), ...
                    'points', coup_points, ...
                    'testfuncs_vals', coup_atom_testfunc_vals, ...
                    'costfunc_vals', coup_costfunc_vals);
            end
        end

        function initializeSimplicialTestFuncs(obj)
            % Initialize some quantities related to the simplicial test functions of the marginals and the simplicial test functions on
            % the quality space

            marg_num = length(obj.MarginalWeights);

            quality_testfuncs = obj.Quality.SimplicialTestFuncs;
            quality_vert_num = size(quality_testfuncs.Vertices, 1);
            tri_num = size(quality_testfuncs.Triangles, 1);
            obj.Storage.QualityTriNum = tri_num;

            % store the number of vertices in the test functions for each marginal
            obj.Storage.MargVertNumList = zeros(marg_num, 1);

            % store the indices of the decision variables that correspond to the intercepts, the coefficients of the test functions for
            % the marginals, and the coefficients of the test functions on the quality space
            obj.Storage.DeciVarInterceptIndices = zeros(marg_num, 1);
            obj.Storage.DeciVarMargTestFuncIndices = cell(marg_num, 1);
            obj.Storage.DeciVarQualityTestFuncIndices = cell(marg_num, 1);

            if quality_testfuncs.Quadratic
                obj.Storage.DeciVarQualityQuadIndices = zeros(marg_num, 1);
            else
                obj.Storage.DeciVarQualityQuadIndices = [];
            end

            ind_counter = 0;

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                obj.Storage.MargVertNumList(marg_id) = size(marg.SimplicialTestFuncs.Vertices, 1);

                if ~obj.Options.sparse_parametrization
                    % if the constant intercept is present, the first test function for the marginal is removed for identification
                    % purposes
                    obj.Storage.DeciVarInterceptIndices(marg_id) = ind_counter + 1;
                    ind_counter = ind_counter + 1;
                    obj.Storage.DeciVarMargTestFuncIndices{marg_id} = ind_counter + (1:obj.Storage.MargVertNumList(marg_id) - 1)';
                    ind_counter = ind_counter + obj.Storage.MargVertNumList(marg_id) - 1;
                else
                    % otherwise, use the full set of thet test functions and omit the constant intercept
                    obj.Storage.DeciVarMargTestFuncIndices{marg_id} = ind_counter + (1:obj.Storage.MargVertNumList(marg_id))';
                    ind_counter = ind_counter + obj.Storage.MargVertNumList(marg_id);
                end

                % note that the first test function on the quality space is removed for identification purposes
                obj.Storage.DeciVarQualityTestFuncIndices{marg_id} = ind_counter + (1:quality_vert_num - 1)';
                ind_counter = ind_counter + quality_vert_num - 1;

                if quality_testfuncs.Quadratic
                    obj.Storage.DeciVarQualityQuadIndices(marg_id) = ind_counter + 1;
                    ind_counter = ind_counter + 1;
                end
            end

            obj.Storage.DeciVarLength = ind_counter;
            obj.Storage.TotalMargVertNum = sum(obj.Storage.MargVertNumList);

            % struct to store information for the global minimization oracle
            GM = struct;

            % when solving the global minimization problem, split the atoms in the type space into batches to avoid excessive memory
            % usage
            if tri_num <= 100
                marg_vert_batch_size = 1e4;
            elseif tri_num <= 200
                marg_vert_batch_size = 5e3;
            elseif tri_num <= 500
                marg_vert_batch_size = 2e3;
            elseif tri_num <= 1000
                marg_vert_batch_size = 1e3;
            elseif tri_num <= 2000
                marg_vert_batch_size = 5e2;
            elseif tri_num <= 5000
                marg_vert_batch_size = 2e2;
            else
                marg_vert_batch_size = 1e2;
            end

            GM.Marginals = cell(marg_num, 1);

            for marg_id = 1:marg_num
                GM_marg = struct;
                marg_vert_num = obj.Storage.MargVertNumList(marg_id);

                % compute the number of batches and the indices of the vertices assigned to each batch
                GM_marg.BatchNum = ceil(marg_vert_num / marg_vert_batch_size);
                GM_marg.BatchLists = cell(GM_marg.BatchNum, 1);

                for batch_id = 1:GM_marg.BatchNum
                    batch_start = (batch_id - 1) * marg_vert_batch_size + 1;
                    batch_end = min(batch_id * marg_vert_batch_size, marg_vert_num);
                    GM_marg.BatchLists{batch_id} = (batch_start:batch_end)';
                end

                GM.Marginals{marg_id} = GM_marg;
            end

            % consider three cases for the location of the optimizer: face, edge, and vertex

            % quantities related to the face case
            % there are three parts involved in the computation of the two coefficients: a fixed part, a part that depends on the
            % vertex on the type space (known in advance but computed at runtime to avoid excessive memory usage), and a part that 
            % depends on the coefficients of the test functions on the quality space; the two coefficients are stored separately
            face_num = tri_num;
            face_fixed = cell(2, 1);
            face_fixed{1} = zeros(face_num, 1);
            face_fixed{2} = zeros(face_num, 1);
            face_marg = cell(2, 1);
            face_marg{1} = zeros(face_num, 2);
            face_marg{2} = zeros(face_num, 2);
            face_tf = cell(2, 1);
            face_tf1 = cell(face_num, 1);
            face_tf2 = cell(face_num, 1);

            for face_id = 1:face_num
                tri_indices = quality_testfuncs.Triangles(face_id, :)';

                % the coordinates of the third vertex
                z3 = quality_testfuncs.Vertices(tri_indices(3), :)';

                % 2-by-2 matrix formed by the first two vertices minus the third vertex
                Z = quality_testfuncs.Vertices(tri_indices(1:2), :) - z3';

                face_fixed_full = -(Z * Z') \ (Z * z3); % 2-by-1 vector
                face_fixed{1}(face_id) = face_fixed_full(1);
                face_fixed{2}(face_id) = face_fixed_full(2);

                face_marg_full = (Z * Z') \ Z; % 2-by-2 matrix
                face_marg{1}(face_id, :) = face_marg_full(1, :);
                face_marg{2}(face_id, :) = face_marg_full(2, :);

                face_tf_full = sparse(Z * Z') \ sparse( ...
                    [1; 1; 2; 2], ...
                    [tri_indices(1); tri_indices(3); tri_indices(2); tri_indices(3)], ...
                    [1; -1; 1; -1], ...
                    2, quality_vert_num); % 2-by-quality_vert_num sparse matrix
                face_tf1{face_id} = face_tf_full(1, :);
                face_tf2{face_id} = face_tf_full(2, :);
            end

            face_tf{1} = vertcat(face_tf1{:});
            face_tf{2} = vertcat(face_tf2{:});

            GM.CaseFace = struct;
            GM.CaseFace.Fixed = face_fixed;
            GM.CaseFace.Marginal = face_marg;
            GM.CaseFace.TestFunc = face_tf;


            % quantities related to the edge case
            % similar to the face case, there are three parts involved in the computation of the coefficient: a fixed part, a part that
            % depends on the vertex on the type space, and a part that depends on the coefficients of the test functions on the quality
            % space

            % first, list the edges
            edge_list = zeros(tri_num * 3, 2);
            edge_counter = 0;

            for tri_id = 1:tri_num
                % sort the indices into ascending order to avoid counting duplicates
                tri_indices = sort(quality_testfuncs.Triangles(tri_id, :), 'ascend');

                if ~ismember(tri_indices([1, 2]), edge_list, 'rows')
                    edge_counter = edge_counter + 1;
                    edge_list(edge_counter, :) = tri_indices([1, 2]);
                end

                if ~ismember(tri_indices([1, 3]), edge_list, 'rows')
                    edge_counter = edge_counter + 1;
                    edge_list(edge_counter, :) = tri_indices([1, 3]);
                end

                if ~ismember(tri_indices([2, 3]), edge_list, 'rows')
                    edge_counter = edge_counter + 1;
                    edge_list(edge_counter, :) = tri_indices([2, 3]);
                end
            end

            edge_list = edge_list(1:edge_counter, :);
            edge_num = edge_counter;

            edge_fixed = zeros(edge_num, 1);
            edge_marg = zeros(edge_num, 2);
            edge_tf = cell(edge_num, 1);

            for edge_id = 1:edge_num
                z1 = quality_testfuncs.Vertices(edge_list(edge_id, 1), :)';
                z2 = quality_testfuncs.Vertices(edge_list(edge_id, 2), :)';
                inv_normsq = 1 / sum((z1 - z2).^2);

                edge_fixed(edge_id) = -inv_normsq * (z1 - z2)' * z2;
                edge_marg(edge_id, :) = inv_normsq * (z1 - z2)';
                edge_tf{edge_id} = sparse( ...
                    [1; 1], ...
                    [edge_list(edge_id, 1); edge_list(edge_id, 2)], ...
                    [1; -1] * inv_normsq, ...
                    1, quality_vert_num);
            end

            GM.CaseEdge = struct('Edges', edge_list, ...
                'Fixed', edge_fixed, ...
                'Marginal', edge_marg, ...
                'TestFunc', vertcat(edge_tf{:}));

            % in the vertex case the cost can be directly evaluated, no pre-computation is needed

            obj.Storage.GlobalMin = GM;
            obj.Storage.SimplicialTestFuncsInitialized = true;

            % updating the simplicial test functions will invalidate all quantities in the runtime environment, thus all variables in
            % the runtime environment need to be flushed
            obj.Runtime = [];
        end

        function coup_cell = generateDiscreteCoupling(obj)
            % Generate a feasible dual solution by solving the discretized version of the problem
            % Output:
            %   coup_cell: cell array containing structs with fields encoding the couplings between a marginal and the measure on the 
            %   quality space

            marg_num = length(obj.MarginalWeights);
            
            if ~obj.Quality.SimplicialTestFuncs.Quadratic
                quality_vertices = obj.Quality.SimplicialTestFuncs.Vertices;
            else
                quality_vertices = [obj.Quality.SimplicialTestFuncs.Vertices; obj.Quality.SimplicialTestFuncs.QuadraticExtraPoint'];
            end

            quality_vert_num = size(quality_vertices, 1);
            marg_vert_num_list = obj.Storage.MargVertNumList;
            decivar_num = sum(marg_vert_num_list) + (quality_vert_num - 1) * marg_num;

            % we build the dual version of the problem (maximization) since it is easier to code
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
                
                % the integrals are placed into the objective vector; note that the first test function for the marginal is removed for
                % identification
                objective_cell{marg_id} = [1; marg_integral(2:end); zeros(quality_vert_num - 1, 1)];

                % matrix for the equality constraints requiring that the transfer functions must sum up to 0; note that the first test 
                % function on the quality space is removed for identification
                A_eq_cell{marg_id} = [sparse(quality_vert_num - 1, marg_vert_num), speye(quality_vert_num - 1)];

                % matrix for the inequality constraints will have all possible combinations of the vertices in the test functions for 
                % the marginal and the vertices in the test functions on the quality space
                [marg_grid, q_grid] = meshgrid(1:marg_vert_num, 1:quality_vert_num);
                marg_grid_indices = marg_grid(:);
                q_grid_indices = q_grid(:);
                A_marg_full = sparse( ...
                    1:marg_vert_num * quality_vert_num, ...
                    marg_grid_indices, ...
                    1, ...
                    marg_vert_num * quality_vert_num, marg_vert_num);
                A_quality_full = sparse( ...
                    1:marg_vert_num * quality_vert_num, ...
                    q_grid_indices, ...
                    1, ...
                    marg_vert_num * quality_vert_num, quality_vert_num);

                % remove the first test function for the marginal and the first test function on the quality space; also prepend a
                % column of 1s
                A_ineq_cell{marg_id} = [sparse(ones(marg_vert_num * quality_vert_num, 1)), ...
                    A_marg_full(:, 2:end), A_quality_full(:, 2:end)];

                % the right-hand side of the inequality constraints are computed from the corresponding coordinates of the vertices
                marg_grid_pts = marg_vertices(marg_grid_indices, :);
                quality_grid_pts = quality_vertices(q_grid_indices, :);
                rhs_ineq_cell{marg_id} = obj.MarginalWeights(marg_id) * sum(quality_grid_pts ...
                    .* (quality_grid_pts - 2 * marg_grid_pts), 2);

                constr_num_list(marg_id) = marg_vert_num * quality_vert_num;
                marg_indices_cell{marg_id} = marg_grid_indices;
                quality_indices_cell{marg_id} = q_grid_indices;
                quality_points_cell{marg_id} = quality_grid_pts;
            end
            
            % build a LP model in gurobi
            model = struct;
            model.modelsense = 'max';
            model.objcon = 0;

            % the coefficients corresponding to the first test function of each marginal is not included in the decision variables for
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
            model.sense = [repmat('=', length(rhs_eq), 1); repmat('<', length(model.rhs) - length(rhs_eq), 1)];

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
                coup_probs = probs(constr_counter + (1:constr_num_list(marg_id)));

                % retain only those atoms with positive probabilities
                pos_list = coup_probs > 0;
                coup_probs = coup_probs(pos_list);
                coup_atom_num = length(coup_probs);
                coup_marg_indices = marg_indices_cell{marg_id}(pos_list, :);
                coup_quality_indices = quality_indices_cell{marg_id}(pos_list);
                coup_quality_points = quality_points_cell{marg_id}(pos_list, :);
                coup_testfunc_vals = sparse(1:coup_atom_num, coup_quality_indices, 1, coup_atom_num, quality_vert_num);

                if obj.Quality.SimplicialTestFuncs.Quadratic
                    extra_point_list = coup_quality_indices == quality_vert_num;

                    if any(extra_point_list)
                        coup_testfunc_vals(extra_point_list, :) = ...
                            repmat(obj.Quality.SimplicialTestFuncs.QuadraticExtraPointEval, sum(extra_point_list), 1); %#ok<SPRIX>
                    end
                end

                coup_cost = rhs_ineq_cell{marg_id}(pos_list, :);

                coup_cell{marg_id} = struct( ...
                    'probabilities', coup_probs, ...
                    'vertex_indices', coup_marg_indices, ...
                    'points', coup_quality_points, ...
                    'testfuncs_vals', coup_testfunc_vals, ...
                    'costfunc_vals', coup_cost);

                % update the constraint counter
                constr_counter = constr_counter + constr_num_list(marg_id);
            end
        end

        function setLSIPSolutions(obj, ...
                primal_sol, ...
                dual_sol, ...
                LSIP_UB, ...
                LSIP_LB)
            % Set the primal and dual solution of the LSIP problem in the runtime environment. This is to allow certain quantities to
            % be computed using stored versions of the primal and dual solutions without executing the cutting-plane algorithm again. 
            % Inputs: 
            %   primal_sol: struct containing information about the primal solution of the LSIP problem
            %   dual_sol: struct containing information about the dual solution of the LSIP problem
            %   LSIP_UB: the upper bound for the optimal value of the LSIP problem
            %   LSIP_LB: the lower bound for the optimal value of the LSIP problem

            if isempty(obj.Runtime)
                obj.Runtime = struct;
            end

            obj.Runtime.PrimalSolution = primal_sol;
            obj.Runtime.DualSolution = dual_sol;
            obj.Runtime.LSIP_UB = LSIP_UB;
            obj.Runtime.LSIP_LB = LSIP_LB;
        end

        function OT_info_cell = performReassembly(obj)
            % Perform reassembly by computing semi-discrete optimal transport
            % Output:
            %   OT_info_cell: cell array where each cell is a cell array containing optimal transport-related information that can be 
            %   saved and loaded later

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- semi-discrete OT starts ---\n');
            end
            
            marg_num = length(obj.MarginalWeights);

            if ~isfield(obj.Storage, 'OT') || isempty(obj.Storage.OT)
                obj.Storage.OT = struct;
            end

            obj.Storage.OT.filtered_atom_list = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg_angle_num = obj.Options.OT.angle_num{marg_id};
                marg_pp_angle_num = obj.Options.OT.pp_angle_num{marg_id};
                marg_normalize_gradient = obj.Options.OT.normalize_gradient{marg_id};

                marg_options = obj.Options.OT.optimization_options{marg_id};

                % the atoms and the corresponding probabilities of the discretized marginal are exactly given by the test functions and
                % their respective integrals; this is only valid due to the quadratic structure of the cost function
                marg = obj.Marginals{marg_id};
                marg_atoms = marg.SimplicialTestFuncs.Vertices;
                marg_probs = marg.SimplicialTestFuncs.Integrals;

                % filter out the atoms with probabilities below the given threshold
                negligible_indices = marg_probs < obj.Options.OT.prob_thres{marg_id};

                % compute the pairwise distance between the negligible atoms and the remaining atoms and then compute the nearest
                % neighbor of each negligible atom
                dist_mat = pdist2(marg_atoms(negligible_indices, :), marg_atoms(~negligible_indices, :));
                [~, nearest_neighbor_indices] = min(dist_mat, [], 2);
                
                atom_list = (1:size(marg_atoms, 1))';
                remaining_atom_list = atom_list(~negligible_indices);
                atom_list(negligible_indices) = remaining_atom_list(nearest_neighbor_indices);

                filtered_atoms = marg_atoms(remaining_atom_list, :);
                [~, ~, atom_list] = unique(atom_list, 'sorted');

                filtered_probs = accumarray(atom_list, marg_probs, size(remaining_atom_list));

                marg.computeOptimalTransport( ...
                    filtered_atoms, ...
                    filtered_probs, ...
                    marg_angle_num, ...
                    [], ...
                    marg_pp_angle_num, ...
                    marg_options, ...
                    [], ...
                    marg_normalize_gradient);

                if obj.Options.display
                    fprintf('%s: marginal %d done\n', class(obj), marg_id);
                end

                % logging
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, '%s: marginal %d done\n', class(obj), marg_id);
                end

                obj.Storage.OT.filtered_atom_list{marg_id} = atom_list;
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- semi-discrete OT ends ---\n\n');
                fclose(log_file);
            end

            obj.Storage.OTComputed = true;

            OT_info_cell = obj.saveOptimalTransportInfo();
        end

        function OT_info_cell = saveOptimalTransportInfo(obj)
            % Retrieve computed semi-discrete optimal transport-related information of each marginal
            % Output:
            %   OT_info_cell: cell array where each cell is a cell array containing optimal transport-related information that can be 
            %   saved and loaded later

            marg_num = length(obj.MarginalWeights);

            OT_info_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                OT_info_cell{marg_id} = marg.saveOptimalTransportInfo();

                if isfield(obj.Storage, 'OT') && isfield(obj.Storage.OT, 'filtered_atom_list')
                    OT_info_cell{marg_id}.filtered_atom_list = obj.Storage.OT.filtered_atom_list{marg_id};
                end
            end
        end

        function loadOptimalTransportInfo(obj, OT_info_cell)
            % Load semi-discrete optimal transport-related information into each marginal
            % Input:
            %   OT_info_cell: cell array where each cell is a cell array containing optimal transport-related information that can be 
            %   saved and loaded later

            marg_num = length(obj.MarginalWeights);

            if ~isfield(obj.Storage, 'OT')
                obj.Storage.OT = struct;
            end

            for marg_id = 1:marg_num
                if isfield(OT_info_cell{marg_id}, 'filtered_atom_list')
                    if ~isfield(obj.Storage.OT, 'filtered_atom_list')
                        obj.Storage.OT.filtered_atom_list = cell(marg_num, 1);
                    end

                    obj.Storage.OT.filtered_atom_list{marg_id} = OT_info_cell{marg_id}.filtered_atom_list;

                    OT_info_cell{marg_id} = rmfield(OT_info_cell{marg_id}, 'filtered_atom_list');
                end

                marg = obj.Marginals{marg_id};
                marg.loadOptimalTransportInfo(OT_info_cell{marg_id});
            end

            obj.Storage.OTComputed = true;
        end

        function vals = evaluateOptParametricFunc(obj, ...
                marg_id, ...
                pts, ...
                batch_size)
            % Evaluate the parametric function on the support of a marginal parametrized by simplicial test functions that is optimized
            % using the cutting plane algorithm. This corresponds to the primal solution of the LSIP problem. The computation is done
            % in batches if necessary.
            % Inputs:
            %   marg_id: the index of the marginal
            %   pts: two-column matrix containing the input points
            %   batch_size: the maximum number of input points to be handled at the same time in the vectorized procedure (default is 
            %   1e4)
            % Output:
            %   vals: vector containing the computed function values

            if ~isfield(obj.Runtime, 'PrimalSolution') || isempty(obj.Runtime.PrimalSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            marg = obj.Marginals{marg_id};
            primal_sol = obj.Runtime.PrimalSolution{marg_id};
            vals = marg.evaluateWeightedSumOfSimplicialTestFuncs(pts, primal_sol.Coefficients, batch_size) + primal_sol.Constant;
        end

        function vals_mat = evaluateOptTransferFuncs(obj, ...
                pts, ...
                ref_pt, ...
                batch_size)
            % Evaluate the transfer functions resulted from optimized parametric functions. These transfer functions is part of
            % approximate matching equilibria. The computation is done in batches if necessary.
            % Inputs:
            %   pts: two-column matrix containing the input points
            %   ref_pt: two-element vector indicating a reference point where the transfer function will evaluate to 0 (default is the 
            %   first vertex characterizing the quality space)
            %   batch_size: the maximum number of input points to be handled at the same time in the vectorized procedure (default is 
            %   10001)
            % Output:
            %   vals_mat: matrix where each column contains the computed transfer function corresponding to a marginal at the input
            %   points

            if ~isfield(obj.Runtime, 'PrimalSolution') || isempty(obj.Runtime.PrimalSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 10001;
            end

            if ~exist('ref_pt', 'var') || isempty(ref_pt)
                ref_pt = obj.Quality.Vertices(1, :)';
            else
                assert(obj.checkIfInsideQualitySpace(ref_pt'), 'the reference point is not in the quality space');
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
                    cur_pts = pts(((batch_id - 1) * batch_size + 1):min(batch_id * batch_size, pt_num), :);

                    CPWA_val_mat = (2 * marg_weight) * cur_pts * marg_vertices' + marg_coefficients';
                    vals_cell{batch_id} = sum(cur_pts .^ 2, 2) .* marg_weight - marg_constant - max(CPWA_val_mat, [], 2);

                    if batch_id == 1
                        % store the values of the transfer functions before shifting evaluated at the reference point
                        ref_vals = vals_cell{batch_id}(1);

                        % shift the transfer functions to make the values vanish at the reference point, then remove the reference 
                        % point
                        vals_cell{batch_id} = vals_cell{batch_id}(2:end) - ref_vals;
                    else
                        % shift the transfer functions to make the values vanish at the reference point
                        vals_cell{batch_id} = vals_cell{batch_id} - ref_vals;
                    end
                end

                vals_mat(:, marg_id) = vertcat(vals_cell{:});
            end

            % the last transfer function is obtained from the balance condition, the resulting transfer functions are guaranteed to add
            % up to 0
            vals_mat(:, end) = -sum(vals_mat(:, 1:end - 1), 2);
        end

        function LB = getMTLowerBound(obj)
            % Retrieve the computed lower bound for the matching for teams problem
            % Output:
            %   LB: the computed lower bound

            if ~isfield(obj.Runtime, 'PrimalSolution') || isempty(obj.Runtime.PrimalSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            LB = -obj.Runtime.LSIP_UB;
        end

        function distri = getOptDiscreteQualityDistr(obj)
            % Retrieve the optimized discrete quality distribution. This is a part of an approximate matching equilibrium.
            % Output:
            %   distri: struct containing fields Probabilities and Atoms

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            candidate = obj.Runtime.DualSolution{obj.Options.discmeas_cand_index};

            distri = struct;
            distri.Probabilities = candidate.Probabilities;
            distri.Atoms = candidate.QualityAtoms;
        end

        function samps = randSampleFromOptJointDistr(obj, ...
                samp_num, ...
                rand_stream, ...
                batch_size)
            % Generate independent random sample from the joint distribution consisting of the candidate discrete measure on the 
            % quality space, the corresponding coupling of the discretized marginals resulted from binding, the coupled continuous 
            % marginals, and the corresponding continuous measure in the quality space
            % Inputs:
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the continuous distribution on the quality space (default is
            %   1e4)
            % Output:
            %   samps: struct containing the following fields:
            %       DiscreteIndices: matrix containing the indices of the atoms in the discretimized marginals; each row represents a 
            %       sample and each column represents a marginal
            %       DiscretePoints: cell array containing the coordinates of samples from the coupling of the discretized marginals
            %       DiscreteQualities: two-column matrix containing the coordinates of samples from the discrete measure on the
            %       quality space
            %       ContinuousPoints: cell array containing the coordinates of samples from the continuous marginals
            %       ContinuousQualities: two-column matrix containing the coordinates of samples from the continuous measure on the 
            %       quality space

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            marg_num = length(obj.MarginalWeights);
            samps = obj.doRandSampleFromPartialReassembly((1:marg_num)', samp_num, rand_stream);

            % first, compute the weights in the quadratic minimization problem
            weights = zeros(samp_num, 2);

            for marg_id = 1:marg_num
                weights = weights + obj.MarginalWeights(marg_id) * samps.ContinuousPoints{marg_id};
            end

            % then, the samples from the continuous measure on the quality space are obtained by solving the quadratic minimization
            % problem over the quality space
            samps.ContinuousQualities = obj.computeQuadMin(weights, batch_size);
        end

        function coup = randSampleFromOptDiscreteCoupling(obj, ...
                marg_id, ...
                samp_num, ...
                rand_stream)
            % Generate independent random sample from the coupling between the discrete measure on the quality space and a continuous
            % marginal
            % Inputs:
            %   marg_id: index of the marginal
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream)
            % Output:
            %   coup: struct containing the following fields:
            %       Qualities: two-column matrix containing the coordinates of samples from the discrete measure on the quality space
            %       AgentTypes: two-column matrix containing the coordinates of samples from the continuous marginal

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            samps = obj.doRandSampleFromPartialReassembly(marg_id, samp_num, rand_stream);

            coup = struct;
            coup.AgentTypes = samps.ContinuousPoints{marg_id};

            % get the corresponding samples from the discrete measure on the quality space
            coup.Qualities = samps.DiscreteQualities;
        end

        function ql_samps = randSampleFromOptContinuousQualityDistr(obj, ...
                samp_num, ...
                rand_stream, ...
                batch_size)
            % Generate independent random sample from the continuous measure on the quality space
            % Inputs:
            %   marg_id: index of the marginal
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the continuous distribution on the quality space (default is
            %   1e4)
            % Output:
            %   ql_samps: two-column matrix containing the coordinates of samples from the discrete measure on the quality space

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            samps = obj.randSampleFromOptJointDistr(samp_num, rand_stream, batch_size);
            ql_samps = samps.ContinuousQualities;
        end

        function coup = randSampleFromOptContinuousCoupling(obj, ...
                marg_id, ...
                samp_num, ...
                rand_stream, ...
                batch_size)
            % Generate independent random sample from the coupling between the continuous measure on the quality space and a continuous
            % marginal
            % Inputs:
            %   marg_id: index of the marginal
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the continuous distribution on the quality space (default is
            %   1e4)
            % Output:
            %   coup: struct containing the following fields:
            %       Qualities: two-column matrix containing the coordinates of samples from the continuous measure on the quality space
            %       AgentTypes: two-column matrix containing the coordinates of samples from the continuous marginal

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            samps = obj.randSampleFromOptJointDistr(samp_num, rand_stream, batch_size);
            coup = struct;
            coup.Qualities = samps.ContinuousQualities;
            coup.AgentTypes = samps.ContinuousPoints{marg_id};
        end

        function [UB_disc, UB_cont, samps] = getMTUpperBounds(obj, ...
                samp_num, rand_stream, batch_size)
            % Compute two upper bounds for the matching for teams problem based on the discrete measure on the quality space and based
            % on the continuous measure on the quality space. The bounds are approximated by Monte Carlo integration.
            % Inputs:
            %   samp_num: number of samples for Monte Carlo integration
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the continuous distribution on the quality space (default is
            %   1e4)
            % Output:
            %   UB_disc: upper bound for the matching for teams problem based on the discrete measure on the quality space
            %   UB_cont: upper bound for the matching for teams problem based on the continuous measure on the quality space
            %   samps: the Monte Carlo samples used in the approximation of the bounds

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            marg_num = length(obj.MarginalWeights);

            samps = obj.randSampleFromOptJointDistr(samp_num, rand_stream, batch_size);

            UB_disc_mat = zeros(samp_num, marg_num);
            UB_cont_mat = zeros(samp_num, marg_num);

            for marg_id = 1:marg_num
                scaled_marg_samps = 2 * samps.ContinuousPoints{marg_id};
                ql_disc = samps.DiscreteQualities;
                ql_cont = samps.ContinuousQualities;

                % compute each term that contributes to the upper bounds via Monte Carlo integration
                UB_disc_mat(:, marg_id) = sum((ql_disc - scaled_marg_samps) .* ql_disc, 2);
                UB_cont_mat(:, marg_id) = sum((ql_cont - scaled_marg_samps) .* ql_cont, 2);
            end

            UB_disc_list = UB_disc_mat * obj.MarginalWeights + obj.QuadraticConstant;
            UB_cont_list = UB_cont_mat * obj.MarginalWeights + obj.QuadraticConstant;

            UB_disc = mean(UB_disc_list);
            UB_cont = mean(UB_cont_list);

            samps.UBDiscrete = UB_disc_list;
            samps.UBContinuous = UB_cont_list;
        end

        function [UB_disc_list, ...
                UB_cont_list, ...
                samps, ...
                samp_histpdf] ...
                = getMTUpperBoundsWRepetition(obj, ...
                samp_num, ...
                rep_num, ...
                rand_stream, ...
                batch_size, ...
                hist_edge_x, ...
                hist_edge_y)
            % Compute two upper bounds for the matching for teams problem with repetition of Monte Carlo integration.
            % Inputs:
            %   samp_num: number of samples for Monte Carlo integration
            %   rep_num: number of repetitions
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream)
            %   batch_size: batch size when computing the samples from the continuous distribution on the quality space (default is
            %   1e4)
            %   hist_edge_x: edges of bins on the x-axis for 2D pdf estimation (default is [])
            %   hist_edge_y: edges of bins on the y-axis for 2D pdf estimation (default is [])
            % Output:
            %   UB_disc_list: vector containing the computed upper bounds for the matching for teams problem based on the discrete 
            %   measure on the quality space
            %   UB_cont_list: vector containing the computed upper bounds for the matching for teams problem based on the continuous 
            %   measure on the quality space
            %   samps: one particular set of Monte Carlo samples used in the approximation of the bounds
            %   samp_histpdf: 2D pdf estimation via the histogram of the generated samples from the continuous quality distribution;
            %   only computed if both hist_edge_x and hist_edge_y are set

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('batch_size', 'var') || isempty(batch_size)
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

                fprintf(log_file, '--- Monte Carlo sampling starts ---\n');
            end

            % initialize the 2D histogram if it needs to be computed
            if ~isempty(hist_edge_x) && ~isempty(hist_edge_y)
                pdf_computed = true;

                samp_histpdf = zeros(length(hist_edge_x) - 1, length(hist_edge_y) - 1);
            else
                pdf_computed = false;
            end

            for rep_id = 1:rep_num
                [UB_disc_list(rep_id), ...
                    UB_cont_list(rep_id), ...
                    samps] ...
                    = obj.getMTUpperBounds( ...
                    samp_num, ...
                    rand_stream, ...
                    batch_size);

                % display output
                if obj.Options.display
                    fprintf('%s: Monte Carlo sampling repetition %3d done\n', class(obj), rep_id);
                end

                % write log
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, '%s: Monte Carlo sampling repetition %3d done\n', class(obj), rep_id);
                end

                % update the 2D histogram
                if pdf_computed
                    samp_histpdf = samp_histpdf + histcounts2(samps.ContinuousQualities(:, 1), samps.ContinuousQualities(:, 2), ...
                        hist_edge_x, hist_edge_y, 'Normalization', 'pdf');
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- Monte Carlo sampling ends ---\n\n');
                fclose(log_file);
            end

            % normalize the 2D histogram
            if pdf_computed
                samp_histpdf = samp_histpdf / rep_num;
            end
        end

        function [W2OT, W2OT_info_cell] = computeW2OptimalCouplings(obj)
            % Compute the Wasserstein-2 optimal couplings between the computed the discrete measure given by the LSIP dual solution and
            % the marginals.
            % Outputs:
            %   W2OT: struct containing Wasserstein-2 optimal transport-related information that can be saved and loaded later 
            %   W2OT_info_cell: cell array where each cell is a cell array containing Wasserstein-2 optimal transport-related
            %   information about a marginal that can be saved and loaded later

            if ~obj.checkIfWasserstein2OTSupported()
                error('Wasserstein-2 optimal transport is not supported by the marginals');
            end

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- semi-discrete Wasserstein-2 OT starts ---\n');
            end
            
            marg_num = length(obj.MarginalWeights);

            % retrieve the discrete measure from the dual LSIP solution
            discrete_meas = obj.getOptDiscreteQualityDistr();
            disc_probs = discrete_meas.Probabilities;
            disc_atoms = discrete_meas.Atoms;

            % some atoms might be very close to each other; they will be combined
            [~, uind, umap] = unique(round(disc_atoms, 6), 'rows', 'stable');
            disc_atoms = disc_atoms(uind, :);
            disc_atom_num = size(disc_atoms, 1);
            disc_probs = accumarray(umap, disc_probs, [disc_atom_num, 1]);

            % atoms with probabilities that are too small are removed to avoid numerical issues
            small_prob_list = disc_probs < obj.Options.dual_prob_thres;
            disc_atoms = disc_atoms(~small_prob_list, :);
            disc_probs = disc_probs(~small_prob_list);
            disc_probs = disc_probs / sum(disc_probs);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_options = obj.Options.W2OT.optimization_options{marg_id};

                % compute the Wasserstein-2 optimal transport to the discrete measure
                marg.computeW2OptimalTransport(disc_atoms, disc_probs, [], marg_options);

                if obj.Options.display
                    fprintf('%s: marginal %d done\n', class(obj), marg_id);
                end

                % logging
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, '%s: marginal %d done\n', class(obj), marg_id);
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- semi-discrete Wasserstein-2 OT ends ---\n\n');
                fclose(log_file);
            end

            obj.Runtime.W2OTComputed = true;
            obj.Runtime.W2OT = struct;
            obj.Runtime.W2OT.DiscMeas = struct('Atoms', disc_atoms, 'Probabilities', disc_probs);

            [W2OT, W2OT_info_cell] = obj.saveW2OptimalTransportInfo();
        end

        function [objective, W2_distances, full_info_cell] = evaluateObjectiveViaW2(obj, atoms, probs, optimization_options)
            % Compute the Wasserstein-2 distances from the given discrete measure to the marginals in order to evaluate its objective 
            % with respect to the Wasserstein barycenter problem
            % Inputs:
            %   atoms: two-column matrix containing the atoms
            %   probs: vector containing the probabilities of the atoms
            %   optimization_options: struct to cell array of structs containing the parameters for the MATLAB fminunc solver
            % Outputs:
            %   objective: the evaluated objective value
            %   W2_distances: vector containing the Wasserstein-2 distances from the given discrete measure to each of the marginals
            %   full_info_cell: the Wasserstein-2 OT related info from the computation with respect to each marginal

            if ~obj.checkIfWasserstein2OTSupported()
                error('Wasserstein-2 optimal transport is not supported by the marginals');
            end

            if ~exist('optimization_options', 'var') || isempty(optimization_options)
                optimization_options = cell(length(obj.MarginalWeights), 1);
            elseif ~iscell(optimization_options)
                optim_options = optimization_options;
                optimization_options = cell(length(obj.MarginalWeights), 1);
                optimization_options(:, :) = {optim_options};
            end

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- semi-discrete Wasserstein-2 OT starts ---\n');
            end
            
            marg_num = length(obj.MarginalWeights);
            W2_distances = zeros(marg_num, 1);
            full_info_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_options = optimization_options{marg_id};

                % compute the Wasserstein-2 optimal transport to the discrete measure
                marg.computeW2OptimalTransport(atoms, probs, [], marg_options);

                if obj.Options.display
                    fprintf('%s: marginal %d done\n', class(obj), marg_id);
                end

                % logging
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, '%s: marginal %d done\n', class(obj), marg_id);
                end

                full_info_cell{marg_id} = marg.W2OT;
                W2_distances(marg_id) = marg.W2OT.Cost;
            end
            
            objective = obj.MarginalWeights' * W2_distances;

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- semi-discrete Wasserstein-2 OT ends ---\n\n');
                fclose(log_file);
            end
        end

        function [objective, W2_distances] = evaluateObjectiveViaW2SAA(obj, atoms, probs, samp_num, rand_stream)
            % Compute the Wasserstein-2 distances from the given discrete measure to the marginals via sample average approximation in 
            % order to evaluate its objective with respect to the Wasserstein barycenter problem
            % Inputs:
            %   atoms: two-column matrix containing the atoms
            %   probs: vector containing the probabilities of the atoms
            %   samp_num: number of samples to generate from each of the marginals for approximation (default is size(atoms, 1))
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream) 
            % Outputs:
            %   objective: the evaluated objective value
            %   W2_distances: vector containing the Wasserstein-2 distances from the given discrete measure to each of the marginals

            if ~exist('samp_num', 'var') || isempty(samp_num)
                samp_num = size(atoms, 1);
            end

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- SAA Wasserstein-2 OT starts ---\n');
            end
            
            marg_num = length(obj.MarginalWeights);
            W2_distances = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};

                marg_samples = marg.randSample(samp_num, rand_stream);

                dist_mat = pdist2(marg_samples, atoms, 'squaredeuclidean');

                cp_options = struct('display', false);
                OT = OTDiscrete(ones(samp_num, 1) / samp_num, probs, dist_mat, cp_options);
                [hcoup_indices, ~] = OT.generateHeuristicCoupling();
                OT.run(hcoup_indices, 1e-5);
                W2_distances(marg_id) = -OT.Runtime.LSIP_LB;

                if obj.Options.display
                    fprintf('%s: marginal %d done\n', class(obj), marg_id);
                end

                % logging
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, '%s: marginal %d done\n', class(obj), marg_id);
                end
            end
            
            objective = obj.MarginalWeights' * W2_distances;

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- SAA Wasserstein-2 OT ends ---\n\n');
                fclose(log_file);
            end
        end

        function [W2OT, W2OT_info_cell] = saveW2OptimalTransportInfo(obj)
            % Retrieve computed semi-discrete Wasserstein-2 optimal transport-related information of each marginal
            % Outputs:
            %   W2OT: struct containing Wasserstein-2 optimal transport-related information that can be saved and loaded later 
            %   W2OT_info_cell: cell array where each cell is a cell array containing Wasserstein-2 optimal transport-related
            %   information about a marginal that can be saved and loaded later

            if ~obj.checkIfWasserstein2OTSupported()
                error('Wasserstein-2 optimal transport is not supported by the marginals');
            end

            marg_num = length(obj.MarginalWeights);

            W2OT = obj.Runtime.W2OT;
            W2OT_info_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                W2OT_info_cell{marg_id} = marg.saveW2OptimalTransportInfo();
            end
        end

        function loadW2OptimalTransportInfo(obj, W2OT, W2OT_info_cell)
            % Load semi-discrete Wasserstein-2 optimal transport-related information into each marginal
            % Inputs:
            %   W2OT: struct containing Wasserstein-2 optimal transport-related information that can be saved and loaded later 
            %   W2OT_info_cell: cell array where each cell is a cell array containing Wasserstein-2 optimal transport-related
            %   information about a marginal that can be saved and loaded later

            if ~obj.checkIfWasserstein2OTSupported()
                error('Wasserstein-2 optimal transport is not supported by the marginals');
            end

            marg_num = length(obj.MarginalWeights);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg.loadW2OptimalTransportInfo(W2OT_info_cell{marg_id});
            end

            obj.Runtime.W2OT = W2OT;
            obj.Runtime.W2OTComputed = true;
        end

        function UB = getMTUpperBoundViaW2DiscreteCoupling(obj)
            % Compute an upper bound for the matching for teams problem via Wasserstein-2 optimal couplings based on the discrete
            % quality measure.
            % Output:
            %   UB: the computed upper bound for the matching for teams problem based on the discrete quality measure

            if ~obj.checkIfWasserstein2OTSupported()
                error('Wasserstein-2 optimal transport is not supported by the marginals');
            end

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            % make sure that the semi-discrete Wasserstein-2 OT problems are solved 
            if ~isfield(obj.Runtime, 'W2OTComputed') || ~obj.Runtime.W2OTComputed
                obj.computeW2OptimalCouplings();
            end

            UB = 0;

            marg_num = length(obj.MarginalWeights);

            for marg_id = 1:marg_num
                UB = UB + obj.MarginalWeights(marg_id) * obj.Marginals{marg_id}.W2OT.Cost;
            end
        end

        function [UB, samps] = getMTUpperBoundViaW2Coupling(obj, ...
                samp_num, rand_stream)
            % Compute an upper bound for the matching for teams problem via Wasserstein-2 optimal couplings based on Monte Carlo 
            % integration. 
            % Inputs: 
            %   samp_num: number of samples for Monte Carlo integration
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream)
            % Output:
            %   UB: the computed upper bound for the matching for teams problem based on Monte Carlo integration
            %   samps: struct containg fields Coupling, Barycenter, and Cost where Coupling is a cell array containing coupled samples 
            %   from the marginals, Barycenter is a two-column matrix containing the corresponding samples from the Wasserstein 
            %   barycenter, and Cost is a vector containing the corresponding values of the cost function

            if ~obj.checkIfWasserstein2OTSupported()
                error('Wasserstein-2 optimal transport is not supported by the marginals');
            end

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            marg_num = length(obj.MarginalWeights);

            coupling_samps = obj.randSampleFromW2Couplings(samp_num, rand_stream);

            barycenters = zeros(samp_num, 2);

            for marg_id = 1:marg_num
                barycenters = barycenters + obj.MarginalWeights(marg_id) * coupling_samps{marg_id};
            end

            costs = zeros(samp_num, 1);

            for marg_id = 1:marg_num
                costs = costs + obj.MarginalWeights(marg_id) * sum((coupling_samps{marg_id} - barycenters) .^ 2, 2);
            end

            samps = struct;
            samps.Coupling = coupling_samps;
            samps.Barycenter = barycenters;
            samps.Cost = costs;

            UB = mean(costs);
        end

        function [UB_list, ...
                samps, ...
                samp_histpdf] ...
                = getMTUpperBoundViaW2CouplingWRepetition(obj, ...
                samp_num, ...
                rep_num, ...
                rand_stream, ...
                hist_edge_x, ...
                hist_edge_y)
            % Compute an upper bound for the matching for teams problem via Wasserstein-2 optimal couplings with repetition of Monte 
            % Carlo integration.
            % Inputs:
            %   samp_num: number of samples for Monte Carlo integration
            %   rep_num: number of repetitions
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream)
            %   hist_edge_x: edges of bins on the x-axis for 2D pdf estimation (default is [])
            %   hist_edge_y: edges of bins on the y-axis for 2D pdf estimation (default is [])
            % Output:
            %   UB_list: vector containing the computed upper bounds for the matching for teams problem
            %   samps: one particular set of Monte Carlo samples used in the approximation of the bounds
            %   samp_histpdf: 2D pdf estimation via the histogram of the generated samples from the approximated Wasserstein 
            %   barycenter; only computed if both hist_edge_x and hist_edge_y are set 

            if ~obj.checkIfWasserstein2OTSupported()
                error('Wasserstein-2 optimal transport is not supported by the marginals');
            end

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~exist('hist_edge_x', 'var') || ~exist('hist_edge_y', 'var')
                hist_edge_x = [];
                hist_edge_y = [];
                samp_histpdf = [];
            end

            UB_list = zeros(rep_num, 1);

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- Monte Carlo sampling starts ---\n');
            end

            % initialize the 2D histogram if it needs to be computed
            if ~isempty(hist_edge_x) && ~isempty(hist_edge_y)
                pdf_computed = true;

                samp_histpdf = zeros(length(hist_edge_x) - 1, length(hist_edge_y) - 1);
            else
                pdf_computed = false;
            end

            for rep_id = 1:rep_num
                [UB_list(rep_id), samps] = obj.getMTUpperBoundViaW2Coupling(samp_num, rand_stream);

                % display output
                if obj.Options.display
                    fprintf('%s: Monte Carlo sampling repetition %3d done\n', class(obj), rep_id);
                end

                % write log
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, '%s: Monte Carlo sampling repetition %3d done\n', class(obj), rep_id);
                end

                % update the 2D histogram
                if pdf_computed
                    samp_histpdf = samp_histpdf + histcounts2(samps.Barycenter(:, 1), samps.Barycenter(:, 2), ...
                        hist_edge_x, hist_edge_y, 'Normalization', 'pdf');
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- Monte Carlo sampling ends ---\n\n');
                fclose(log_file);
            end

            % normalize the 2D histogram
            if pdf_computed
                samp_histpdf = samp_histpdf / rep_num;
            end
        end
        
        function EB = getMTErrorBoundBasedOnOT(obj, LSIP_tolerance)
            % Compute the error bound for the objective value of the matching for teams problem based on the Lipschitz constant of the 
            % cost functions and the optimal transport distances between the marginals and their discretizations as well as the 
            % supremum optimal transport distances between discrete measures on the quality space and the candidate
            % Input:
            %   LSIP_tolerance: tolerance value used in the computation of the LSIP problem
            % Output:
            %   EB: the computed error bound

            % make sure that the semi-discrete OT problems are solved
            if ~isfield(obj.Storage, 'OTComputed') || ~obj.Storage.OTComputed
                obj.performReassembly();
            end

            marg_num = length(obj.MarginalWeights);
            EB = LSIP_tolerance;

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                EB = EB + obj.MarginalWeights(marg_id) * (2 * obj.Quality.MaxNorm * marg.OT.Cost);

                if marg_id ~= obj.Options.discmeas_cand_index
                    EB = EB + 2 * obj.MarginalWeights(marg_id) * (obj.Quality.MaxNorm + marg.Supp.MaxNorm) ...
                        * 2 * obj.Quality.SimplicialTestFuncs.MeshSize;
                end
            end
        end

        function EB = getMTTheoreticalErrorBound(obj, LSIP_tolerance)
            % Compute the theoretical error bound for the objective value of the matching for teams problem via the Lipschitz constant 
            % of the cost functions and the mesh sizes of the simplicial covers for the marginals
            % Input:
            %   LSIP_tolerance: tolerance value used in the computation of the LSIP problem
            % Output:
            %   EB: the computed error bound

            marg_num = length(obj.MarginalWeights);

            EB = LSIP_tolerance;

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                EB = EB + 2 * obj.MarginalWeights(marg_id) * obj.Quality.MaxNorm * 2 * marg.SimplicialTestFuncs.MeshSize;

                if marg_id ~= obj.Options.discmeas_cand_index
                    EB = EB + 2 * obj.MarginalWeights(marg_id) * (obj.Quality.MaxNorm + marg.Supp.MaxNorm) ...
                        * 2 * obj.Quality.SimplicialTestFuncs.MeshSize;
                end
            end
        end
    end

    methods(Access = protected)

        function inside = doCheckIfInsideQualitySpace(obj, ...
                pts)
            % Check if the points are inside the quality space, which is the union of polytopes
            % Input:
            %   pts: two-column matrix containing the input points
            % Output:
            %   inside: boolean vector indicating whether each of the points is inside the quality space

            pt_num = size(pts, 1);
            poly_num = length(obj.Quality.Polytopes);
            inside = any(reshape(all(reshape(obj.Quality.Hyperplane.w * pts' - obj.Quality.Hyperplane.b <= ...
                MatchTeam2DWassersteinBarycenter.INSIDE_TOLERANCE, obj.Quality.Hyperplane.num, poly_num * pt_num), 1), ...
                poly_num, pt_num), 1)'; 
        end

        function initializeBeforeRun(obj)
            % Initialize the algorithm by computing some static quantities

            if ~obj.Storage.SimplicialTestFuncsInitialized
                obj.initializeSimplicialTestFuncs();
            end
        end

        function prepareRuntime(obj)
            % Prepare the runtime environment by initializing some variables
            obj.Runtime = struct;

            marg_num = length(obj.MarginalWeights);

            % the warm-start basis for the constraints are stored in a cell array
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

        function updateRuntimeAfterLP(obj, ...
                result)
            % Update the runtime environment after solving each LP
            % Input:
            %   result: struct produced by gurobi

            % update the warm-start basis used by gurobi
            if isfield(result, 'vbasis') && ~isempty(result.vbasis) && isfield(result, 'cbasis') && ~isempty(result.cbasis)
                quality_vert_num = size(obj.Quality.SimplicialTestFuncs.Vertices, 1) + obj.Quality.SimplicialTestFuncs.Quadratic;

                obj.Runtime.vbasis = result.vbasis;
                obj.Runtime.cbasis = result.cbasis;

                % update the equality and inequality parts of the basis separately
                obj.Runtime.cbasis_eq = obj.Runtime.cbasis(1:quality_vert_num - 1);
                cbasis_ineq = obj.Runtime.cbasis(quality_vert_num:end);

                marg_num = length(obj.MarginalWeights);
                cut_num_list = obj.Runtime.CutNumList;
                cut_counter = 0;

                for marg_id = 1:marg_num
                    obj.Runtime.cbasis_ineq_cell{marg_id} = cbasis_ineq(cut_counter + (1:cut_num_list(marg_id)));
                    cut_counter = cut_counter + cut_num_list(marg_id);
                end
            end
        end

        function [min_pts, ...
                min_vals] ...
                = computeQuadMin(obj, ...
                weights, ...
                batch_size)
            % Solve the minimization of z' * z - 2 * weights' * z where z is in the quality space. Computation is done in batches if
            % necessary to avoid excessive memory usage.
            % Inputs:
            %   weights: two-column matrix where each row contains a weight
            %   batch_size: the maximum number of inputs to be computed in a batch (default is 1e4)
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
                    batch_weight_indices = ((batch_id - 1) * batch_size + 1):min(batch_id * batch_size, remain_num);
                    [min_pts(remain_indices(batch_weight_indices), :), ...
                        min_vals(remain_indices(batch_weight_indices))] ...
                        = obj.doComputeQuadMinOutside( ...
                        weight_remain(batch_weight_indices, :));
                end
            end
        end

        function [min_pts, ...
                min_vals] ...
                = doComputeQuadMinOutside(obj, ...
                weights)
            % Solve the minimization of z' * z - 2 * weights' * z where z is in the quality space and weight is outside the quality
            % space (thus some projection is necessary)
            % Input:
            %   weights: two-column matrix where each row contains a weight
            % Outputs:
            %   min_pts: the computed minimizers
            %   min_vals: the computed minimal values

            weight_num = size(weights, 1);
            min_pts = zeros(weight_num, 2);
            
            % shifting the right-hand side values of the constraints characterizing the polytopes and transform the problem into a
            % Euclidean norm minimization problem
            ineq_rhs_mat = obj.Storage.QuadMin.IneqRHS' - weights * obj.Storage.QuadMin.IneqWeights';
            edgecand_num = size(ineq_rhs_mat, 2);

            % express each candidate point as the affine combination of two end points of an edge of a polytope (if weight is in [0, 1]
            % then the candidate is valid, else it is invalid)
            convcomb_mat = weights * obj.Storage.QuadMin.ConvCombCoefs' + obj.Storage.QuadMin.ConvCombIntercepts';

            % the x- and y-coordinates of the candidate points stored in two matrices
            edgecand_x_mat_shifted = ineq_rhs_mat .* obj.Storage.QuadMin.UnscaledEdgeCandidates(:, 1)';
            edgecand_y_mat_shifted = ineq_rhs_mat .* obj.Storage.QuadMin.UnscaledEdgeCandidates(:, 2)';
            edgecand_x_mat = edgecand_x_mat_shifted + weights(:, 1);
            edgecand_y_mat = edgecand_y_mat_shifted + weights(:, 2);
            
            % the corresponding objective values
            edgecand_val_mat = edgecand_x_mat_shifted .^ 2 + edgecand_y_mat_shifted .^ 2 - sum(weights .^ 2, 2);

            % remove those invalid candidates
            edgecand_val_mat(convcomb_mat < 0 | convcomb_mat > 1) = inf;

            cand_val_mat = [edgecand_val_mat, obj.Storage.QuadMin.VerticesSS' - weights * (2 * obj.Quality.Vertices')];

            % choose a valid candidate with minimum objective value
            [min_vals, min_inds] = min(cand_val_mat, [], 2);
            edge_min_attained_list = min_inds <= edgecand_num;
            edge_min_attained_indices = find(edge_min_attained_list);
            min_attained_linind = sub2ind(size(edgecand_val_mat), edge_min_attained_indices, min_inds(edge_min_attained_indices));
            min_pts(edge_min_attained_indices, 1) = edgecand_x_mat(min_attained_linind);
            min_pts(edge_min_attained_indices, 2) = edgecand_y_mat(min_attained_linind);
            min_pts(~edge_min_attained_list, :) = obj.Quality.Vertices(min_inds(~edge_min_attained_list) - edgecand_num, :);
        end

        function weighted_sum = computeWeightedSumOfVertices(obj, ...
                indices, ...
                batch_size)
            % Compute weighted sum of indices in the triangulation of the marginals, where the weights are specified in
            % obj.MarginalWeights. Computation is done in batches if necessary. 
            % Inputs:
            %   indices: matrix of indices where each row corresponds to an input and each column corresponds to a marginal; each index
            %   corresponds to the index of a vertex in the triangulation of a marginal
            %   batch_size: the maximum number of inputs to be handled in a vectorized procedure (default is 1e4)
            % Output:
            %   weighted_sum: two-column matrix where each row represents a weighted sum

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            input_num = size(indices, 1);
            marg_num = length(obj.MarginalWeights);

            if input_num <= batch_size
                % if the number of input is less than the batch size, do the computation directly

                % first add the offsets to the indices so that they become row indices in the matrix containing all vertices
                vertex_indices = indices' + obj.Storage.MargVertNumOffsets;

                % sparse combination matrix to add up the coordinates of the vertices
                comb_mat = sparse(repelem((1:input_num)', marg_num, 1), ...
                    (1:input_num * marg_num)', ...
                    1, ...
                    input_num, ...
                    input_num * marg_num);

                weighted_sum = full(comb_mat * obj.Storage.WeightedMargVertices(vertex_indices, :));
            else
                % if more than one batch is needed, pre-compute the combination matrix
                batch_num = ceil(input_num / batch_size);
                weighted_sum_cell = cell(batch_num, 1);
                comb_mat = sparse(repelem((1:batch_size)', marg_num, 1), ...
                    (1:batch_size * marg_num)', ...
                    1, ...
                    batch_size, ...
                    batch_size * marg_num);

                for batch_id = 1:batch_num
                    if batch_id < batch_num
                        vertex_indices = indices(((batch_id - 1) * batch_size + 1):batch_id * batch_size, :)' ...
                            + obj.Storage.MargVertNumOffsets;
                        weighted_sum_cell{batch_id} = full(comb_mat * obj.Storage.WeightedMargVertices(vertex_indices, :));
                    else
                        last_batch_size = input_num - (batch_num - 1) * batch_size;
                        vertex_indices = indices(((batch_num - 1) * batch_size + 1):end, :)' + obj.Storage.MargVertNumOffsets;
                        weighted_sum_cell{batch_id} = full(comb_mat(1:last_batch_size, 1:(last_batch_size * marg_num)) ...
                            * obj.Storage.WeightedMargVertices(vertex_indices, :));
                    end
                end

                weighted_sum = vertcat(weighted_sum_cell{:});
            end
        end

        function [vals, ...
                inside] = doEvaluateQualitySimplicialTestFuncs(obj, ...
                pts)
            % Function that actually evaluates the test functions on the quality space at the given locations
            % Input:
            %   pts: two-column matrix containing the input points
            % Output:
            %   vals: sparse matrix containing the computed function values where each row corresponds to an input and each column
            %   corresponds to a test function
            %   inside: sparse boolean matrix indicating whether an input point is inside each of the triangles; each row corresponds
            %   to an input and each column corresponds to a triangle

            input_num = size(pts, 1);
            vert_num = size(obj.Quality.SimplicialTestFuncs.Vertices, 1);
            tri_num = size(obj.Quality.SimplicialTestFuncs.Triangles, 1);

            % apply the inverse transformation to recover the coefficients in the affine combination (with respect to the three 
            % vertices of each triangle)
            coef_mat = obj.Quality.SimplicialTestFuncs.InvTransMat * [pts'; ones(1, input_num)];

            % allow a small tolerance to avoid numerical precision issues
            inside = reshape(all(reshape(sparse(coef_mat(:) > -1e-10), 3, tri_num * input_num), 1), tri_num, input_num)';

            % normalize the matrix such that each row sums up to 1; this is to deal with the case where a point belongs to more than 
            % one triangles; in such a case, we compute the test function with respect to each triangle that the point belongs to and 
            % then take the average
            inside_norm = inside ./ sum(inside, 2);

            % compute the weights in the convex combination
            conv_coef_mat = sparse(coef_mat .* repelem(inside_norm', 3, 1));

            % compute the input indices and the vertex indices corresponding to each non-zero entry in the convex combination 
            % coefficient matrix
            input_id_mat = repmat((1:input_num), tri_num * 3, 1);
            triangles = obj.Quality.SimplicialTestFuncs.Triangles';
            vert_id_mat = repmat(triangles(:), 1, input_num);
            [conv_coef_r, conv_coef_c, conv_coef_vals] = find(conv_coef_mat);
            conv_coef_id = sub2ind([tri_num * 3, input_num], conv_coef_r, conv_coef_c);

            % sum up the coefficients and return a sparse matrix
            vals = accumarray([input_id_mat(conv_coef_id), vert_id_mat(conv_coef_id)], conv_coef_vals, [input_num, vert_num], ...
                [], [], true);
        end

        function model = generateInitialMinModel(obj)
            % Generate the initial linear programming model for gurobi
            % Output:
            %   model: struct containing the linear programming model in gurobi

            model = struct;
            marg_num = length(obj.MarginalWeights);
            marg_vert_num_list = obj.Storage.MargVertNumList;
            quality_vert_num = size(obj.Quality.SimplicialTestFuncs.Vertices, 1) + obj.Quality.SimplicialTestFuncs.Quadratic;
            decivar_num = obj.Storage.DeciVarLength;

            objective_cell = cell(marg_num, 1);
            A_eq_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};

                % the integrals are placed into the objective vector; note that the first test function for the marginal is removed for
                % identification
                marg_integral = marg.SimplicialTestFuncs.Integrals;

                if ~obj.Options.sparse_parametrization
                    objective_cell{marg_id} = [-1; -marg_integral(2:end); zeros(quality_vert_num - 1, 1)];
                else
                    objective_cell{marg_id} = [-marg_integral; zeros(quality_vert_num - 1, 1)];
                end

                % matrix for the equality constraints requiring that the transfer functions must sum up to 0; note that the first test 
                % function on the quality space is removed for identification
                A_eq_cell{marg_id} = [sparse(quality_vert_num - 1, length(marg_integral)), speye(quality_vert_num - 1)];
            end

            % since the cutting plane algorithm assumes that the problem is a minimization problem, we need to transform our 
            % maximization problem into a minimization problem
            model.modelsense = 'min';

            % the constant in the objective is the sum of the quadratic constants that do not affect the matching
            model.objcon = -obj.QuadraticConstant;

            % the coefficients corresponding to the first test function of each marginal is not included in the decision variables for
            % identification purposes
            model.obj = vertcat(objective_cell{:});

            model.lb = -inf(decivar_num, 1);
            model.ub = inf(decivar_num, 1);

            % store the equality constraints as fields of the model
            model.A_eq = horzcat(A_eq_cell{:});
            model.rhs_eq = zeros(quality_vert_num - 1, 1);

            % store the inequality constraints into cell arrays as fields of the model
            model.A_ineq_cell = cell(marg_num, 1);
            model.rhs_ineq_cell = cell(marg_num, 1);

            for marg_id = 1:marg_num
                model.A_ineq_cell{marg_id} = sparse(0, marg_vert_num_list(marg_id) - 1);
                model.rhs_ineq_cell{marg_id} = zeros(0, 1);
            end

            % the constraints in the model need to be assembled from these fields
            model.A = [model.A_eq; blkdiag(model.A_ineq_cell{:})];
            model.rhs = [model.rhs_eq; vertcat(model.rhs_ineq_cell{:})];
            model.sense = [repmat('=', length(model.rhs_eq), 1); repmat('>', length(model.rhs) - length(model.rhs_eq), 1)];
        end

        function [min_lb, ...
                optimizers] ...
                = callGlobalMinOracle(obj, ...
                vec)
            % Given a decision vector, call the global minimization oracle to approximately determine the "most violated" constraints
            % and return a lower bound for the minimal value
            % Input:
            %   vec: vector corresponding to the current LP solution
            % Outputs:
            %   min_lb: lower bound for the global minimization problem
            %   optimizers: cell array containing structs with fields
            %   vertex_indices and points that correspond to the approximate optimizers of the global minimization problems

            marg_num = length(obj.MarginalWeights);
            intecept_indices = obj.Storage.DeciVarInterceptIndices;
            marg_tfuncs_indices = obj.Storage.DeciVarMargTestFuncIndices;
            q_tfuncs_indices = obj.Storage.DeciVarQualityTestFuncIndices;

            min_vals_list = zeros(marg_num, 1);
            optimizers = cell(marg_num, 1);

            for marg_id = 1:marg_num
                % disassemble the vector resulted from solving the LP problem
                if ~obj.Options.sparse_parametrization
                    % if the constant intercept is present, add back the omitted first test function on the support of the marginal
                    objective_const = vec(intecept_indices(marg_id));
                    marg_vert_vals = [0; vec(marg_tfuncs_indices{marg_id})];
                else
                    % otherwise, set the constant intercept to 0 and take the full set of test functions
                    objective_const = 0;
                    marg_vert_vals = vec(marg_tfuncs_indices{marg_id});
                end

                % the first test function is added back to the list of test functions on the quality space
                quality_vert_vals = [0; vec(q_tfuncs_indices{marg_id})];

                % compute a pool of approximate solutions
                if ~obj.Quality.SimplicialTestFuncs.Quadratic
                    quality_quad_coef = 0;
                else
                    quality_quad_coef = vec(obj.Storage.DeciVarQualityQuadIndices(marg_id));
                    quality_vert_vals = quality_vert_vals - quality_quad_coef * obj.Quality.SimplicialTestFuncs.QuadraticValues;
                end

                [pool_inds, ...
                    pool_pts, ...
                    pool_testfuncs, ...
                    pool_vals, ...
                    pool_costs, ...
                    min_val] ...
                    = obj.computeGlobalMin( ...
                    marg_id, ...
                    quality_quad_coef, ...
                    quality_vert_vals, ...
                    marg_vert_vals, ...
                    objective_const); %#ok<ASGLU>

                if ~isempty(min_val)
                    min_vals_list(marg_id) = min_val;
                else
                    min_vals_list(marg_id) = 0;
                end

                optimizers{marg_id} = struct( ...
                    'vertex_indices', pool_inds, ...
                    'points', pool_pts, ...
                    'testfuncs_vals', pool_testfuncs, ...
                    'costfunc_vals', pool_costs, ...
                    'min_val', min_vals_list(marg_id));
            end

            min_lb = sum(min_vals_list);
        end

        function [pool_inds, ...
                pool_pts, ...
                pool_testfuncs, ...
                pool_vals, ...
                pool_costs, ...
                min_val] ...
                = computeGlobalMin(obj, ...
                marg_id, ...
                quality_quad_coef, ...
                quality_vert_vals, ...
                marg_vert_vals, ...
                objective_const)
            % Solve the global minimization problem 
            % Input:
            %   marg_id: the index of the marginal with respect to which the global minimization problem needs to be solved
            %   quailty_quad_coef: coefficient of the quadratic test function on the quality space
            %   quality_vert_vals: values of the simplicial test functions on the quality spaces at the vertices
            %   marg_vert_vals: values of the simplicial test functions for the marginals at the vertices
            %   objective_const: constant part of the objective function
            % Outputs:
            %   pool_inds: matrix containing approximate optimizers of the global minimization problem in the form of vertex indices in
            %   the triangulations of the input measures; the algorithm will retain up to obj.GlobalOptions.pool_size approximate
            %   optimizers
            %   pool_pts: two-column matrix containing the approximate optimizers in the quality space
            %   pool_vals: the corresponding objective values of the approximate optimizers
            %   pool_costs: the corresponding values of the cost function evaluated at the approximate optimizers
            %   min_val: the computed minimal value

            % open the log file
            if ~isempty(obj.GlobalOptions.log_file)
                log_file = fopen(obj.GlobalOptions.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- Function call starts ---\n');
            end


            boundary_buffer = obj.GlobalOptions.boundary_buffer;

            marg_weight = obj.MarginalWeights(marg_id);
            marg_vertices = obj.Marginals{marg_id}.SimplicialTestFuncs.Vertices;

            quality_triangles = obj.Quality.SimplicialTestFuncs.Triangles;
            quality_vertices = obj.Quality.SimplicialTestFuncs.Vertices;
            quality_vertices1 = quality_vertices(:, 1);
            quality_vertices2 = quality_vertices(:, 2);
            quality_vert_num = size(obj.Quality.SimplicialTestFuncs.Vertices, 1);

            if obj.Quality.SimplicialTestFuncs.Quadratic
                quality_vert_quad_vals = obj.Quality.SimplicialTestFuncs.QuadraticValues;
            end

            coef1 = marg_weight / (marg_weight - quality_quad_coef);
            coef2 = 1 / 2 / (marg_weight - quality_quad_coef);

            GM = obj.Storage.GlobalMin;
            GM_marg = GM.Marginals{marg_id};

            % outputs from batches are stored in cell arrays
            pool_costs_cell = cell(GM_marg.BatchNum, 1);
            pool_vals_cell = cell(GM_marg.BatchNum, 1);
            pool_pts_cell = cell(GM_marg.BatchNum, 1);
            pool_inds_cell = cell(GM_marg.BatchNum, 1);
            pool_testfuncs_cell = cell(GM_marg.BatchNum, 1);

            for batch_id = 1:GM_marg.BatchNum
                marg_batch_vertices = marg_vertices(GM_marg.BatchLists{batch_id}, :);
                marg_batch_vert_vals = marg_vert_vals(GM_marg.BatchLists{batch_id});
                marg_batch_vert_num = length(GM_marg.BatchLists{batch_id});

                % the vertex case: compute the objective function at each vertex
                case_vert_vert_quad_vals = sum(quality_vertices.^2, 2);
                case_vert_costs = marg_weight * case_vert_vert_quad_vals - 2 * marg_weight * quality_vertices * marg_batch_vertices';
                case_vert_vals = case_vert_costs - quality_quad_coef * case_vert_vert_quad_vals ...
                    - quality_vert_vals - marg_batch_vert_vals' - objective_const;

                % filter out all the combinations with non-negative objectives
                [case_vert_neg_r, case_vert_neg_c] = find(case_vert_vals < obj.GlobalOptions.objective_threshold);
                case_vert_neg_ind = sub2ind([quality_vert_num, marg_batch_vert_num], case_vert_neg_r, case_vert_neg_c);

                % reformat the values into a vector for sorting
                case_vert_costs_vec = case_vert_costs(case_vert_neg_ind);
                case_vert_vals_vec = case_vert_vals(case_vert_neg_ind);
                case_vert_q_inds = case_vert_neg_r;
                case_vert_m_inds = GM_marg.BatchLists{batch_id}(case_vert_neg_c);
                case_vert_num = length(case_vert_vals_vec);

                if marg_weight - quality_quad_coef <= 1e-6
                    % if the total quadratic coefficient is non-positive (or approximately so), it is sufficient to consider only the
                    % vertex case

                    pool_size = min(obj.GlobalOptions.pool_size, length(case_vert_vals_vec));

                    [case_vert_vals_vec, s_ind] = sort(case_vert_vals_vec, 'ascend');
                    pool_sort = s_ind(1:pool_size);
                    pool_vals_cell{batch_id} = case_vert_vals_vec(1:pool_size);
                    pool_inds_cell{batch_id} = case_vert_m_inds(pool_sort);
                    pool_q_inds = case_vert_q_inds(pool_sort);
                    pool_pts_cell{batch_id} = quality_vertices(pool_q_inds, :);
                    pool_costs_cell{batch_id} = case_vert_costs_vec(pool_sort);

                    % depending on whether the quadratic test function is enabled, there will be an extra test function
                    pool_testfuncs_cell{batch_id} = sparse( ...
                        1:pool_size, ...
                        pool_q_inds, ...
                        1, ...
                        pool_size, quality_vert_num + obj.Quality.SimplicialTestFuncs.Quadratic);
                    continue;
                end

                % the face case: compute the potentially minimal values within the interior of each face
                face_num = size(quality_triangles, 1);

                % evaluate the two coefficients and store the results in two matrices
                A1 = GM.CaseFace.Fixed{1} + coef1 * GM.CaseFace.Marginal{1} * marg_batch_vertices' ...
                    + coef2 * full(GM.CaseFace.TestFunc{1} * quality_vert_vals);
                A2 = GM.CaseFace.Fixed{2} + coef1 * GM.CaseFace.Marginal{2} * marg_batch_vertices' ...
                    + coef2 * full(GM.CaseFace.TestFunc{2} * quality_vert_vals);

                % only treat those coefficients above obj.GlobalOptions.boundary_buffer as valid
                [case_face_valid_r, case_face_valid_c] = find((A1 >= boundary_buffer) & (A2 >= boundary_buffer) ...
                    & (A1 + A2 <= 1 - boundary_buffer));
                case_face_valid_indices = sub2ind([face_num, marg_batch_vert_num], case_face_valid_r, case_face_valid_c);
                case_face_coef1 = A1(case_face_valid_indices);
                case_face_coef2 = A2(case_face_valid_indices);
                case_face_coef_mat = [case_face_coef1, case_face_coef2, 1 - (case_face_coef1 + case_face_coef2)];
                case_face_m_inds = GM_marg.BatchLists{batch_id}(case_face_valid_c);
                case_face_tri_id = quality_triangles(case_face_valid_r, :);
                case_face_pts = [sum(case_face_coef_mat .* quality_vertices1(case_face_tri_id')', 2), ...
                    sum(case_face_coef_mat .* quality_vertices2(case_face_tri_id')', 2)];
                case_face_pts_quad_vals = sum(case_face_pts.^2, 2);
                case_face_costs_vec = marg_weight * case_face_pts_quad_vals ...
                    - 2 * marg_weight * sum(case_face_pts .* marg_batch_vertices(case_face_valid_c, :), 2);
                case_face_vals_vec = case_face_costs_vec - quality_quad_coef * case_face_pts_quad_vals ...
                    - marg_batch_vert_vals(case_face_valid_c) ...
                    - sum(case_face_coef_mat .* quality_vert_vals(case_face_tri_id')', 2) - objective_const;

                % filter out the points with non-negative objectives
                case_face_neg_indices = case_face_vals_vec < obj.GlobalOptions.objective_threshold;
                case_face_coef_mat = case_face_coef_mat(case_face_neg_indices, :);
                case_face_m_inds = case_face_m_inds(case_face_neg_indices);
                case_face_tri_id = case_face_tri_id(case_face_neg_indices, :);
                case_face_pts = case_face_pts(case_face_neg_indices, :);
                case_face_costs_vec = case_face_costs_vec(case_face_neg_indices);
                case_face_vals_vec = case_face_vals_vec(case_face_neg_indices);
                case_face_num = length(case_face_vals_vec);

                % the edge case: compute the potentially minimal values within the relative interior of each edge
                edge_list = GM.CaseEdge.Edges;
                edge_num = size(edge_list, 1);

                % evaluate the coefficient and store the results in a matrix
                A1 = GM.CaseEdge.Fixed + coef1 * GM.CaseEdge.Marginal * marg_batch_vertices' ...
                    + coef2 * full(GM.CaseEdge.TestFunc * quality_vert_vals);

                % only treat those coefficients above obj.GlobalOptions.boundary_buffer as valid
                [case_edge_valid_r, case_edge_valid_c] = find((A1 >= boundary_buffer) & (A1 <= 1 - boundary_buffer));
                case_edge_valid_indices = sub2ind([edge_num, marg_batch_vert_num], case_edge_valid_r, case_edge_valid_c);
                case_edge_coef = A1(case_edge_valid_indices);
                case_edge_coef_mat = [case_edge_coef, 1 - case_edge_coef];
                case_edge_m_inds = GM_marg.BatchLists{batch_id}(case_edge_valid_c);
                case_edge_edge_id = edge_list(case_edge_valid_r, :);
                case_edge_pts = [sum(case_edge_coef_mat .* quality_vertices1(case_edge_edge_id')', 2), ...
                    sum(case_edge_coef_mat .* quality_vertices2(case_edge_edge_id')', 2)];
                case_edge_pts_quad_vals = sum(case_edge_pts.^2, 2);
                case_edge_costs_vec = marg_weight * case_edge_pts_quad_vals ...
                    - 2 * marg_weight * sum(case_edge_pts .* marg_batch_vertices(case_edge_valid_c, :), 2);
                case_edge_vals_vec = case_edge_costs_vec - quality_quad_coef * case_edge_pts_quad_vals ...
                    - marg_batch_vert_vals(case_edge_valid_c) ...
                    - sum(case_edge_coef_mat .* quality_vert_vals(case_edge_edge_id')', 2) - objective_const;

                % filter out the points with non-negative objectives
                case_edge_neg_indices = case_edge_vals_vec < obj.GlobalOptions.objective_threshold;
                case_edge_coef_mat = case_edge_coef_mat(case_edge_neg_indices, :);
                case_edge_m_inds = case_edge_m_inds(case_edge_neg_indices);
                case_edge_edge_id = case_edge_edge_id(case_edge_neg_indices, :);
                case_edge_pts = case_edge_pts(case_edge_neg_indices, :);
                case_edge_costs_vec = case_edge_costs_vec(case_edge_neg_indices);
                case_edge_vals_vec = case_edge_vals_vec(case_edge_neg_indices);
                case_edge_num = length(case_edge_vals_vec);

                % combine all the candidate from the three cases and sort the objective values
                pool_size = min(obj.GlobalOptions.pool_size, case_vert_num + case_edge_num + case_face_num);

                [~, s_ind] = sort([case_vert_vals_vec; case_edge_vals_vec; case_face_vals_vec], 'ascend');

                % split the sorting indices back into the three cases
                pool_sort = s_ind(1:pool_size);
                pool_sort_vert = pool_sort(pool_sort <= case_vert_num);
                pool_sort_edge = pool_sort(pool_sort > case_vert_num & pool_sort <= case_vert_num + case_edge_num) - case_vert_num;
                pool_sort_face = pool_sort(pool_sort > case_vert_num + case_edge_num) - (case_vert_num + case_edge_num);


                % prepare the outputs for the vertex case
                pool_vert_costs = case_vert_costs_vec(pool_sort_vert);
                pool_vert_vals = case_vert_vals_vec(pool_sort_vert);
                pool_vert_q_inds = case_vert_q_inds(pool_sort_vert);
                pool_vert_pts = quality_vertices(pool_vert_q_inds, :);
                pool_vert_inds = case_vert_m_inds(pool_sort_vert);

                % depending on whether the quadratic test function is enabled, there will be an extra test function
                pool_vert_testfuncs = sparse( ...
                    1:length(pool_sort_vert), ...
                    pool_vert_q_inds, ...
                    1, ...
                    length(pool_sort_vert), quality_vert_num + obj.Quality.SimplicialTestFuncs.Quadratic);

                % prepare the outputs for the edge case
                pool_edge_costs = case_edge_costs_vec(pool_sort_edge);
                pool_edge_vals = case_edge_vals_vec(pool_sort_edge);
                pool_edge_pts = case_edge_pts(pool_sort_edge, :);
                pool_edge_inds = case_edge_m_inds(pool_sort_edge);
                pool_edge_coefs = case_edge_coef_mat(pool_sort_edge, :);
                pool_edge_edge_id = case_edge_edge_id(pool_sort_edge, :);

                % depending on whether the quadratic test function is enabled, there will be an extra test function
                pool_edge_testfuncs = sparse( ...
                    repmat((1:length(pool_sort_edge))', 2, 1), ...
                    pool_edge_edge_id(:), ...
                    pool_edge_coefs(:), ...
                    length(pool_sort_edge), quality_vert_num);

                if obj.Quality.SimplicialTestFuncs.Quadratic
                    pool_edge_quad_vals = sum(pool_edge_pts.^2, 2) ...
                        - sum(quality_vert_quad_vals(pool_edge_edge_id')' .* pool_edge_coefs, 2);
                    pool_edge_testfuncs = [pool_edge_testfuncs, sparse(pool_edge_quad_vals)]; %#ok<AGROW>
                end


                % prepare the outputs for the face case
                pool_face_costs = case_face_costs_vec(pool_sort_face);
                pool_face_vals = case_face_vals_vec(pool_sort_face);
                pool_face_pts = case_face_pts(pool_sort_face, :);
                pool_face_inds = case_face_m_inds(pool_sort_face);
                pool_face_coefs = case_face_coef_mat(pool_sort_face, :);
                pool_face_tri_id = case_face_tri_id(pool_sort_face, :);

                % depending on whether the quadratic test function is enabled, there will be an extra test function
                pool_face_testfuncs = sparse( ...
                    repmat((1:length(pool_sort_face))', 3, 1), ...
                    pool_face_tri_id(:), ...
                    pool_face_coefs(:), ...
                    length(pool_sort_face), quality_vert_num);

                if obj.Quality.SimplicialTestFuncs.Quadratic
                    pool_face_quad_vals = sum(pool_face_pts.^2, 2) ...
                        - sum(quality_vert_quad_vals(pool_face_tri_id')' .* pool_face_coefs, 2);
                    pool_face_testfuncs = [pool_face_testfuncs, sparse(pool_face_quad_vals)]; %#ok<AGROW>
                end


                % combine the outputs in all three cases
                pool_costs_cell{batch_id} = [pool_vert_costs; pool_edge_costs; pool_face_costs];
                pool_vals_cell{batch_id} = [pool_vert_vals; pool_edge_vals; pool_face_vals];
                pool_pts_cell{batch_id} = [pool_vert_pts; pool_edge_pts; pool_face_pts];
                pool_inds_cell{batch_id} = [pool_vert_inds; pool_edge_inds; pool_face_inds];
                pool_testfuncs_cell{batch_id} = [pool_vert_testfuncs; pool_edge_testfuncs; pool_face_testfuncs];
            end

            % aggregate the batches
            min_inds = vertcat(pool_inds_cell{:});
            min_pts = vertcat(pool_pts_cell{:});
            min_testfuncs = vertcat(pool_testfuncs_cell{:});
            min_vals = vertcat(pool_vals_cell{:});
            min_costs = vertcat(pool_costs_cell{:});
            min_val = min(min_vals);

            % sort the optimizers from the batches
            [min_vals, sorted_order] = sort(min_vals, 1, 'ascend');
            min_costs = min_costs(sorted_order);
            min_inds = min_inds(sorted_order);
            min_pts = min_pts(sorted_order, :);
            min_testfuncs = min_testfuncs(sorted_order, :);

            pool_size = min(obj.GlobalOptions.pool_size, length(min_vals));
            pool_inds = min_inds(1:pool_size);
            pool_pts = min_pts(1:pool_size, :);
            pool_testfuncs = min_testfuncs(1:pool_size, :);
            pool_vals = min_vals(1:pool_size);
            pool_costs = min_costs(1:pool_size);

            if obj.GlobalOptions.display
                fprintf('%s: Found %3d candidates, minimum value = %10.4f\n', class(obj), pool_size, min_val);
            end

            if ~isempty(obj.GlobalOptions.log_file)
                fprintf(log_file, '%s: Found %3d candidates, minimum value = %10.4f\n', class(obj), pool_size, min_val);
            end

            % close the log file
            if ~isempty(obj.GlobalOptions.log_file)
                fprintf(log_file, '--- Function call ends ---\n\n');
                fclose(log_file);
            end
        end

        function updateLSIPUB(obj, ...
                min_lb, ...
                optimizers) 
            % Update the LSIP upper bound after each call to the global minimization oracle
            % Inputs:
            %   min_lb: the lower bound for the global minimization problem
            %   optimizers: a set of approximate optimizers of the global minimization problem

            obj.Runtime.LSIP_UB = min(obj.Runtime.LSIP_UB, obj.Runtime.LSIP_LB - min_lb);

            marg_num = length(obj.MarginalWeights);
            obj.Runtime.GlobalMin = struct;
            obj.Runtime.GlobalMin.MinVals = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                obj.Runtime.GlobalMin.MinVals(marg_id) = optimizers{marg_id}.min_val;
            end
        end
        
        function addConstraints(obj, ...
                optimizers)
            % Given a collection of approximate optimizers from the global minimization oracle, generate and add the corresponding
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
                    % numerically unstable; should be avoided
                    if obj.Quality.SimplicialTestFuncs.Quadratic
                        simp_testfunc_vals = obj.evaluateQualitySimplicialTestFuncs(optimizers_i.points);
                        quality_testfunc_vals = [simp_testfunc_vals, ...
                            sum(optimizers_i.points.^2, 2) - simp_testfunc_vals * obj.Quality.SimplicialTestFuncs.QuadraticValues];
                    else
                        quality_testfunc_vals = obj.evaluateQualitySimplicialTestFuncs(optimizers_i.points);
                    end
                end

                % first generate a matrix containing all test functions (each row corresponds to an approximate optimizer, each column 
                % corresponds to a test function)
                A_marg_full = sparse( ...
                    (1:constr_num)', ...
                    optimizers_i.vertex_indices, ...
                    1, ...
                    constr_num, obj.Storage.MargVertNumList(marg_id));

                if ~obj.Options.sparse_parametrization
                    % if the constant intercept is present, remove the first test function for identification, then prepend a column of
                    % 1 and add the values of the simplicial test functions on the quality space on the right
                    A_new = [sparse(ones(constr_num, 1)), A_marg_full(:, 2:end), sparse(quality_testfunc_vals(:, 2:end))];
                else
                    % otherwise, omit the constant intercept, use the full set of test functions on the support of the marginal and
                    % remove the first test function on the quality space for identification
                    A_new = [A_marg_full, sparse(quality_testfunc_vals(:, 2:end))];
                end

                rhs_new = optimizers_i.costfunc_vals;

                % add the newly generated constraints to the end
                obj.Runtime.CurrentLPModel.A_ineq_cell{marg_id} = [obj.Runtime.CurrentLPModel.A_ineq_cell{marg_id}; -A_new];
                obj.Runtime.CurrentLPModel.rhs_ineq_cell{marg_id} = [obj.Runtime.CurrentLPModel.rhs_ineq_cell{marg_id}; -rhs_new];

                % add the added indices to the runtime environment
                obj.Runtime.CutIndices{marg_id} = [obj.Runtime.CutIndices{marg_id}; optimizers_i.vertex_indices];
                obj.Runtime.CutPoints{marg_id} = [obj.Runtime.CutPoints{marg_id}; optimizers_i.points];
                obj.Runtime.CutNumList(marg_id) = obj.Runtime.CutNumList(marg_id) + constr_num;

                % update the logical vector indicating the initial constraints
                if isfield(obj.Runtime, 'InitialConstraintsList') && ~isempty(obj.Runtime.InitialConstraintsList)
                    obj.Runtime.InitialConstraintsList{marg_id} = [obj.Runtime.InitialConstraintsList{marg_id}; false(constr_num, 1)];
                end
            end

            % the constraints in the model need to be assembled again
            obj.Runtime.CurrentLPModel.A = [obj.Runtime.CurrentLPModel.A_eq; blkdiag(obj.Runtime.CurrentLPModel.A_ineq_cell{:})];
            obj.Runtime.CurrentLPModel.rhs = [obj.Runtime.CurrentLPModel.rhs_eq; vertcat(obj.Runtime.CurrentLPModel.rhs_ineq_cell{:})];
            obj.Runtime.CurrentLPModel.sense = [repmat('=', length(obj.Runtime.CurrentLPModel.rhs_eq), 1); ...
                repmat('>', length(obj.Runtime.CurrentLPModel.rhs) - length(obj.Runtime.CurrentLPModel.rhs_eq), 1)];

            if ~isempty(obj.Runtime.vbasis) && ~isempty(obj.Runtime.cbasis)
                obj.Runtime.CurrentLPModel.vbasis = obj.Runtime.vbasis;
                
                % set the warm-start basis of the new constraints to 0
                for marg_id = 1:marg_num
                    optimizers_i = optimizers{marg_id};
                    constr_num = size(optimizers_i.vertex_indices, 1);
                    obj.Runtime.cbasis_ineq_cell{marg_id} = [obj.Runtime.cbasis_ineq_cell{marg_id}; zeros(constr_num, 1)];
                end

                obj.Runtime.CurrentLPModel.cbasis = [obj.Runtime.cbasis_eq; vertcat(obj.Runtime.cbasis_ineq_cell{:})];
            else
                if isfield(obj.Runtime.CurrentLPModel, 'vbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'vbasis');
                end

                if isfield(obj.Runtime.CurrentLPModel, 'cbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'cbasis');
                end
            end

            % if obj.Runtime.NumOfInitialConstraints is not set, it means that this is the first call to obj.addConstraints which
            % generates the initial constraints; this number is stored in the runtime environment
            if ~isfield(obj.Runtime, 'NumOfInitialConstraints') || isempty(obj.Runtime.NumOfInitialConstraints)
                obj.Runtime.NumOfInitialConstraints = zeros(marg_num, 1);
                obj.Runtime.InitialConstraintsList = cell(marg_num, 1);

                for marg_id = 1:marg_num
                    optimizers_i = optimizers{marg_id};
                    constr_num = size(optimizers_i.vertex_indices, 1);
                    obj.Runtime.NumOfInitialConstraints(marg_id) = constr_num;
                    obj.Runtime.InitialConstraintsList{marg_id} = true(constr_num, 1);
                end
            end
        end

        function reduceConstraints(obj, result)
            % Remove some of the constraints to speed up the LP solver
            % Input:
            %   result: the output from the gurobi LP solver

            if isinf(obj.Options.reduce.thres) || obj.Runtime.iter <= 0 || obj.Runtime.iter > obj.Options.reduce.max_iter ...
                    || mod(obj.Runtime.iter, obj.Options.reduce.freq) ~= 0
                return;
            end

            marg_num = length(obj.MarginalWeights);
            quality_vert_num = size(obj.Quality.SimplicialTestFuncs.Vertices, 1) + obj.Quality.SimplicialTestFuncs.Quadratic;
            cut_num_list = obj.Runtime.CutNumList;

            cut_slack = result.slack(quality_vert_num:end);
            flagged_indices = cut_slack < -obj.Options.reduce.min_slack;

            if obj.Options.reduce.preserve_init_constr && isfield(obj.Runtime, 'InitialConstraintsList') ...
                    && ~isempty(obj.Runtime.InitialConstraintsList)
                flagged_indices = flagged_indices & ~vertcat(obj.Runtime.InitialConstraintsList{:});
            end

            cut_slack_thres_quantile = quantile(-cut_slack(flagged_indices), obj.Options.reduce.thres_quantile);

            cut_counter = 0;

            for marg_id = 1:marg_num
                % the list of constraints to be kept (here since the directions of the inequalities are all >=, the slackness is 
                % non-positive; the thresholds specifies the maximum absolute value of slackness)
                slack = cut_slack(cut_counter + (1:cut_num_list(marg_id)));
                keep_list = slack >= -obj.Options.reduce.thres & slack >= -cut_slack_thres_quantile;

                % keep the initial constraints depending on obj.Options.reduce.preserve_init_constr
                if obj.Options.reduce.preserve_init_constr
                    keep_list(1:obj.Runtime.NumOfInitialConstraints(marg_id)) = true;
                end

                % update all variables
                obj.Runtime.CutIndices{marg_id} = obj.Runtime.CutIndices{marg_id}(keep_list);
                obj.Runtime.CutPoints{marg_id} = obj.Runtime.CutPoints{marg_id}(keep_list, :);
                obj.Runtime.CutNumList(marg_id) = sum(keep_list);
                obj.Runtime.CurrentLPModel.A_ineq_cell{marg_id} = obj.Runtime.CurrentLPModel.A_ineq_cell{marg_id}(keep_list, :);
                obj.Runtime.CurrentLPModel.rhs_ineq_cell{marg_id} = obj.Runtime.CurrentLPModel.rhs_ineq_cell{marg_id}(keep_list);

                if ~isempty(obj.Runtime.cbasis)
                    obj.Runtime.cbasis_ineq_cell{marg_id} = obj.Runtime.cbasis_ineq_cell{marg_id}(keep_list);
                end

                if isfield(obj.Runtime, 'InitialConstraintsList') && ~isempty(obj.Runtime.InitialConstraintsList)
                    obj.Runtime.InitialConstraintsList{marg_id} = obj.Runtime.InitialConstraintsList{marg_id}(keep_list);
                end

                cut_counter = cut_counter + cut_num_list(marg_id);
            end

            if ~isempty(obj.Runtime.cbasis)
                obj.Runtime.CurrentLPModel.cbasis = [obj.Runtime.cbasis_eq; vertcat(obj.Runtime.cbasis_ineq_cell{:})];
            end

            % the constraints in the model need to be assembled again
            obj.Runtime.CurrentLPModel.A = [obj.Runtime.CurrentLPModel.A_eq; blkdiag(obj.Runtime.CurrentLPModel.A_ineq_cell{:})];
            obj.Runtime.CurrentLPModel.rhs = [obj.Runtime.CurrentLPModel.rhs_eq; vertcat(obj.Runtime.CurrentLPModel.rhs_ineq_cell{:})];
            obj.Runtime.CurrentLPModel.sense = [repmat('=', length(obj.Runtime.CurrentLPModel.rhs_eq), 1); ...
                repmat('>', length(obj.Runtime.CurrentLPModel.rhs) - length(obj.Runtime.CurrentLPModel.rhs_eq), 1)];
        end

        function primal_sol = buildPrimalSolution(obj, result, ~)
            % Given the output from gurobi and a lower bound for the optimal value of the global minimization oracle, build the
            % corresponding primal solution
            % Inputs:
            %   result: output of the gurobi LP solver
            %   violation: a lower bound for the global minimization problem
            % Output:
            %   primal_sol: the constructed potential functions on the support of the input measures as well as the transfer functions 
            %   on the quality space

            marg_num = length(obj.MarginalWeights);
            vec = result.x;
            primal_sol = cell(marg_num, 1);

            for marg_id = 1:marg_num
                primal_sol{marg_id} = struct;

                if ~obj.Options.sparse_parametrization
                    primal_sol{marg_id}.Constant = vec(obj.Storage.DeciVarInterceptIndices(marg_id)) ...
                        - obj.Runtime.GlobalMin.MinVals(marg_id);
    
                    % add back the first test function for the marginal
                    primal_sol{marg_id}.Coefficients = [0; vec(obj.Storage.DeciVarMargTestFuncIndices{marg_id})];
                else
                    % if the constant intercept does not appear in the parametrization, make the coefficient corresponding to the first
                    % test function equal to 0 for consistency with the other case
                    primal_sol{marg_id}.Coefficients = vec(obj.Storage.DeciVarMargTestFuncIndices{marg_id});
                    primal_sol{marg_id}.Constant = primal_sol{marg_id}.Coefficients(1);
                    primal_sol{marg_id}.Coefficients = primal_sol{marg_id}.Coefficients - primal_sol{marg_id}.Constant;
                end

                % add back the first test function on the quality space
                primal_sol{marg_id}.TransFuncCoefficients = [0; vec(obj.Storage.DeciVarQualityTestFuncIndices{marg_id})];

                if obj.Quality.SimplicialTestFuncs.Quadratic
                    primal_sol{marg_id}.TransFuncQuadCoefficient = vec(obj.Storage.DeciVarQualityQuadIndices(marg_id));
                end
            end
        end

        function dual_sol = buildDualSolution(obj, result)
            % Given the output from gurobi, build the corresponding dual solution
            % Input:
            %   result: output of the gurobi LP solver
            % Output:
            %   dual_sol: the constructed discrete probability measure for the relaxed problem coupled with the corresponding discrete 
            %   probability measure on the quality space

            marg_num = length(obj.MarginalWeights);
            quality_vert_num = size(obj.Quality.SimplicialTestFuncs.Vertices, 1) + obj.Quality.SimplicialTestFuncs.Quadratic;
            cut_num_list = obj.Runtime.CutNumList;
            cut_dual = result.pi(quality_vert_num:end);

            dual_sol = cell(marg_num, 1);
            cut_counter = 0;

            for marg_id = 1:marg_num
                marg_cut_dual = cut_dual(cut_counter + (1:cut_num_list(marg_id)));
                pos_list = marg_cut_dual >= obj.Options.dual_prob_thres;
                dual_sol{marg_id} = struct;
                dual_sol{marg_id}.Probabilities = marg_cut_dual(pos_list);
                dual_sol{marg_id}.VertexIndices = obj.Runtime.CutIndices{marg_id}(pos_list);
                dual_sol{marg_id}.QualityAtoms = obj.Runtime.CutPoints{marg_id}(pos_list, :);

                % normalize the probabilities to resolve small numerical inaccuracies
                dual_sol{marg_id}.Probabilities = dual_sol{marg_id}.Probabilities / sum(dual_sol{marg_id}.Probabilities);

                cut_counter = cut_counter + cut_num_list(marg_id);
            end
        end

        function prepareDiscreteOT(obj)
            % Compute discrete to discrete optimal transport between the candidate discrete measure on the quality space and the rest
            % of the discrete measures on the quality space

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            cand_id = obj.Options.discmeas_cand_index;
            cand = obj.Runtime.DualSolution{cand_id};
            [~, cand_uind, cand_umap] = unique(round(cand.QualityAtoms, 6), 'rows', 'stable');
            cand_atoms = cand.QualityAtoms(cand_uind, :);
            cand_atom_num = size(cand_atoms, 1);
            cand_probs = accumarray(cand_umap, cand.Probabilities, [cand_atom_num, 1]);

            marg_num = length(obj.MarginalWeights);

            coup_cell = cell(marg_num, 1);
            discreteOT_cost_list = zeros(marg_num, 1);

            for marg_id = 1:marg_num
                marg_atom_num = obj.Storage.MargVertNumList(marg_id);

                if marg_id == cand_id
                    % simply store the coupling between the candidate measure and the discretized marginal
                    coup_cell{marg_id} = sparse( ...
                        cand_umap, ...
                        cand.VertexIndices, ...
                        cand.Probabilities, ...
                        cand_atom_num, marg_atom_num);
                    discreteOT_cost_list(marg_id) = 0;

                    continue;
                end

                dual_sol = obj.Runtime.DualSolution{marg_id};
                [~, uind, umap] = unique(round(dual_sol.QualityAtoms, 6), 'rows', 'stable');
                quality_atoms = dual_sol.QualityAtoms(uind, :);
                quality_atom_num = size(quality_atoms, 1);

                % construct the sparse matrix representing the coupling between the discrete measure on the quality space and the
                % discretized marginal
                q2m_coup_mat = sparse( ...
                    umap, ...
                    dual_sol.VertexIndices, ...
                    dual_sol.Probabilities, ...
                    quality_atom_num, marg_atom_num);

                % compute the optimal coupling between the candidate discrete measure and this discrete measure on the quality space
                dist_mat = pdist2(cand_atoms, quality_atoms, 'euclidean');
                quality_probs = accumarray(umap, dual_sol.Probabilities, [quality_atom_num, 1]);

                if size(dist_mat, 1) * size(dist_mat, 2) > 1e6
                    % if there are too many atoms in the discrete measures, a direct computation of discrete OT may cause memory 
                    % throttling; thus, we resort to a constraint generation scheme
                    cp_options = struct('display', false);
                    OT = OTDiscrete(cand_probs, quality_probs, dist_mat, cp_options);
                    [hcoup_indices, ~] = OT.generateHeuristicCoupling();
                    OT.run(hcoup_indices, 1e-6);
                    coup = OT.Runtime.DualSolution;
                    coup_atoms = coup.CoupIndices;
                    coup_probs = coup.Probabilities;
                    discreteOT_cost_list(marg_id) = -OT.Runtime.LSIP_LB;
                else
                    [coup_atoms, coup_probs, ...
                        discreteOT_cost_list(marg_id)] ...
                        = discrete_OT( ...
                        cand_probs, ...
                        quality_probs, ...
                        dist_mat);
                end
                
                % construct the sparse matrix representing the coupling between the candidate discrete measure and this discrete
                % measure on the quality space
                q2q_coup_mat = sparse( ...
                    coup_atoms(:, 1), ...
                    coup_atoms(:, 2), ...
                    coup_probs, ...
                    cand_atom_num, ...
                    quality_atom_num);

                % the coupling between the discrete candidate measure and the discretized marginal is formed by composing the two
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

            % store the computed discrete OT costs in the runtime environment
            obj.Runtime.DiscreteOTCosts = discreteOT_cost_list;

            % store the computed couplings in the runtime environment
            obj.Runtime.DiscreteCouplings = coup_cell;
            obj.Runtime.DiscreteCouplingsComputed = true;
        end

        function samps = doRandSampleFromPartialReassembly(obj, ...
                marg_to_reassemble, ...
                samp_num, ...
                rand_stream)
            % Generate independent random samples from a partial reassembly where only couplings of certain marginals are sampled
            % Inputs:
            %   marg_to_reassemble: vector containing the indices of the marginals whose corresponding couplings will be sampled
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling
            % Outputs:
            %   samps: struct with the following fields:
            %       DiscreteIndices: matrix containing the indices of the atoms in the discretimized marginals; each row represents a 
            %       sample and each column represents a marginal
            %       DiscretePoints: cell array containing the coordinates of samples from the coupling of the discretized marginals
            %       ContinuousPoints: cell array containing the coordinates of samples from the continuous marginals (only those
            %       included in marg_to_reassmeble will be sampled)

            if ~isfield(obj.Runtime, 'DualSolution') || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            % first make sure that the semi-discrete OT problems are solved
            if ~isfield(obj.Storage, 'OTComputed') || ~obj.Storage.OTComputed
                obj.performReassembly();
            end

            % then make sure that the couplings between the candidate discrete measure and the discretized marginals are computed
            if ~isfield(obj.Runtime, 'DiscreteCouplingsComputed') || ~obj.Runtime.DiscreteCouplingsComputed
                obj.prepareDiscreteOT();
            end

            marg_num = length(obj.MarginalWeights);
            cand = obj.Runtime.DiscreteQualityCandidate;
            cand_atom_num = length(cand.Probabilities);

            % generate random indices of the atoms in the candidate discrete measure according to the probabilities
            disc_atom_index_samps = randsample(rand_stream, cand_atom_num, samp_num, true, cand.Probabilities);

            samps = struct;
            samps.DiscreteQualities = cand.Atoms(disc_atom_index_samps, :);
            samps.DiscreteIndices = zeros(samp_num, marg_num);
            samps.DiscretePoints = cell(marg_num, 1);
            samps.ContinuousPoints = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};
                marg_atoms = marg.SimplicialTestFuncs.Vertices;

                % sample the indices of atoms in the discretized marginal
                for cand_atom_id = 1:cand_atom_num
                    fill_list = disc_atom_index_samps == cand_atom_id;

                    [cond_atoms, ~, cond_probs] = find(obj.Runtime.DiscreteCouplings{marg_id}(cand_atom_id, :)');
                    
                    if length(cond_probs) > 1
                        % normalize the conditional probabilities
                        cond_probs = cond_probs / sum(cond_probs);
    
                        cond_samp_num = sum(fill_list);
                        samps.DiscreteIndices(fill_list, marg_id) = randsample(rand_stream, cond_atoms, cond_samp_num, ....
                            true, cond_probs);
                    else
                        % set all samples to be equal to that atom
                        samps.DiscreteIndices(fill_list, marg_id) = cond_atoms;
                    end
                end

                % store the coordinates of the samples with discrete marginals
                samps.DiscretePoints{marg_id} = marg_atoms(samps.DiscreteIndices(:, marg_id), :);

                if ~ismember(marg_id, marg_to_reassemble)
                    % skip those continuous marginals that are not required to be sampled
                    continue;
                end

                samps.ContinuousPoints{marg_id} = zeros(samp_num, 2);

                % count the number of samples coupled with each of the atoms in the discretized marginal
                atom_num_list = accumarray(samps.DiscreteIndices(:, marg_id), 1, [size(marg_atoms, 1), 1]);

                % retrieve the mapping from atoms with negligible probabilities to the remaining atoms
                filtered_atom_list = obj.Storage.OT.filtered_atom_list{marg_id};

                % generate from the conditional distributions
                cont_samp_cell = marg.conditionalRandSample(filtered_atom_list, atom_num_list, rand_stream);

                % fill in the coupled samples from the continuous marginals
                for atom_id = 1:length(atom_num_list)
                    samps.ContinuousPoints{marg_id}(samps.DiscreteIndices(:, marg_id) == atom_id, :) = cont_samp_cell{atom_id};
                end
            end
        end

        function samps = randSampleFromW2Couplings(obj, ...
                samp_num, rand_stream)
            % Generate independent random samples from a probability measure that glues together the Wasserstein-2 optimal couplings 
            % between the optimal discrete quality measure and the marginals.
            % Inputs: 
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling
            % Outputs:
            %   samps: cell array containing coupled samples from the marginals where each cell corresponds to a marginal

            % make sure that the semi-discrete Wasserstein-2 OT problems are solved 
            if ~isfield(obj.Runtime, 'W2OTComputed') || ~obj.Runtime.W2OTComputed
                obj.computeW2OptimalCouplings();
            end

            marg_num = length(obj.MarginalWeights);
            n = length(obj.Runtime.W2OT.DiscMeas.Probabilities);
            
            % generate random indices of the atoms according to the probabilities
            disc_atom_index_samps = randsample(rand_stream, ...
                n, ...
                samp_num, ...
                true, ...
                obj.Runtime.W2OT.DiscMeas.Probabilities);

            samps = cell(marg_num, 1);

            for marg_id = 1:marg_num
                marg = obj.Marginals{marg_id};

                samps{marg_id} = zeros(samp_num, 2);

                % count the number of samples coupled with each of the atoms in the discretized marginal
                atom_num_list = accumarray(disc_atom_index_samps, 1, [n, 1]);

                % generate from the conditional distributions
                cond_samp_cell = marg.conditionalRandSampleFromW2Coupling((1:n)', atom_num_list, rand_stream);

                % fill in the coupled samples from the continuous marginals
                for atom_id = 1:length(atom_num_list)
                    samps{marg_id}(disc_atom_index_samps == atom_id, :) = cond_samp_cell{atom_id};
                end
            end
        end

        function supported = checkIfWasserstein2OTSupported(obj)
            % Check whether all the marginals support Wasserstein-2 optimal transport.
            % Output:
            %   supported: boolean indicating whether all the marginals support Wasserstein-2 optimal transport

            supported = true;

            for marg_id = 1:length(obj.MarginalWeights)
                if ~isa(obj.Marginals{marg_id}, 'ProbMeas2DConvexPolytopeWithW2OT')
                    supported = false;
                    break;
                end
            end
        end

        function LP_options_runtime_new = handleLPErrors(obj, LP_result, LP_options_runtime, LP_trial_num) %#ok<INUSD>
            % Handle numerical errors that occurred while solving LP
            % Inputs: 
            %   LP_result: struct returned by the gurobi function representing the result from solving LP
            %   LP_options_runtime: struct containing the current options for solving LP
            %   LP_trial_num: integer representing the number of trials so far
            % Output:
            %   LP_options_runtime_new: struct containing the updated options for solving LP

            LP_options_runtime_new = LP_options_runtime;

            if LP_trial_num == 1
                % if the LP solver has failed once (reaching the time limit without converging), retry after  setting higher numeric 
                % focus, turning off presolve, and removing the existing bases
                LP_options_runtime_new.TimeLimit = LP_options_runtime.TimeLimit * 2;
                LP_options_runtime_new.NumericFocus = 3;
                LP_options_runtime_new.Quad = 1;
                LP_options_runtime_new.Presolve = 0;

                if isfield(obj.Runtime.CurrentLPModel, 'cbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'cbasis');
                end

                if isfield(obj.Runtime.CurrentLPModel, 'vbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'vbasis');
                end
            elseif LP_trial_num == 2 && isfield(LP_options_runtime, 'Method') && LP_options_runtime.Method == 2 ...
                && isfield(LP_options_runtime, 'Crossover') && LP_options_runtime.Crossover == 0
                % if the solver is using the barrier algorithm with crossover disabled, switch to the dual-simplex algorithm
                LP_options_runtime_new.Method = 1;
                LP_options_runtime_new = rmfield(LP_options_runtime_new, 'Crossover');
            elseif LP_trial_num == 3 && isfield(LP_options_runtime, 'Method') && LP_options_runtime.Method == 1
                % if the dual-simplex algorithm also fails, switch to the primal-simplex algorithm
                LP_options_runtime_new.Method = 0;
            else
                while true
                    % do nothing
                    warning('waiting for intervention...');
                    pause(10);
                    
                    continue_flag = false;

                    if continue_flag
                        break; %#ok<UNRCH>
                    end
                end
            end
        end
    end
end