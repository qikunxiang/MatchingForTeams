classdef ProbMeas2DMixNorm < ProbMeas2DConvexPolytope ...
        & HasTractableQuadraticIntegrals
    % Class for mixture of bivariate normal probability measures truncated to a union of triangles. The triangles need to satisfy the 
    % following condition: if two triangles A and B are not disjoint, then the intersection of A and B needs to be a face of both A and
    % B, i.e., the intersection of A and B is either an edge of A or a vertex of A, and is also either an edge of B or a vertex of B. 

    methods(Access = public)
        function obj = ProbMeas2DMixNorm(vertices, triangles, mixnorm)
            % Constructor function
            % Inputs:
            %   vertices: two-column matrix containing the vertices of all triangles
            %   triangles: three-column matrix indicating the triangulation
            %   mixnorm: struct with fields weights and components where weights is a vector containing positive weights of the mixture
            %   component and components is a cell array containing mean vectors and covariant matrices

            % call the superclass constructor to initialize the support
            obj@ProbMeas2DConvexPolytope(vertices);

            % check if there are duplicate vertices
            assert(size(unique(round(vertices, 6), 'rows'), 1) == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not appear in any triangle)
            assert(length(unique(triangles(:))) == size(vertices, 1), 'there are redundant vertices');

            % check that the triangles are non-empty (have non-zero area) and satisfy the condition
            [empty, is_simpcover] = check_triangulation_simplicial_cover(vertices, triangles);

            assert(all(~empty), 'there are empty triangles');
            assert(is_simpcover, 'the condition of simplicial cover is violated');

            obj.Supp.TriangleVertices = vertices;
            obj.Supp.Triangles = triangles;

            % number of mixture components
            comp_num = length(mixnorm.weights);

            % check the weights of the mixture components
            assert(all(mixnorm.weights > 0) && abs(sum(mixnorm.weights) - 1) < 1e-14, 'mixture weights must be positive and sum to 1');

            assert(length(mixnorm.components) == comp_num, 'mixture components mis-specified');

            tri_num = size(triangles, 1);

            % evaluate the integral of the unnormalized density within each triangular region
            tri_integral = zeros(tri_num, comp_num);
            
            % evaluate the integrals of x1, x2, x1^2, x2^2, and x1*x2 with respect to the mixture of bivariate Gaussian measure within 
            % each triangular region
            tri_lin_integral_x1 = zeros(tri_num, comp_num);
            tri_lin_integral_x2 = zeros(tri_num, comp_num);
            tri_quad_integral_x1sq = zeros(tri_num, comp_num);
            tri_quad_integral_x2sq = zeros(tri_num, comp_num);
            tri_quad_integral_x1x2 = zeros(tri_num, comp_num);

            for comp_id = 1:comp_num
                mean_vec = mixnorm.components{comp_id}.mean_vec;
                cov_mat = mixnorm.components{comp_id}.cov_mat;

                for tri_id = 1:tri_num
                    tri_verts = vertices(triangles(tri_id, :), :)';
                    
                    int_vals = ...
                        ProbMeas2DMixNorm.gaussianIntegralTriangle( ...
                        mean_vec, cov_mat, tri_verts, ...
                        [0, 0, 0, 0, 0, 1; ...
                        0, 0, 0, 1, 0, 0; ...
                        0, 0, 0, 0, 1, 0; ...
                        1, 0, 0, 0, 0, 0; ...
                        0, 1, 0, 0, 0, 0; ...
                        0, 0, 1, 0, 0, 0]', ...
                        [true; false; false; false; false; false]);
                    tri_integral(tri_id, comp_id) = int_vals(1);
                    tri_lin_integral_x1(tri_id, comp_id) = int_vals(2);
                    tri_lin_integral_x2(tri_id, comp_id) = int_vals(3);
                    tri_quad_integral_x1sq(tri_id, comp_id) = int_vals(4);
                    tri_quad_integral_x2sq(tri_id, comp_id) = int_vals(5);
                    tri_quad_integral_x1x2(tri_id, comp_id) = int_vals(6);
                end
            end

            % calculate the normalizing constant and normalize the density
            norm_const = sum(tri_integral * mixnorm.weights);

            % calculate the expected value of x1, x2, x1^2, x2^2, x1*x2
            obj.FirstMomentVec = zeros(2, 1);
            obj.FirstMomentVec(1) = sum(tri_lin_integral_x1 * mixnorm.weights) / norm_const;
            obj.FirstMomentVec(2) = sum(tri_lin_integral_x2 * mixnorm.weights) / norm_const;
            obj.SecondMomentMat = zeros(2, 2);
            obj.SecondMomentMat(1, 1) = sum(tri_quad_integral_x1sq * mixnorm.weights) / norm_const;
            obj.SecondMomentMat(2, 2) = sum(tri_quad_integral_x2sq * mixnorm.weights) / norm_const;
            obj.SecondMomentMat(1, 2) = sum(tri_quad_integral_x1x2 * mixnorm.weights) / norm_const;
            obj.SecondMomentMat(2, 1) = obj.SecondMomentMat(1, 2);

            % throughout this class, we evaluate densities and integrals of the non-truncated mixture of bivariate Gaussian and then
            % scale by the normalizing constant
            obj.Dens = struct('MixWeights', mixnorm.weights, 'NormConst', norm_const);

            obj.Dens.MixComponents = cell(comp_num, 1);

            for comp_id = 1:comp_num
                obj.Dens.MixComponents{comp_id} = struct( ...
                    'MeanVec', mixnorm.components{comp_id}.mean_vec, ...
                    'CovMat', mixnorm.components{comp_id}.cov_mat);
                obj.Dens.MixComponents{comp_id}.CovMatChol = chol(obj.Dens.MixComponents{comp_id}.CovMat);
            end

            obj.Dens.EdgeIneq = struct;
            
            ineq_w_cell = cell(tri_num, 1);
            ineq_b_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_vert_id = triangles(tri_id, :)';

                tri_vert_mat = vertices(tri_vert_id, :);

                % sort the three vertices counterclockwise
                ccorder = convhull(tri_vert_mat(:, 1), tri_vert_mat(:, 2), 'Simplify', true);

                % there are now 4 rows, the first and the last are identical
                tri_vert_cc = tri_vert_mat(ccorder, :);

                tri_vert1 = tri_vert_cc(1:3, :);
                tri_vert2 = tri_vert_cc(2:4, :);

                ineq_w_cell{tri_id} = [tri_vert2(:, 2) - tri_vert1(:, 2), tri_vert1(:, 1) - tri_vert2(:, 1)];
                ineq_b_cell{tri_id} = tri_vert2(:, 2) .* tri_vert1(:, 1) - tri_vert1(:, 2) .* tri_vert2(:, 1);

                % correct numerical inaccuracies and make sure that the vertices themselves are contained in the triangular region
                ineq_b_cell{tri_id} = max(ineq_b_cell{tri_id}, max(ineq_w_cell{tri_id} * tri_vert1', [], 2));
            end
            
            obj.Dens.EdgeIneq.w = vertcat(ineq_w_cell{:});
            obj.Dens.EdgeIneq.b = vertcat(ineq_b_cell{:});

            % prepare information needed for rejection sampling; since we use the non-truncated mixture of bivariate Gaussians as the
            % proposal, the multiplier is simply the inverse of the normalization constant
            obj.RejSamp = struct;
            obj.RejSamp.Multiplier = 1 / obj.Dens.NormConst;
        end

        function dens = densityFunction(obj, pts)
            % Compute the probability density function
            % Input: 
            %   pts: two-column matrix containing the input points
            % Output: 
            %   dens: vector containing the computed densities

            input_num = size(pts, 1);
            comp_num = length(obj.Dens.MixWeights);
            components = obj.Dens.MixComponents;

            dens_mat = zeros(input_num, comp_num);

            for comp_id = 1:comp_num
                % compute the density with respect to each mixture component
                dens_mat(:, comp_id) = mvnpdf(pts, components{comp_id}.MeanVec', components{comp_id}.CovMat);
            end

            inside = obj.checkIfInsideSupport(pts);

            dens = (dens_mat * obj.Dens.MixWeights) .* inside / obj.Dens.NormConst;
        end
    
        function inside_mat = checkIfInsideTriangularRegion(obj, pts)
            % Check if each of the input points is inside each of the triangular regions
            % Input: 
            %   pts: two-column matrix containing the inputs
            % Output:
            %   inside_mat: logical matrix where each row represents an input point and each column represents a triangle region

            region_num = size(obj.Supp.Triangles, 1);
            input_num = size(pts, 1);

            region_ineq_mat = obj.Dens.EdgeIneq.w * pts' - obj.Dens.EdgeIneq.b <= ProbMeas2DConvexPolytope.INSIDE_TOLERANCE;
            inside_mat = reshape(all(reshape(region_ineq_mat(:), 3, region_num * input_num), 1)', region_num, input_num)';
        end

        function pts_san = sanitizePoints(obj, pts, tolerance)
            % Fix small numerical inaccuracies in the points such that points that are slightly outside the support are moved into the 
            % support
            % Inputs: 
            %   pts: two-column matrix containing the input points
            %   tolerance: numerical tolerance values for constraint violations (default is 1e-6)
            % Output:
            %   pts_san: two-column matrix containing the sanitized points

            if ~exist('tolerance', 'var') || isempty(tolerance)
                tolerance = 1e-6;
            end

            region_num = size(obj.Supp.Triangles, 1);
            input_num = size(pts, 1);
            constr_vio = obj.Dens.EdgeIneq.w * pts' - obj.Dens.EdgeIneq.b;

            vio_mat = reshape(max(reshape(constr_vio(:), 3, region_num * input_num), [], 1)', region_num, input_num)';

            if max(min(vio_mat, [], 2)) > tolerance
                error('constraint violations exceed the tolerance');
            end

            pts_san = pts;

            for input_id = 1:input_num
                % find the triangle that the point is closest to
                [min_vio, min_id] = min(vio_mat(input_id, :));

                if min_vio <= 0
                    continue;
                end

                pt = pts(input_id, :)';
                tri = obj.Supp.Triangles(min_id, :);
                tri_verts = obj.Supp.TriangleVertices(tri, :);

                % represent the point as a convex combination of the vertices of the triangle
                coefs = [tri_verts'; ones(1, 3)] \ [pt; 1];

                % make all the coefficients positive
                coefs = max(coefs, 1e-6);

                % normalize the coefficients to sum to 1
                coefs = coefs / sum(coefs);

                % update the point
                pts_san(input_id, :) = tri_verts' * coefs;
            end
        end
    end

    methods(Access = protected)
        function inside = doCheckIfInsideSupport(obj, pts)
            % Check if the input points are inside the support of the probability measure
            % Input:
            %   pts: two-column matrix containing the input points
            % Output: 
            %   inside: boolean vector indicating whether each point is inside the support

            tri_num = size(obj.Supp.Triangles, 1);
            input_num = size(pts, 1);

            tri_ineq_mat = obj.Dens.EdgeIneq.w * pts' - obj.Dens.EdgeIneq.b <= ProbMeas2DConvexPolytope.INSIDE_TOLERANCE;
            tri_inside_mat = reshape(all(reshape(tri_ineq_mat(:), 3, tri_num * input_num), 1)', tri_num, input_num);
            inside = any(tri_inside_mat, 1)';
        end

        function raw_samps = sampleFromProposalDistribution(obj, samp_num, rand_stream)
            % Generate random samples from a proposal distribution for rejection sampling. For this probability measure, the proposal 
            % distribution is the non-truncated mixture of bivariate Gaussians with the same parameters.
            % Inputs:
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling

            comp_weights = obj.Dens.MixWeights;
            components = obj.Dens.MixComponents;
            comp_num = length(comp_weights);

            % first generate the component labels
            label_list = randsample(rand_stream, comp_num, samp_num, true, comp_weights);

            raw_samps = zeros(samp_num, 2);

            for comp_id = 1:comp_num
                component = components{comp_id};
                comp_id_list = label_list == comp_id;

                raw_samps(comp_id_list, :) = randn(rand_stream, sum(comp_id_list), 2) * component.CovMatChol + component.MeanVec';
            end
        end

        function accept_probs = rejSampAcceptProbFunction(obj, raw_samps)
            % Computes the acceptance probability when using rejection sampling to generate independent samples from the probability
            % measure
            % Input:
            %   raw_samps: two-column matrix containing the raw samples generated from the proposal distribution before rejection
            % Output:
            %   accept_probs: vector containing the computed acceptance probability

            accept_probs = obj.checkIfInsideSupport(raw_samps);
        end

        function prep = prepareOptimalTransport(obj, atoms, probs, angle_num)
            % Prepare quantities for the computation of semi-discrete optimal transport.
            % Inputs:
            %   atoms: two-column matrix containing the locations of the atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in the discrete measure
            %   angle_num: number of angles used in the integration

            prep = prepareOptimalTransport@ProbMeas2DConvexPolytope(obj, atoms, probs, angle_num);

            angle_list = prep.angle_list;
            angle_num = length(angle_list);

            comp_num = length(obj.Dens.MixWeights);
            atom_num = size(atoms, 1);
            tri_num = size(obj.Supp.Triangles, 1);
            prep.radial = cell(atom_num, 1);

            % coefficients in the left-hand side of the inequalities characterizing the triangles with repect to the distance away from
            % the atom; 
            % if this value is positive, it means that the ray is "leaving" the closed half-space, i.e., when the radial distance is
            % large enough the points on the ray will eventually be outside the closed half-space; 
            % if this value is non-positive, it means that the ray is "entering" or "remaining in" the closed half-space, i.e., when
            % the radial distance is large enough the points on the ray will eventually be inside the closed half-space
            ineq_coef = obj.Dens.EdgeIneq.w * [cos(angle_list)'; sin(angle_list)'];

            for atom_id = 1:atom_num
                % differences between the left-hand side of the inequalities characterizing the triangles and the right-hand side 
                % evaluated at the atom; 
                % if this value is positive it means that the atom is outside the closed half-space;
                % if this value is non-positive it means that the atom is inside the closed half-space
                ineq_diff = obj.Dens.EdgeIneq.w * prep.atoms(atom_id, :)' - obj.Dens.EdgeIneq.b;
                
                % radial distance of the point of intersection of the line going through the atom and a boundary of a triangle
                ineq_isec = -ineq_diff ./ ineq_coef;

                % 0/0 may occur, which means that the atom is on the boundary of the closed half-space and the ray is parallel to the
                % boundary; in this case, we set the intersection point to be at infinity
                ineq_isec_nan = isnan(ineq_isec);
                ineq_isec(ineq_isec_nan) = inf;
                
                % lower boundaries of the intervals
                itvl_lb = zeros(tri_num * 3, angle_num);

                % upper boundaries of the intervals
                itvl_ub = inf(tri_num * 3, angle_num);

                % boolean matrix indicating whether each atom is outside each closed half-space
                ineq_diff_pos = ineq_diff > 0;

                % boolean matrix indicating whether each ray is "leaving" each closed half-space
                ineq_coef_pos = ineq_coef > 0;

                % boolean matrix indicating if an atom is inside a closed half-space and its ray is "leaving" the closed half-space
                mat_npos_pos = (~ineq_diff_pos) & ineq_coef_pos;

                % boolean matrix indicating if an atom is outside a closed half-space and its ray is "entering" the closed half-space
                mat_pos_npos = ineq_diff_pos & (~ineq_coef_pos);

                % boolean matrix indicating if an atom is outside a closed half-space and its ray is not "entering" the closed
                % half-space
                mat_pos_nneg = ineq_diff_pos & (ineq_coef >= 0);

                % if an atom is outside a closed half-space and its ray is "entering" the closed half-space, then the intersection of 
                % the ray and the closed half-space is a ray emitting from the point of intersection
                itvl_lb(mat_pos_npos) = ineq_isec(mat_pos_npos);

                % if an atom is inside a closed half-space and its ray is "leaving" the closed half-space, then the intersection of the
                % ray and the closed half-space is a line segment between the atom and the point of intersection
                itvl_ub(mat_npos_pos) = ineq_isec(mat_npos_pos);

                % if an atom is outside a closed half-space and its ray is not "entering" the closed half-space, then the intersection
                % of the ray and the closed half-space is empty
                itvl_ub(mat_pos_nneg) = 0;

                % taking the intersection of the three intervals to get the intersection of the ray and the triangle
                tri_itvl_lb = reshape(max(reshape(itvl_lb(:), 3, tri_num * angle_num), [], 1), tri_num, angle_num);
                tri_itvl_ub = reshape(min(reshape(itvl_ub(:), 3, tri_num * angle_num), [], 1), tri_num, angle_num);

                % for intervals that are empty (or have empty interior), we set both their lower bounds and upper bounds to 0
                tri_itvl_empty = tri_itvl_lb >= tri_itvl_ub;
                tri_itvl_lb(tri_itvl_empty) = 0;
                tri_itvl_ub(tri_itvl_empty) = 0;


                % prepare quantities used to compute the integrals

                % unit vectors in each direction
                unit_v = [cos(angle_list), sin(angle_list)];

                % multiplicative coefficients and intercepts of the radial distance
                rho_coef = zeros(comp_num, angle_num);
                rho_intercept = zeros(comp_num, angle_num);

                % the multiplicative constant when evaluating the integrals of r -> r^2 * p(r); each row corresponds to a mixture 
                % component, each column corresponds to an angle
                quad_mult = zeros(comp_num, angle_num);

                % the multiplicative constant in front of the term involving normcdf
                quad_cdf_mult = zeros(comp_num, angle_num);

                % the multiplicative constant when evaluating the integrals of r -> r * p(r); each row corresponds to a mixture 
                % component, each column corresponds to an angle
                lin_mult = zeros(comp_num, angle_num);

                for comp_id = 1:comp_num
                    mean_vec = obj.Dens.MixComponents{comp_id}.MeanVec;
                    cov_mat = obj.Dens.MixComponents{comp_id}.CovMat;

                    diff_vec = atoms(atom_id, :)' - mean_vec;
                    a1 = sum(unit_v' .* (cov_mat \ unit_v'), 1)';
                    a2 = sum(unit_v' .* (cov_mat \ diff_vec), 1)';
                    a3 = diff_vec' * (cov_mat \ diff_vec);
                    a1_sqrt = sqrt(a1);
                    a2_over_a1_sqrt = a2 ./ a1_sqrt;

                    rho_coef(comp_id, :) = a1_sqrt';
                    rho_intercept(comp_id, :) = a2_over_a1_sqrt';

                    mult_const = 1 / sqrt(2 * pi) / sqrt(det(cov_mat)) * exp(-1/2 * (a3 - a2.^2 ./ a1)) * obj.Dens.MixWeights(comp_id);

                    quad_mult(comp_id, :) = (mult_const ./ a1 ./ a1_sqrt)';
                    quad_cdf_mult(comp_id, :) = (a2.^2 ./ a1 + 1)';
                    lin_mult(comp_id, :) = (-mult_const ./ a1)';
                end

                prep.radial{atom_id} = struct( ...
                    'interval_lb', tri_itvl_lb, ...
                    'interval_ub', tri_itvl_ub, ...
                    'radial_coefficient', rho_coef, ...
                    'radial_intercept', rho_intercept, ...
                    'intquad_mult', quad_mult / obj.Dens.NormConst, ...
                    'intquad_cdf_mult', quad_cdf_mult, ...
                    'intlin_mult', lin_mult / obj.Dens.NormConst);
            end
        end

        function [val1, val2] = computeInnerIntegral(obj, ...
                atom_id, prep, batch_indices, upper_limits)
            % Evaluates inner integrals along a list of angles in the polar coordinate system. The integrands are r -> r * p(r) and 
            % r -> r^2 * p(r) where p(r) denotes the probability density function.
            % Inputs:
            %   atom_id: the index of the atom (which is treated as the origin)
            %   prep: struct returned by the function obj.prepareOptimalTransport
            %   batch_indices: vector indicating the indices of the pre-specified angles where the inner integral is evaluated
            %   upper_limits: the upper integral limits specified in a vector (the lower integral limits is always 0)
            % Outputs:
            %   val1: the integral of r -> r * p(r)
            %   val2: the integral of r -> r^2 * p(r)

            radial = prep.radial{atom_id};
            rho_coef = radial.radial_coefficient(:, batch_indices);
            rho_intercept = radial.radial_intercept(:, batch_indices);
            quad_mult = radial.intquad_mult(:, batch_indices);
            quad_cdf_mult = radial.intquad_cdf_mult(:, batch_indices);
            lin_mult = radial.intlin_mult(:, batch_indices);
            
            % intersect the intervals (representing line segments in every direction) with (0, upper_limits), and making sure that 
            % empty intervals are represented with equal lower and upper limits
            itg_ub = min(radial.interval_ub(:, batch_indices), upper_limits');
            itg_lb = min(radial.interval_lb(:, batch_indices), itg_ub);

            % evaluate the integrals
            comp_num = length(obj.Dens.MixWeights);
            [intv_num, angle_num] = size(itg_ub);

            intquad_mat = zeros(intv_num, angle_num);
            intlin_mat = zeros(intv_num, angle_num);

            for comp_id = 1:comp_num
                rho_coef_comp = rho_coef(comp_id, :);
                rho_intercept_comp = rho_intercept(comp_id, :);

                % scale and shift the upper bounds and lower bounds
                ub_ss = rho_coef_comp .* itg_ub + rho_intercept_comp;
                lb_ss = rho_coef_comp .* itg_lb + rho_intercept_comp;

                ub_cdf = normcdf(ub_ss);
                lb_cdf = normcdf(lb_ss);
                ub_pdf = normpdf(ub_ss);
                lb_pdf = normpdf(lb_ss);

                ub_ss_alt = -rho_coef_comp .* itg_ub + rho_intercept_comp;
                lb_ss_alt = -rho_coef_comp .* itg_lb + rho_intercept_comp;

                % compute the integrals with each mixture component
                intquad_mat = intquad_mat + (quad_cdf_mult(comp_id, :) .* (ub_cdf - lb_cdf) ...
                    + (ub_ss_alt .* ub_pdf) - (lb_ss_alt .* lb_pdf)) .* quad_mult(comp_id, :);
                intlin_mat = intlin_mat + ((ub_pdf - lb_pdf) + rho_intercept_comp .* (ub_cdf - lb_cdf)) .* lin_mult(comp_id, :);
            end

            val1 = sum(intlin_mat, 1)';
            val2 = sum(intquad_mat, 1)';
        end
     
        function postProcessOptimalTransport(obj, pp_angle_num)
            % Post-processing after computing semi-discrete optimal transport including computing the contour of each cell in the
            % Laguerre diagram and preparing for conditional rejection sampling
            % Input:
            %   pp_angle_num: number of angles used in post-processing

            atom_num = size(obj.OT.DiscMeas.Atoms, 1);
            
            % cell array containing information about rejection sampling from the conditional distributions (conditional on each atom)
            obj.OT.CondRejSamp = cell(atom_num, 1);

            hyp_a = (obj.OT.Weights - obj.OT.Weights') / 2;

            % take an evenly-spaced 10000 angles for this step
            prep = obj.prepareOptimalTransport(obj.OT.DiscMeas.Atoms, obj.OT.DiscMeas.Probs, pp_angle_num);
            angle_samp = prep.angle_list;

            for atom_id = 1:atom_num
                [rho_temp, hyp_full] = ProbMeas2DConvexPolytope.hyperbolic_Dirichlet_tessellation_polar( ...
                    angle_samp, ...
                    hyp_a(:, atom_id), ...
                    prep.hyp.c(:, atom_id), ...
                    prep.hyp.intercept(:, atom_id), ...
                    prep.lin{atom_id}.c, ...
                    prep.lin{atom_id}.intercept, ...
                    true);

                radial = prep.radial{atom_id};

                obj.OT.CondRejSamp{atom_id} = struct;
                itvl_ub_max = max(radial.interval_ub, [], 1)';

                % take the distance to be the smaller of the distance calculated from the Laguerre diagram and the intersection with 
                % the triangles
                rho_i = min(rho_temp, itvl_ub_max);

                % add 10% of the standard deviation to make sure that the circle covers the whole region
                std_rho = std(rho_i);
                radius = max(rho_i) + std_rho * 0.1;

                obj.OT.CondRejSamp{atom_id}.Radius = radius;
                obj.OT.CondRejSamp{atom_id}.NeighborAtoms = [find(any(hyp_full < radius, 1))'; atom_id];
                
                % conservative estimate of the maximum density in the cell in each direction
                comp_num = length(obj.Dens.MixWeights);
                circ_center = obj.OT.DiscMeas.Atoms(atom_id, :)';
                circum_pts =  circ_center' + radius * [cos(angle_samp), sin(angle_samp)];
                max_density = 0;

                for comp_id = 1:comp_num
                    mean_vec = obj.Dens.MixComponents{comp_id}.MeanVec;
                    cov_mat = obj.Dens.MixComponents{comp_id}.CovMat;
                    cov_chol = obj.Dens.MixComponents{comp_id}.CovMatChol;

                    if sum((circ_center - mean_vec).^2) <= radius^2
                        % if the mean vector of the mixture component is inside the circle, the maximum is attained at the mean vector
                        min_dist_sq = 0;
                    else
                        % otherwise, choose the point that minimizes the Mahalanobis distance to the mean vector in a grid on the 
                        % circumference, then minus a 2% margin
                        min_dist_sq = min(sum(((circum_pts - mean_vec') / cov_chol) .^ 2, 2)) * 0.98;
                    end

                    % evaluate the estimated un-normalized density
                    max_density = max_density + 1 / (2 * pi) / sqrt(det(cov_mat)) * exp(-1/2 * min_dist_sq) ...
                        * obj.Dens.MixWeights(comp_id);
                end

                % first divide by the normalizing constant to get the estimated maximum marginal density, then divide by the
                % probability of the atom to get the estimated maximum conditional density
                max_density = max_density / obj.Dens.NormConst / obj.OT.DiscMeas.Probs(atom_id);

                obj.OT.CondRejSamp{atom_id}.MaxDensity = max_density;
                obj.OT.CondRejSamp{atom_id}.Multiplier = max_density * pi * radius ^ 2;
            end
        end

        function raw_samps = sampleFromCondProposalDistribution(obj, ...
                atom_id, samp_num, rand_stream)
            % Generate random samples from a proposal distribution for conditional rejection sampling. For this procedure, the proposal
            % distribution is uniform on a circle.
            % Inputs:
            %   atom_id: the index of the corresponding atom in the discrete measure
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object used for sampling

            % randomly generate points uniformly in the circle
            dists_from_center = sqrt(rand(rand_stream, samp_num, 1)) * obj.OT.CondRejSamp{atom_id}.Radius;
            angles = rand(rand_stream, samp_num, 1) * 2 * pi;

            % raw samples before the rejection step
            raw_samps = obj.OT.DiscMeas.Atoms(atom_id, :) + dists_from_center .* [cos(angles), sin(angles)];
        end

        function accept_probs = condRejSampAcceptProbFunction(obj, ...
                atom_id, raw_samps)
            % Compute the acceptance probability in conditional rejection sampling from a cell
            % Inputs: 
            %   atom_id: the index of the atom (corresponding to the cell that is being sampled)
            %   raw_samps: two-column matrix containing the raw samples before the rejection step 
            % Output:
            %   accept_probs: vector containing the computed acceptance probabilities

            raw_density = obj.densityFunction(raw_samps);

            condrejsamp = obj.OT.CondRejSamp{atom_id};
            
            radius = condrejsamp.Radius;
            multiplier = condrejsamp.Multiplier;
            neighbors = condrejsamp.NeighborAtoms;

            neighbor_atoms = obj.OT.DiscMeas.Atoms(neighbors, :);

            weighted_dist = sqrt((raw_samps(:, 1) - neighbor_atoms(:, 1)') .^ 2 + (raw_samps(:, 2) - neighbor_atoms(:, 2)') .^ 2) ...
                - obj.OT.Weights(neighbors)';

            accept_probs = raw_density .* all(weighted_dist(:, neighbors == atom_id) <= weighted_dist, 2) ...
                / (obj.OT.DiscMeas.Probs(atom_id) * multiplier / (pi * radius ^ 2));
        end
    
        function integrals = doSetSimplicialTestFuncs(obj, ...
                vertices, triangles_cell, check_simpcover)
            % Initialize the test functions with respect to a simplicial cover and compute the integrals of the test functions with
            % respect to the probability measure
            % Inputs: 
            %   vertices: two-column matrix containing the vertices used in the triangulation; the first rows must be identical to the
            %   vertices of the support
            %   triangles_cell: cell array where the number of cells is equal to the number of triangular regions; each cell contains a 
            %   three-column matrix representing the triangulation within a triangular region
            %   check_simpcover: check whether the given triangulation forms a simplicial cover (default is true)
            % Output:
            %   integrals: vector containing integrals of the test functions with respect to the probability measure; the number of 
            %   test functions is equal to the number of vertices in the triangulation

            if ~exist('check_simpcover', 'var') || isempty(check_simpcover)
                check_simpcover = true;
            end

            region_num = size(obj.Supp.Triangles, 1);
            assert(length(triangles_cell) == region_num, 'triangulation mis-specified');

            % check if the first vertices are identical to the vertices in
            % the triangular regions
            assert(max(max(abs(vertices(1:size(obj.Supp.TriangleVertices, 1), :) - obj.Supp.TriangleVertices))) < 1e-14, ...
                'the vertices do not begin with the vertices of the triangular regions');

            % check if there are duplicate vertices
            assert(size(unique(round(vertices, 6), 'rows'), 1) == size(vertices, 1), 'there are duplicate vertices');

            % check if there are redundant vertices (those that do not appear in any triangle)
            triangles_agg = vertcat(triangles_cell{:});
            assert(length(unique(triangles_agg(:))) == size(vertices, 1), 'there are redundant vertices');

            region_inside_mat = obj.checkIfInsideTriangularRegion(vertices);

            for region_id = 1:region_num
                % check if all vertices are contained in the triangular region
                region_vert_indices = unique(triangles_cell{region_id}(:));
                assert(all(region_inside_mat(region_vert_indices, region_id)), ...
                    'there are vertices used in the triangulation of this region that is outside of the region');
            end

            if check_simpcover
                % check if the triangulation forms a simplicial cover
                [empty, is_simpcover] = check_triangulation_simplicial_cover(vertices, triangles_agg);
    
                if any(empty)
                    error('some triangles are empty');
                end
    
                if ~is_simpcover
                    error('the triangulation does not form a simplicial cover');
                end
            end
            
            obj.SimplicialTestFuncs = struct;
            obj.SimplicialTestFuncs.Vertices = vertices;
            obj.SimplicialTestFuncs.Triangles = triangles_agg;

            vert_num = size(vertices, 1);
            tri_num = size(triangles_agg, 1);

            components = obj.Dens.MixComponents;
            comp_num = length(components);

            % evaluate the integral of three interpolation functions with respect to the probability measure within each triangular
            % region via an affine transformation which places the vertices on (0, 0), (1, 0), and (0, 1)
            tri_integrals = zeros(tri_num, 3);

            % compute the inverse transformation matrix which transform a coordinate into weights in the affine combination
            invtrans_cell = cell(tri_num, 1);

            for tri_id = 1:tri_num
                tri_verts = vertices(triangles_agg(tri_id, :), :);

                % compute the integrals of the three interpolation functions
                for comp_id = 1:comp_num
                    tri_integrals(tri_id, :) = tri_integrals(tri_id, :) + obj.Dens.MixWeights(comp_id) ...
                        * ProbMeas2DMixNorm.gaussianIntegralTriangle( ...
                        components{comp_id}.MeanVec, ...
                        components{comp_id}.CovMat, ...
                        tri_verts', ...
                        [0, 0, 0, -1, 0, 1; ...
                        0, 0, 0, 1, -1, 0; ...
                        0, 0, 0, 0, 1, 0]', ...
                        true(3, 1))';
                end

                % divide by the normalization constant
                tri_integrals(tri_id, :) = tri_integrals(tri_id, :) / obj.Dens.NormConst;

                % compute the inverse transformation matrix
                invtrans_cell{tri_id} = [tri_verts'; ones(1, 3)] \ eye(3);
            end

            integrals = accumarray(triangles_agg(:), tri_integrals(:), [vert_num, 1]);
            obj.SimplicialTestFuncs.Integrals = integrals;

            if abs(sum(integrals) - 1) > 1e-10
                error('integrals are incorrectly computed');
            end


            obj.SimplicialTestFuncs.InvTransMat = vertcat(invtrans_cell{:});
        end
    end

    methods(Static, Access = protected)
        function vals = gaussianIntegralTriangle(mean_vec, cov_mat, ...
                z_mat, polynomial_coefs, after_transform)
            % Integrate bivariate second-order polynomials inside a triangle over a bivariate Gaussian measure
            % Inputs:
            %   mean_vec: the mean vector (2-by-1) of the Gaussian measure
            %   cov_mat: the covariance matrix (2-by-2) of the Gaussian measure
            %   z_mat: 2-by-3 matrix containing the three vertices of the triangle as columns
            %   polynomial_coefs: matrix with six rows where each column represents a polynomial integrand and each row represents a
            %   coefficient in the polynomial: c1 * x1^2 + c2 * x2^2 + c3 * x1 * x2 + c4 * x1 + c5 * x2 + c6
            %   after_transform: logical vector indicating whether the coefficients of each polynomial is applicable to the variables 
            %   before transformation or after transforming to the triangle with vertices [0; 0], [1; 0], and [1; 1]
            % Output: 
            %   vals: vector containing the values of the integrals; each entry corresponds to a polynomial

            % compute the affine transformation that will transform the triangle formed by the columns of z_mat to the triangled with
            % vertices [0; 0], [1; 0], and [1; 1]
            trans_all = [0, 1, 1; 0, 0, 1] / [z_mat; ones(1, 3)];
            trans_mat = trans_all(:, 1:2);
            trans_shift = trans_all(:, 3);
            trans_mat_inv = inv(trans_mat);

            % compute the mean vector and the covariance matrix of the Gaussian measure after the affine transformation
            trans_mean_vec = trans_mat * mean_vec + trans_shift;
            trans_cov_mat = trans_mat * cov_mat * trans_mat';

            vals = zeros(size(polynomial_coefs, 2), 1);

            for poly_id = 1:size(polynomial_coefs, 2)
                if ~after_transform(poly_id)
                    % if the coefficients apply to the variables before transformation, we need to compute the coefficients after 
                    % transformation
                    coefs = polynomial_coefs(:, poly_id);
                    coefs_mat = [coefs(1), coefs(3)/2; coefs(3)/2, coefs(2)];
                    coefs_vec = coefs(4:5);
                    coefs_const = coefs(6);

                    tcoefs_mat = trans_mat_inv' * coefs_mat * trans_mat_inv;
                    tcoefs_vec = trans_mat_inv' * coefs_vec - 2 * tcoefs_mat * trans_shift;
                    tcoefs_const = trans_shift' * tcoefs_mat * trans_shift - coefs_vec' * trans_mat_inv * trans_shift ...
                        + coefs_const; %#ok<MINV>
                    coefs = [tcoefs_mat(1, 1); tcoefs_mat(2, 2); 2 * tcoefs_mat(1, 2); tcoefs_vec; tcoefs_const];
                else
                    % otherwise, we can directly integrate over the polynomial with the given coefficients
                    coefs = polynomial_coefs(:, poly_id);
                end

                outer_integrand = @(x1)(ProbMeas2DMixNorm.outerGaussianIntegrand(trans_mean_vec, trans_cov_mat, coefs, x1));

                vals(poly_id) = integral(outer_integrand, 0, 1) / sqrt(trans_cov_mat(1, 1));
            end
        end

        function vals = outerGaussianIntegrand(mean_vec, cov_mat, c, x1)
            % Evaluates the outer integrand when integrating a bivariate second-order polynomial inside a triangle over a bivariate
            % Gaussian measure
            % Inputs: 
            %   mean_vec: the mean vector of the transformed bivariate Gaussian
            %   cov_mat: the covariance matrix of the transformed bivariate Gaussian
            %   c: 6-by-1 vector containing the coefficients in the 
            %   polynomial: c1 * x1^2 + c2 * x2^2 + c3 * x1 * x2 + c4 * x1 + c5 * x2 + c6
            %   x1: value of the first variable where the inner integral is evaluated
            % Output:
            %   vals: the computed integrands

            % compute the conditional mean and variance of the second variable conditional on the first variable
            cond_m = mean_vec(2) + cov_mat(1, 2) / cov_mat(1, 1) ...
                * (x1 - mean_vec(1));
            cond_s = sqrt(cov_mat(2, 2) - cov_mat(1, 2)^2 / cov_mat(1, 1));

            vals = (c(2) * (cond_m.^2 + cond_s^2) + (c(3) * x1 + c(5)) .* cond_m + c(1) * x1.^2 + c(4) * x1 + c(6)) ...
                .* (normcdf((x1 - cond_m) / cond_s) - normcdf(- cond_m / cond_s)) + (c(2) * cond_s * cond_m ...
                + (c(3) * x1 + c(5)) * cond_s) .* normpdf(- cond_m / cond_s) - (c(2) * (cond_s * cond_m + cond_s * x1) ...
                + (c(3) * x1 + c(5)) * cond_s) .* normpdf((x1 - cond_m) / cond_s);

            vals = vals .* normpdf((x1 - mean_vec(1)) / sqrt(cov_mat(1, 1)));
        end
    end
end