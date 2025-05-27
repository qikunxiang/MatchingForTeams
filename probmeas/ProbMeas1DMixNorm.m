classdef ProbMeas1DMixNorm < ProbMeas1DInterval
    % Class for mixture of normal probability measures in 1D truncated to a
    % compact interval

    methods(Access = public)
        function obj = ProbMeas1DMixNorm(trunc_lb, trunc_ub, mixnorm)
            % Constructor function
            % Inputs:
            %   trunc_lb: lower end of the support interval
            %   trunc_ub: upper end of the support interval
            %   mixnorm: struct with fields weights, mean_list, and
            %   std_list where weights is a vector containing positive
            %   weights of the mixture component, mean_list is a vector
            %   containing the mean of the mixture components, and
            %   std_list is a vector containing the standard deviation of
            %   the mixture components

            % call the superclass constructor to initialize the support
            obj@ProbMeas1DInterval(trunc_lb, trunc_ub);

            comp_num = length(mixnorm.weights);
            obj.Dens.NumOfComponents = comp_num;

            assert(abs(sum(mixnorm.weights) - 1) < 1e-10, ...
                'the weights must sum up to 1');

            obj.Dens.Weights = mixnorm.weights / sum(mixnorm.weights);

            assert(length(mixnorm.mean_list) == comp_num, ...
                'list of mean mis-specified');
            assert(length(mixnorm.std_list) == comp_num, ...
                'list of standard deviation mis-specified');

            obj.Dens.MeanList = mixnorm.mean_list;
            obj.Dens.StdList = mixnorm.std_list;

            % the CDF of the untruncated normal mixture evaluated at the
            % lower end point of the support interval
            obj.Dens.UntruncatedCDFAtLB = obj.untruncatedCDF( ...
                obj.Supp.LowerBound);

            % the normalizing constant, i.e., the integral of the indicator
            % function of the support interval with respect to the
            % untruncated normal mixture
            obj.Dens.NormConst = obj.untruncatedCDF( ...
                obj.Supp.UpperBound) - obj.Dens.UntruncatedCDFAtLB;
        end

        function dens = densityFunction(obj, x)
            % Compute the probability density function
            % Input: 
            %   x: vector containing the input points
            % Output: 
            %   dens: vector containing the computed densities

            inside = obj.checkIfInsideSupport(x);

            dens = zeros(length(x), 1);

            for comp_id = 1:obj.Dens.NumOfComponents
                dens = dens + obj.Dens.Weights(comp_id) ...
                    * normpdf(x, obj.Dens.MeanList(comp_id), ...
                    obj.Dens.StdList(comp_id));
            end

            dens = (dens / obj.Dens.NormConst) .* inside;
        end

        function x = evaluateInverseCDF(obj, probs, tolerance)
            % Evaluate the inverse cumulative density function. The inverse
            % values are approximated via the bisection method up to a
            % given accuracy. 
            % Inputs:
            %   probs: vector containing inputs (must be between 0 and 1)
            %   tolerance: required accuracy for the approximation (default
            %   is 1e-8)
            % Output:
            %   x: vector containing outputs

            if ~exist('tolerance', 'var') || isempty(tolerance)
                tolerance = 1e-8;
            end

            input_num = length(probs);

            prob_rescaled = probs * obj.Dens.NormConst ...
                + obj.Dens.UntruncatedCDFAtLB;

            % since we always return the mid-points of the sub-intervals
            % from the final iteration, error will always be reduced by
            % half before returning
            iter_num = ceil(log2((obj.Supp.UpperBound ...
                - obj.Supp.LowerBound) / tolerance)) - 1;

            itv_lb = obj.Supp.LowerBound * ones(input_num, 1);
            itv_ub = obj.Supp.UpperBound * ones(input_num, 1);

            for iter_id = 1:iter_num
                % evaluate the untruncated CDF at the mid-point of the
                % interval
                itv_mid = (itv_lb + itv_ub) / 2;
                cdf_mid = obj.untruncatedCDF(itv_mid);

                % list of inputs whose inverse belong to the right
                % sub-interval
                right_subitv_list = prob_rescaled > cdf_mid;

                % update the sub-intervals
                itv_lb(right_subitv_list) = itv_mid(right_subitv_list);
                itv_ub(~right_subitv_list) = itv_mid(~right_subitv_list);
            end

            % return the mid-points
            x = (itv_lb + itv_ub) / 2;

            % when evaluating the inverse CDF with inputs close to 0/1, the
            % returned values might not be the lower/upper bound due to the
            % probabilities being "concentrated" away from the boundary; in
            % this case, we force the outputs to be equal to the
            % lower/upper bound
            x(probs <= 1e-15) = obj.Supp.LowerBound;
            x(probs >= 1 - 1e-15) = obj.Supp.UpperBound;
        end

        function cost = computeOTCost(obj)
            % Compute the optimal transport cost to the coupled discrete
            % measure
            % Output:
            %   cost: the computed optimal transport cost

            if ~isfield(obj.OT, 'Coupled') || ~obj.OT.Coupled
                error('must set the coupled discrete measure first');
            end

            % retrieve the atoms in the discrete meaure
            atoms = obj.OT.DiscMeas.Atoms;
            [sorted_atoms, sorted_order] = sort(atoms, 'ascend');
            sorted_probs = obj.OT.DiscMeas.Probs(sorted_order);

            % compute intervals that are coupled with each atom by
            % evaluating the inverse CDF
            itv_right = obj.evaluateInverseCDF(cumsum(sorted_probs));
            itv_left = [obj.Supp.LowerBound; itv_right(1:end - 1)];
            itv_right_dist = abs(itv_right - sorted_atoms);
            itv_left_dist = abs(itv_left - sorted_atoms);
            
            % the point in the interval with the least distance to the atom
            % in the discrete measure
            itv_center = max(min(sorted_atoms, itv_right), itv_left);
            itv_center_dist = abs(itv_center - sorted_atoms);

            % the (discontinuous) piece-wise affine integrand function is
            % described on each (possibly empty) interval
            pwalim_lb = [itv_left; itv_center];
            pwalim_ub = [itv_center; itv_right];
            pwalim_diff = pwalim_ub - pwalim_lb;
            pwaval_lb = [itv_left_dist; itv_center_dist];
            pwaval_ub = [itv_center_dist; itv_right_dist];
            pwaval_diff = pwaval_ub - pwaval_lb;

            % when an interval has length 0 it will result in 0 division
            % error; thus, we change the denominator to 1
            pwalim_diff(pwalim_diff == 0) = 1;
            
            % start to compute the integrals
            piece_integral = zeros(length(pwalim_lb), 1);

            for comp_id = 1:obj.Dens.NumOfComponents
                slopes = pwaval_diff ./ pwalim_diff;
                piece_integral = piece_integral ...
                    + obj.Dens.Weights(comp_id) ...
                    * obj.normalPartialExpectation( ...
                    obj.Dens.MeanList(comp_id), ...
                    obj.Dens.StdList(comp_id), ...
                    pwalim_lb, pwalim_ub, ...
                    slopes, pwaval_lb - slopes .* pwalim_lb);
            end

            % sum up all integrals to compute the optimal transport cost
            cost = sum(piece_integral) / obj.Dens.NormConst;
        end
    end

    methods(Access = protected)

        function probs = untruncatedCDF(obj, x)
            % Compute the untruncated cumulative distribution function,
            % i.e., the same normal mixture before being truncated to the
            % support interval
            % Input:
            %   x: vector containing the input points
            % Output:
            %   probs: vector containing the computed probabilities

            probs = zeros(length(x), 1);

            for comp_id = 1:obj.Dens.NumOfComponents
                probs = probs + obj.Dens.Weights(comp_id) ...
                    * normcdf(x, obj.Dens.MeanList(comp_id), ...
                    obj.Dens.StdList(comp_id));
            end
        end
    
        function integrals = integrateSimplicialTestFuncs(obj, knots)
            % Compute the integrals of the test functions with respect to
            % the probability measure 
            % Input: 
            %   knots: vector containing the knots in the simplicial test
            %   functions
            % Output:
            %   integrals: vector containing integrals of the test
            %   functions with respect to the probability measure; the
            %   number of test functions is equal to the number of knots
            
            knot_num = length(knots);
            knot_diff = diff(knots);

            % integrals of the piece to the left of each knot
            left_integrals = zeros(knot_num, 1);

            % integrals of the piece to the right of each knot
            right_integrals = zeros(knot_num, 1);

            for comp_id = 1:obj.Dens.NumOfComponents
                % the first knot does not have a piece to its left; for the
                % rest of the knots, integrate the affine function that
                % evaluates to 1 at this knot and evaluates to 0 at the
                % prior knot
                left_integrals(2:end) = left_integrals(2:end) ...
                    + obj.Dens.Weights(comp_id) ...
                    * obj.normalPartialExpectation( ...
                    obj.Dens.MeanList(comp_id), ...
                    obj.Dens.StdList(comp_id), ...
                    knots(1:end - 1), knots(2:end), ...
                    1 ./ knot_diff, -knots(1:end - 1) ./ knot_diff);

                % the last knot does not have a piece to its right; for the
                % rest of the knots, integrate the affine function that
                % evaluates to 1 at this knot and evaluates to 0 at the
                % next knot
                right_integrals(1:end - 1) = right_integrals(1:end - 1) ...
                    + obj.Dens.Weights(comp_id) ...
                    * obj.normalPartialExpectation( ...
                    obj.Dens.MeanList(comp_id), ...
                    obj.Dens.StdList(comp_id), ...
                    knots(1:end - 1), knots(2:end), ...
                    -1 ./ knot_diff, knots(2:end) ./ knot_diff);
            end

            integrals = (left_integrals + right_integrals) ...
                / obj.Dens.NormConst;
        end
    end

    methods (Static, Access = public)
        function vals = normalPartialExpectation(norm_mean, norm_std, ...
            interval_lb, interval_ub, slope, intercept)
            % Compute the partial expectation of a normal distribution of
            % the form E[(aX+b)I{k1,k2}] given the parameters of the normal
            % distribution and the values of a, b, k1, and k2.
            % Scalar-valued inputs will be extended to vectors with
            % appropriate lengths.
            % Inputs:
            %   norm_mean: the mean parameter of the normal distribution
            %   norm_std: the standard deviation parameter of the normal
            %   distribution
            %   interval_lb: the left end point of the interval (lower
            %   integration limit), i.e., k1
            %   interval_ub: the right end point of the interval (upper
            %   integration limit), i.e., k2
            %   slope: the coefficient of the linear part, i.e., a
            %   intercept: the intercept of the linaer part, i.e., b
            % Outputs:
            %   vals: the computed partial expectation

            if isscalar(interval_lb) && ~isscalar(interval_ub) > 1
                interval_lb = repmat(interval_lb, length(interval_ub), 1);
            elseif isscalar(interval_ub) && length(interval_lb) > 1
                interval_ub = repmat(interval_ub, length(interval_lb), 1);
            end

            interval_ub = max(interval_lb, interval_ub);

            % standardize the end points
            standardized = ([interval_lb, interval_ub] - norm_mean) ...
                ./ norm_std;

            cdf_mat = normcdf(standardized);

            if all(slope == 0)
                vals = intercept .* (cdf_mat(:, 2) - cdf_mat(:, 1));
                return;
            end

            exp_mat = exp(-0.5 * (standardized .^ 2));

            vals = (slope .* norm_mean + intercept) ...
                .* (cdf_mat(:, 2) - cdf_mat(:, 1)) + (1 / sqrt(2 * pi)) ...
                * slope .* norm_std .* (exp_mat(:, 1) - exp_mat(:, 2));
        end
    end
end