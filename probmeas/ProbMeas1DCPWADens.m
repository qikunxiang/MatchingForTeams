classdef ProbMeas1DCPWADens < ProbMeas1DInterval
    % Class for probability measures with continuous piece-wise affine (CPWA) density function on a one-dimensional interval

    methods(Access = public)
        function obj = ProbMeas1DCPWADens(knots, dens_knots)
            % Constructor function
            % Inputs:
            %   knots: vector containing knots in the CPWA density function (must be in ascending order)
            %   dens_knots: vector containing the density at each knot (thus specifying the entire density function via interpolation)

            % call the superclass constructor to initialize the support
            obj@ProbMeas1DInterval(knots(1), knots(end));

            assert(all(diff(knots) > 0), 'knots are not in ascending order');

            if any(dens_knots < 0)
                error('the density is not non-negative');
            end

            % evaluate the integral of the unnormalized density within each sub-interval
            knot_diff = diff(knots);
            subint_integrals = (dens_knots(1:end - 1) + dens_knots(2:end)) .* knot_diff / 2;

            % calculate the normalizing constant and normalize the density
            norm_const = sum(subint_integrals);
            dens_knots = dens_knots / norm_const;
            obj.Dens = struct;
            obj.Dens.Knots = knots;
            obj.Dens.KnotDiff = knot_diff;
            obj.Dens.KnotDensities = dens_knots;
            obj.Dens.KnotDensityDiff = diff(dens_knots);
            obj.Dens.NormConst = norm_const;

            % prepare quantities for the evaluation of inverse CDF
            subint_integrals = subint_integrals / norm_const;
            obj.Dens.InvCDF = struct;
            obj.Dens.InvCDF.CDFKnots = [0; cumsum(subint_integrals)];

            % coefficients in the quadratic equations required for computing the inverse CDF
            obj.Dens.InvCDF.a = obj.Dens.KnotDiff .* obj.Dens.KnotDensityDiff;
            obj.Dens.InvCDF.b = obj.Dens.KnotDiff .* obj.Dens.KnotDensities(1:end - 1);
        end

        function dens = densityFunction(obj, x)
            % Compute the probability density function
            % Input: 
            %   x: vector containing the input points
            % Output: 
            %   dens: vector containing the computed densities

            inside = obj.checkIfInsideSupport(x);

            dens = (sum(min(max((x - obj.Dens.Knots(1:end - 1)') ./ obj.Dens.KnotDiff', 0), 1) .* obj.Dens.KnotDensityDiff', 2) ...
                + obj.Dens.KnotDensities(1)) .* inside;
        end

        function x = evaluateInverseCDF(obj, probs, batch_size)
            % Evaluate the inverse cumulative density function. Computation is done in batches if necessary.
            % Inputs:
            %   probs: vector containing inputs (must be between 0 and 1)
            %   batch_size: maximum number of inputs to be handled in the same batch (default is 1e4)
            % Output:
            %   x: vector containing outputs

            if ~exist('batch_size', 'var') || isempty(batch_size)
                batch_size = 1e4;
            end

            cdf_knots = obj.Dens.InvCDF.CDFKnots;
            knot_num = length(cdf_knots);

            input_num = size(probs, 1);
            batch_num = ceil(input_num / batch_size);
            vals_cell = cell(batch_num, 1);

            for batch_id = 1:batch_num
                probs_batch = probs(((batch_id - 1) * batch_size + 1):min(batch_id * batch_size, input_num));

                % compute which sub-interval the input belongs to
                int_id = min(sum(probs_batch - cdf_knots' >= 0, 2), knot_num - 1);

                % coefficients in the quadratic equations
                a = obj.Dens.InvCDF.a(int_id);
                b = obj.Dens.InvCDF.b(int_id);

                % divide into the cases where a == 0 and a ~= 0
                a_zero_list = abs(a) < 1e-8;

                % apply the quadratic formula when a ~= 0
                interp_weights = (sqrt(b .^ 2 + 2 * a .* (probs_batch - cdf_knots(int_id))) - b) ./ a;

                % apply the linear formula when a == 0
                interp_weights(a_zero_list) = (probs_batch(a_zero_list) - cdf_knots(int_id(a_zero_list))) ./ b(a_zero_list);

                vals_cell{batch_id} = obj.Dens.Knots(int_id) + interp_weights .* obj.Dens.KnotDiff(int_id);
            end

            x = min(max(vertcat(vals_cell{:}), obj.Supp.LowerBound), obj.Supp.UpperBound);
        end

        function cost = computeOTCost(obj)
            % Compute the optimal transport cost to the coupled discrete measure
            % Output:
            %   cost: the computed optimal transport cost

            if ~isfield(obj.OT, 'Coupled') || ~obj.OT.Coupled
                error('must set the coupled discrete measure first');
            end

            % retrieve the atoms in the discrete meaure
            atoms = obj.OT.DiscMeas.Atoms;
            [sorted_atoms, sorted_order] = sort(atoms, 'ascend');
            sorted_probs = obj.OT.DiscMeas.Probs(sorted_order);

            % compute intervals that are coupled with each atom by evaluating the inverse CDF
            itv_right = obj.evaluateInverseCDF(cumsum(sorted_probs));
            itv_left = [obj.Supp.LowerBound; itv_right(1:end - 1)];
            itv_right_dist = abs(itv_right - sorted_atoms);
            itv_left_dist = abs(itv_left - sorted_atoms);
            
            % the point in the interval with the least distance to the atom in the discrete measure
            itv_center = max(min(sorted_atoms, itv_right), itv_left);
            itv_center_dist = abs(itv_center - sorted_atoms);

            % the (discontinuous) piece-wise affine integrand function is described on each (possibly empty) interval
            pwalim_lb = [itv_left; itv_center];
            pwalim_ub = [itv_center; itv_right];
            pwalim_diff = pwalim_ub - pwalim_lb;
            pwaval_lb = [itv_left_dist; itv_center_dist];
            pwaval_ub = [itv_center_dist; itv_right_dist];
            pwaval_diff = pwaval_ub - pwaval_lb;

            % when an interval has length 0 it will result in 0 division error; thus, we change the denominator to 1
            pwalim_diff(pwalim_diff == 0) = 1;
            
            % prepare quantites used in the integration
            pdfknots_lb = obj.Dens.Knots(1:end - 1);
            pdfknots_ub = obj.Dens.Knots(2:end);
            pdfknots_diff = pdfknots_ub - pdfknots_lb;
            dens_lb = obj.Dens.KnotDensities(1:end - 1);
            dens_ub = obj.Dens.KnotDensities(2:end);
            dens_diff = dens_ub - dens_lb;

            % compute the lower and upper integration limits
            intlim_lb = max(pwalim_lb, pdfknots_lb');
            intlim_ub = max(min(pwalim_ub, pdfknots_ub'), intlim_lb);

            % evaluate the integrals
            D3 = intlim_ub .^ 3 - intlim_lb .^ 3;
            D2 = intlim_ub .^ 2 - intlim_lb .^ 2;
            D1 = intlim_ub .^ 1 - intlim_lb .^ 1;
            T1 = pwaval_lb .* pwalim_ub - pwaval_ub .* pwalim_lb;
            T2 = (dens_lb .* pdfknots_ub - dens_ub .* pdfknots_lb)';

            intval_mat = (dens_diff' .* pwaval_diff .* D3 / 3 + (dens_diff' .* T1 + pwaval_diff .* T2) .* D2 / 2 + T1 .* T2 .* D1) ...
                ./ pwalim_diff ./ pdfknots_diff';

            % sum up all integrals to compute the optimal transport cost
            cost = sum(sum(intval_mat));
        end
    end

    methods(Access = protected)
    
        function integrals = integrateSimplicialTestFuncs(obj, knots)
            % Compute the integrals of the test functions with respect to the probability measure 
            % Input: 
            %   knots: vector containing the knots in the simplicial test functions
            % Output:
            %   integrals: vector containing integrals of the test functions with respect to the probability measure; the number of 
            %   test functions is equal to the number of knots
            
            densities = obj.Dens.KnotDensities;
            density_diff = obj.Dens.KnotDensityDiff;
            density_knots = obj.Dens.Knots;
            density_knot_diff = diff(density_knots);
            knot_diff = diff(knots);

            % intersection of sub-intervals in the simplicial cover and sub-intervals in the density function
            int_lb = max(knots(1:end - 1), density_knots(1:end - 1)');
            int_ub = max(min(knots(2:end), density_knots(2:end)'), int_lb);

            D3 = int_ub .^ 3 - int_lb .^ 3;
            D2 = int_ub .^ 2 - int_lb .^ 2;
            D1 = int_ub - int_lb;
            T1 = densities(1:end - 1) .* density_knots(2:end) - densities(2:end) .* density_knots(1:end - 1);

            int_integrals1 = (density_diff' .* D3 / 3 + (T1' - density_diff' .* knots(1:end - 1)) ...
                .* D2 / 2 - knots(1:end - 1) .* T1' .* D1) ./ density_knot_diff' ./ knot_diff;
            int_integrals_tot = (density_diff' .* D2 / 2 + T1' .* D1) ./ density_knot_diff';

            right_integrals = sum(int_integrals1, 2);
            left_integrals = sum(int_integrals_tot, 2) - right_integrals;

            integrals = [left_integrals; 0] + [0; right_integrals];
        end
    end
end