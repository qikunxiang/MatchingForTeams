function [WB_mean_vec, WB_cov_mat, min_cost, cost_iter] = ...
    WB_elliptical_fixedpoint(mean_vec_cell, cov_mat_cell, weights, tolerance)
% Compute the 2-Wasserstein barycenter of probability measures from a family of elliptical distributions via a fixed-point scheme.
% Inputs: 
%       mean_vec_cell: cell array containing the mean vectors of the input measures
%       cov_mat_cell: cell array containing the covariance matrices of the input measures
%       weights: the weights defining the 2-Wasserstein barycenter
%       tolerance: the algorithm terminates when the violation of the fixed-point criterion is below this tolerance (default is 1e-10)
% Outputs: 
%       WB_mean_vec: the mean vector of the computed 2-Wasserstein barycenter
%       WB_cov_mat: the covariance matrix of the computed 2-Wasserstein barycenter
%       min_cost: the minimized cost value
%       cost_iter: vector containing the cost value at each iteration

if ~exist('tolerance', 'var') || isempty(tolerance)
    tolerance = 1e-10;
end

input_num = length(mean_vec_cell);

assert(length(cov_mat_cell) == input_num || length(weights) == input_num, 'measures and weights mis-match');

dim = length(mean_vec_cell{1});

% the mean vector of the Wasserstein barycenter is simply the barycenter of the mean vectors of the input measures
WB_mean_vec = horzcat(mean_vec_cell{:}) * weights;

% compute the cost related to mean vectors and the trace of the covariance matrices
min_cost = 0;
for input_id = 1:input_num
    min_cost = min_cost + weights(input_id) * (sum((mean_vec_cell{input_id} - WB_mean_vec) .^ 2) + trace(cov_mat_cell{input_id})); 
end

S = eye(dim);
iter = 0;

% usually the iteration terminates after a few iterations
cost_iter = zeros(100, 1);

while true
    S_sqrt = sqrtm(S);
    T = zeros(dim, dim);
    
    for input_id = 1:input_num
        T = T + weights(input_id) * sqrtm(S_sqrt * cov_mat_cell{input_id} * S_sqrt);
    end

    if iter > 0
        cost_iter(iter) = trace(S) - 2 * trace(T);
    end
    
    % check the termination criterion
    if norm(S - T) < tolerance
        break;
    end

    S = (S_sqrt \ T) * (T / S_sqrt);

    % update the iteration counter
    iter = iter + 1;
end

WB_cov_mat = S;
cost_iter = cost_iter(1:iter) + min_cost;
min_cost = min_cost - trace(S);

end

