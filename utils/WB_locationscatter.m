function [WB_mean_vec, WB_cov_mat, min_cost, cost_iter] = ...
    WB_locationscatter(mean_vec_cell, cov_mat_cell, weights, change_tol)
% Compute the 2-Wasserstein barycenter of probability measures from a
% location-scatter family via a fixed-point scheme
% Inputs: 
%       mean_vec_cell: cell array containing the mean vectors of the
%       input measures
%       cov_mat_cell: cell array containing the covariance matrices of the
%       input measures
%       weights: the weights defining the 2-Wasserstein barycenter
%       change_tol: when the change in the cost is less than change_tol, 
%       the algorithm terminates (default is 1e-6)
% Outputs: 
%       WB_mean_vec: the mean vector of the computed 2-Wasserstein
%       barycenter
%       cost: the minimized cost value

if ~exist('change_tol', 'var') || isempty(change_tol)
    change_tol = 1e-6;
end

input_num = length(mean_vec_cell);

assert(length(cov_mat_cell) == input_num ...
    || length(weights) == input_num, 'measures and weights mis-match');

dim = length(mean_vec_cell{1});

% the mean vector of the Wasserstein barycenter is simply the barycenter of
% the mean vectors of the input measures
WB_mean_vec = horzcat(mean_vec_cell{:}) * weights;

% compute the cost related to mean vectors and the trace of the covariance
% matrices
min_cost = 0;
for input_id = 1:input_num
    min_cost = min_cost + weights(input_id) ...
        * (sum((mean_vec_cell{input_id} - WB_mean_vec) .^ 2) ...
        + trace(cov_mat_cell{input_id})); 
end

S = eye(dim);
cost_var = inf;
iter = 0;

% usually the iteration terminates after a few iterations
cost_iter = zeros(100, 1);

while true
    S_sqrt = sqrtm(S);
    T = zeros(dim, dim);
    
    for input_id = 1:input_num
        T = T + weights(input_id) ...
            * sqrtm(S_sqrt * cov_mat_cell{input_id} * S_sqrt);
    end
    
    cost_var_new = trace(S) - 2 * trace(T);

    if iter > 0
        cost_iter(iter) = cost_var_new;
    end
    
    % check the termination criterion
    if abs(cost_var_new - cost_var) < change_tol
        break;
    end

    cost_var = cost_var_new;
    S = (S_sqrt \ T) * (T / S_sqrt);

    % update the iteration counter
    iter = iter + 1;
end

WB_cov_mat = S;
cost_iter = cost_iter(1:iter) + min_cost;
min_cost = min_cost + cost_var;

end

