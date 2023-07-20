function [rho, hyp_full, lin_full] ...
    = hyperbolic_Dirichlet_tessellation_polar(theta, ...
    hyp_a, hyp_c, hyp_intercept, lin_c, lin_intercept, full_info)
% Express a region in hyperbolic Dirichlet tessellation in the polar
% coordinate where the origin is a center. The region is bound by
% hyperbolas and straight lines. The function returns the distance from the
% origin for the given inputs theta in radian. The parameters have been
% preprocessed in the sense that after rotation by hyp_intercept, a
% hyperbola will be in its canonical orientation, and that after rotation
% by lin_intercept, a straight line will be parallel to the x-axis above
% the origin
% Inputs: 
%   theta: the vector-valued input in radian
%   hyp_a: signed a-values of the hyperbolas (a positive value indicates
%   the left branch and a negative value indicates the right branch)
%   hyp_c: c-values of the hyperbolas (i.e., distances from other centers
%   divided by 2)
%   hyp_intercept: rotation angles (in radian) to set the hyperbolas into
%   their canonical orientation
%   lin_c: distances of the straight lines from the origin
%   lin_intercept: rotation angles (in radian) to set the straight lines to
%   be parallel to the x-axis and above the origin
%   full_info: boolean indicating whether to return the full distance
%   information in matrices (default is false)
% Outputs:
%   rho: the distances from the origin along the boundary of the region
%   hyp_full: matrix containing full distance information about boundaries
%   formed by hyperbolas (equal to [] if full_info = false)
%   lin_full: matrix containing full distance information about boundaires
%   formed by hyperplanes (equal to [] if full_info = false)

if ~exist('full_info', 'var') || isempty(full_info)
    full_info = false;
end

% compute b^2
hyp_b_sq = hyp_c .^ 2 - hyp_a .^ 2;

if any(hyp_b_sq <= 0 & hyp_a > 0)
    % if c < a, then the set is empty; if c = a, then the set is a line,
    % which is negligible
    rho = zeros(length(theta), 1);

    hyp_full = [];
    lin_full = [];
else
    list = hyp_b_sq > 0;

    hyp_rho = hyp_b_sq(list)' ./ (hyp_a(list)' - hyp_c(list)' ...
        .* cos(theta - hyp_intercept(list)'));
    hyp_rho(hyp_rho < 0) = inf;

    lin_sin = sin(theta - lin_intercept');
    lin_rho = lin_c' ./ lin_sin;
    lin_rho(lin_sin < 0) = inf;

    rho = min([hyp_rho, lin_rho], [], 2);

    if full_info
        hyp_full = inf(length(theta), length(list));
        hyp_full(:, list) = hyp_rho;

        lin_full = lin_rho;
    else
        hyp_full = [];
        lin_full = [];
    end
end

end

