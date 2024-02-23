classdef (Abstract) HasTractableQuadraticIntegrals < handle
    % Abstract class providing providing a property FirstMomentVec which is
    % equal to the expected value of x as well as a property 
    % SecondMomentMat which is equal to the expected value of x' * x
    
    properties(GetAccess = public, SetAccess = protected)
        % vector containing the first moment
        FirstMomentVec;

        % matrix containing the second moment
        SecondMomentMat;
    end

    methods(Access = public)
        function mean_vec = meanVector(obj)
            % Returns the mean vector of the probability measure
            % Output:
            %   mean_vec: the mean vector of the probability measure
            
            mean_vec = obj.FirstMomentVec;
        end

        function cov_mat = covarianceMatrix(obj)
            % Retruns the covariance matrix of the probability measure
            % Output:
            %   cov_mat: the covariance matrix of the probability measure

            cov_mat = obj.SecondMomentMat ...
                - obj.FirstMomentVec * obj.FirstMomentVec';
        end
    end
end

