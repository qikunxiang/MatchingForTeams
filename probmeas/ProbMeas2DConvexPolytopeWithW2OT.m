classdef (Abstract) ProbMeas2DConvexPolytopeWithW2OT ...
        < ProbMeas2DConvexPolytope
    % Abstract class for 2-dimensional probability measures supported within a convex polytope which allows for the computation of 
    % optimal transport with respect to the squared Euclidean distance, i.e., Wasserstein-2 optimal transport. The abstract methods can
    % be implemented for continuous probability measures. 

    properties(SetAccess = protected, GetAccess = public)
        % struct storing information about the Wasserstein-2 optimal transference plan and the way to sample from the conditional
        % measure 
        W2OT = struct;
    end
    
    methods(Access = public)
        
        function [OT_weights, ...
                OT_cost, ...
                OT_output] ...
                = computeW2OptimalTransport(obj, ...
                atoms, ...
                probs, ...
                init_weights, ...
                optim_options)
            % Compute the semi-discrete optimal transport with respect to
            % the squared Euclidean distance
            % Inputs:
            %   atoms: two-column matrix containing the locations of the atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in the discrete measure
            %   init_weights: the initial weights in the optimization (default is the all-zero vector)
            %   optim_options: struct containing options used in the optimization
            % Outputs:
            %   OT_weights: the optimized weights in the Laguerre diagram
            %   OT_cost: the optimal transport cost
            %   OT_output: struct containing additional information with fields
            %       exitflag: the exit flag of the fminunc function
            %       optim_output: struct output by the fminunc function
            %       final_grad: the gradient at termination

            atom_num = size(atoms, 1);

            if ~exist('init_weights', 'var') || isempty(init_weights)
                init_weights = zeros(atom_num, 1);
            end

            obj.prepareW2OptimalTransport(atoms, probs);

            options = optimoptions('fminunc', ...
                'Display', 'none', ...
                'Algorithm', 'quasi-newton', ...
                'SpecifyObjectiveGradient', true, ...
                'MaxFunctionEvaluations', 1e5, ...
                'MaxIterations', 1e5, ...
                'FunctionTolerance', 0, ...
                'StepTolerance', 0, ...
                'OptimalityTolerance', 1e-5);

            % set the additional optimization options
            if exist('optim_options', 'var') && ~isempty(optim_options)
                optim_fields = fieldnames(optim_options);
                optim_values = struct2cell(optim_options);

                for f_id = 1:length(optim_fields)
                    options.(optim_fields{f_id}) = optim_values{f_id};
                end
            end

            OT_obj_func = @(w)(obj.evaluateW2OTMinObjective(w));

            OT_output = struct;
            [OT_weights, ...
                objective_min, ...
                OT_output.exitflag, ...
                OT_output.optim_output, ...
                OT_output.final_grad] ...
                = fminunc(OT_obj_func, init_weights, options);
            OT_cost = -objective_min;

            if OT_output.exitflag == 1
                obj.W2OT.Coupled = true;
                obj.W2OT.Weights = OT_weights;
                obj.W2OT.Cost = OT_cost;
            else
                error('Error encountered during the computation of Wasserstein-2 optimal transport');
            end

            obj.postProcessW2OptimalTransport();
        end

        function info = saveW2OptimalTransportInfo(obj)
            % Save optimal transport related information in a struct that can be used later
            % Output:
            %   info: struct that copies the fields from obj.W2OT

            info = struct;
            info.atoms = obj.W2OT.DiscMeas.Atoms;
            info.probs = obj.W2OT.DiscMeas.Probs;
            info.bounding_box = struct;
            info.bounding_box.bottom_left = obj.W2OT.BoundingBox.BottomLeft;
            info.bounding_box.top_right = obj.W2OT.BoundingBox.TopRight;
            info.weights = obj.W2OT.Weights;
            info.cost = obj.W2OT.Cost;
        end

        function loadW2OptimalTransportInfo(obj, info)
            % Load optimal transport related information from a previously saved struct
            % Input:
            %   info: struct that will be copied to the fields of obj.W2OT

            obj.W2OT.DiscMeas = struct('Atoms', info.atoms, 'Probs', info.probs);
            obj.W2OT.BoundingBox = struct;
            obj.W2OT.BoundingBox.BottomLeft = info.bounding_box.bottom_left;
            obj.W2OT.BoundingBox.TopRight = info.bounding_box.top_right;
            obj.W2OT.Weights = info.weights;
            obj.W2OT.Cost = info.cost;
            obj.W2OT.Coupled = true;
        end

        function testW2OTObjective(obj, ...
                atoms, ...
                probs, ...
                OT_weights, ...
                grid_pts_x, ...
                grid_pts_y)
            % Test the objective and gradient in the computation of semi-discrete Wasserstein-2 optimal transport by comparing with 
            % two-dimensional integration over a grid
            % Inputs: 
            %   atoms: two-column matrix containing the locations of the atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in the discrete measure
            %   OT_weights: weight vector in the objective (including the first component)
            %   grid_pts_x: grid points on the x-axis for the integration
            %   grid_pts_y: grid points on the y-axis for the integration

            atom_num = size(atoms, 1);

            obj.prepareW2OptimalTransport(atoms, probs);

            % the objective value and gradient computed by the integration in the polar coordinates
            [val, grad] = obj.evaluateW2OTMinObjective(OT_weights);

            [grid_x, grid_y] = meshgrid(grid_pts_x, grid_pts_y);

            % approximate the probability measure with a discrete measure supported on the grid
            grid_density = obj.densityFunction([grid_x(:), grid_y(:)]);
            grid_density = grid_density / sum(grid_density);

            weighted_sqdist = (grid_x(:) - atoms(:, 1)') .^ 2 + (grid_y(:) - atoms(:, 2)') .^ 2 - OT_weights';
            [min_weighted_sqdist, min_id] = min(weighted_sqdist, [], 2);

            val_grid = -(OT_weights' * probs) - min_weighted_sqdist' * grid_density;
            grad_grid = accumarray(min_id, grid_density, [atom_num, 1]) - probs;
            grad_grid = grad_grid(2:end);

            fprintf('objective via integration: %.6f\n', val);
            fprintf('objective via grid: %.6f\n', val_grid);
            fprintf('difference: %e\n', abs(val - val_grid));

            fprintf('gradients:\n');
            display([grad, grad_grid]);

            fprintf('max diff. in gradients: %e\n', max(abs(grad - grad_grid)));
        end

        function testW2OTGradient(obj, atoms, probs, OT_weights, nudge)
            % Check the gradient of the objective function in the computation of semi-discrete Wasserstein-2 optimal transport
            % Inputs: 
            %   atoms: two-column matrix containing the locations of the atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in the discrete measure
            %   OT_weights: weight vector in the objective (including the first component)
            %   nudge: the size of the disturbance in finite differencing (default is 1e-6)

            if ~exist('nudge', 'var') || isempty(nudge)
                nudge = 1e-6;
            end

            atom_num = size(atoms, 1);

            obj.prepareW2OptimalTransport(atoms, probs);

            % the objective value and gradient computed by the integration in the polar coordinates
            [val, grad] = obj.evaluateW2OTMinObjective(OT_weights);

            grad_fd = zeros(atom_num - 1, 1);

            for comp_id = 1:atom_num
                OT_weights_nudged = OT_weights;
                OT_weights_nudged(comp_id) = OT_weights_nudged(comp_id) + nudge;

                val_nudged = obj.evaluateW2OTMinObjective(OT_weights_nudged);

                grad_fd(comp_id - 1) = (val_nudged - val) / nudge;
            end

            fprintf('gradients:\n');
            display([grad, grad_fd]);

            fprintf('max diff. in gradients: %e\n', max(abs(grad - grad_fd)));
        end

        function samp_cell = conditionalRandSampleFromW2Coupling(obj, atom_list, samp_num_list, rand_stream)
            % Randomly generate samples from the conditional distributions given the Wasserstein-2 optimally coupled atoms in the
            % discrete measure 
            % Inputs: 
            %   atom_list: vector containing the atom indices to be sampling from
            %   samp_num_list: vector containing number of samples coupled with each of the atoms in atom_list
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream) 
            % Output:
            %   samp_cell: cell array where each cell contains a two-column array containing the generated samples

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~isfield(obj.W2OT, 'Coupled') || ~obj.W2OT.Coupled
                error('must set the coupled discrete measure first');
            end

            cell_num = length(atom_list);
            
            assert(length(samp_num_list) == cell_num, 'number of atoms mismatch');

            samp_cell = cell(atom_num, 1);

            for cell_id = 1:cell_num
                % call the abstract function to perform the rejection
                % sampling
                samp_cell{cell_id} = obj.condRejSampFromW2Coupling(atom_list(cell_id), samp_num_list(cell_id), rand_stream);
            end
        end

        function [samps, disc_samps] = randSampleW2Coupled(obj, ...
                samp_num, rand_stream)
            % Randomly generate samples from the Wasserstein-2 optimal coupling of the probability measure and the discrete measure
            % Inputs: 
            %   samp_num: number of samples to generate
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream) 
            % Output:
            %   samps: two-column matrix containing the generated samples
            %   disc_samps: vector containing the corresponding sampled atom indices in the discrete measure

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            if ~isfield(obj.W2OT, 'Coupled') || ~obj.W2OT.Coupled
                error('must set the Wasserstein-2 coupled discrete measure first');
            end

            atom_num = length(obj.W2OT.DiscMeas.Probs);

            % generate samples from the discrete measure
            disc_samps = randsample(rand_stream, atom_num, samp_num, true, obj.W2OT.DiscMeas.Probs);

            % count the number of samples for each cell
            samp_num_list = accumarray(disc_samps, ones(samp_num, 1), [atom_num, 1]);

            samp_cell = obj.conditionalRandSampleFromW2Coupling((1:atom_num)', samp_num_list, rand_stream);
            
            samps = vertcat(samp_cell{:});

            % order the samples according the discrete samples
            [~, ordering] = sort(disc_samps, 'ascend');
            samps(ordering, :) = samps;
        end

        function [power_cell_list, ...
                power_cell_triangle_list, ...
                power_cell_indices, ...
                mesh_indices] = ...
                computePowerDiagram(obj, anchors, weights)
            % Compute the contour of each cell in the power diagram with given anchors and weights
            % Inputs: 
            %   anchors: two-column matrix containing the locations of the anchor points
            %   weights: weights in the power diagram
            % Outputs:
            %   power_cell_list: vector containing polyshape objects describing the power cells
            %   power_cell_triangle_list: vector containing polyshape objects describing triangles in the triangulation of the
            %   power cell intersection with the density mesh
            %   power_cell_indices: vector containing the indices of power cells that correspond to the triangles in the power cell
            %   triangulation
            %   mesh_indices: vector containing the indices of mesh triangles that correspond to the the triangles in the power cell 
            %   triangulation

            anchor_num = size(anchors, 1);

            bbox_bottomleft = min(min(obj.Supp.TriangleVertices, [], 1), min(anchors, [], 1))';
            bbox_topright = max(max(obj.Supp.TriangleVertices, [], 1), max(anchors, [], 1))';
            
            % shift the weights to make them all positive so that they can  be interpreted as squared radii of circles in a power 
            % diagram
            squared_radii = weights - min(weights) + 1;

            % invoke the mex function to compute a triangular mesh formed by the intersections of power cells and 
            [output_vertices, ...
                output_triangulation, ...
                power_cell_indices, ...
                mesh_indices] ...
                = mesh_intersect_power_diagram( ...
                obj.Supp.TriangleVertices, ...
                obj.Supp.Triangles, ...
                anchors, ...
                squared_radii, ...
                bbox_bottomleft, ...
                bbox_topright);

            triangle_num = size(output_triangulation, 1);
            triangle_cell = cell(triangle_num, 1);

            for triangle_id = 1:triangle_num
                cell_vertices = output_vertices(output_triangulation(triangle_id, :), :);
                triangle_cell{triangle_id} = polyshape(cell_vertices);
            end

            power_cell_triangle_list = [triangle_cell{:}];

            poly_cell = cell(anchor_num, 1);

            for anchor_id = 1:anchor_num
                poly_cell{anchor_id} = union(power_cell_triangle_list(power_cell_indices == anchor_id), ...
                    'Simplify', true, 'KeepCollinearPoints', false);
            end

            power_cell_list = [poly_cell{:}];
        end
    end

    methods(Access = protected)

        function prepareW2OptimalTransport(obj, atoms, probs)
            % Prepare quantities for the computation of semi-discrete Wasserstein-2 optimal transport.
            % Inputs:
            %   atoms: two-column matrix containing the locations of the atoms in the discrete measure
            %   probs: vector containing the probabilities of the atoms in the discrete measure

            % store the discrete measure 
            obj.W2OT.DiscMeas = struct('Atoms', atoms, 'Probs', probs);

            obj.W2OT.BoundingBox = struct( ...
                'BottomLeft', min(min(atoms, [], 1), min(obj.Supp.TriangleVertices, [], 1))', ...
                'TopRight', max(max(atoms, [], 1), max(obj.Supp.TriangleVertices, [], 1))');
        end

        function [val, grad] = evaluateW2OTMinObjective(obj, OT_weights)
            % The objective function of the minimization problem when solving semi-discrete optimal transport (OT) problem.
            % Inputs:
            %   OT_weights: the input of the objective function corresponding to the dual weights in the power diagram
            % Outputs:
            %   val: the computed objective value
            %   grad: the computed gradient of the objective function

            % shift the weights to make them all positive so that they can be interpreted as squared radii of circles in a power 
            % diagram
            squared_radii = OT_weights - min(OT_weights) + 1;

            % invoke the mex function to compute a triangular mesh formed by the intersections of power cells and 
            [output_vertices, ...
                output_triangulation, ...
                power_cell_indices, ...
                mesh_indices] ...
                = mesh_intersect_power_diagram( ...
                obj.Supp.TriangleVertices, ...
                obj.Supp.Triangles, ...
                obj.W2OT.DiscMeas.Atoms, ...
                squared_radii, ...
                obj.W2OT.BoundingBox.BottomLeft, ...
                obj.W2OT.BoundingBox.TopRight);

            % compute a matrix containing the centers of the power cells as rows
            power_cell_centers = obj.W2OT.DiscMeas.Atoms(power_cell_indices, :);

            [int1, int2] = obj.computeW2InnerIntegral(output_vertices, output_triangulation, power_cell_centers, mesh_indices);

            % regroup the integral values based on the power cell indices; this evaluates the probability within each power cell
            int1 = accumarray(power_cell_indices, int1, [length(OT_weights), 1]);

            val = -(OT_weights' * obj.W2OT.DiscMeas.Probs) - sum(int2) + OT_weights' * int1;
            grad = int1 - obj.W2OT.DiscMeas.Probs;
        end
    end

    methods(Abstract, Access = protected)
        % Compute integrals of constant 1 functions as well as the squared distance functions (from given centers) within triangles 
        % with respect to the probability measure
        [int1, int2] = computeW2InnerIntegral(obj, vertices, triangles, centers, mesh_indices);

        % Post-processing after computing semi-discrete Wasserstein-2 optimal transport including computing the contour of each cell in 
        % the Laguerre diagram and preparing for conditional rejection sampling 
        postProcessW2OptimalTransport(obj);

        % Perform rejection sampling of the conditional distribution
        samps = condRejSampFromW2Coupling(obj, atom_id, samp_num, rand_stream);
    end
end