function [new_coup_atom_indices, new_coup_probs] = ...
    discrete_reassembly_1D(coup_atoms, coup_probs, ...
    new_marg_atoms_cell, new_marg_probs_cell)
% Compute the reassembly of a discrete measure with discrete
% one-dimensional marginals
% Inputs: 
%   coup_atoms: matrix containing the atoms in the original coupling
%   coup_probs: vector containing the corresponding probabilities in the 
%   original coupling
%   new_marg_atoms_cell: cell array containing the atoms in each new 
%   marginal
%   new_marg_probs_cell: cell array containing the probabilities of atoms
%   in each new marginal
% Outputs: 
%   new_coup_atom_indices: matrix containing the atoms in the reassembly
%   represented by the indices of atoms in the new marginals
%   new_coup_probs: vector containing the corresponding probabilities in 
%   the reassembly

marg_num = length(new_marg_atoms_cell);

% check the inputs
assert(size(coup_atoms, 2) == marg_num, ...
    'atoms in the coupling mis-specified');
assert(size(coup_atoms, 1) == length(coup_probs), ...
    'probabilities in the coupling mis-specified');

% normalize the probabilities (in case there is a small numerical error)
coup_probs = coup_probs / sum(coup_probs);

for marg_id = 1:marg_num
    assert(length(new_marg_probs_cell{marg_id}) ...
        == length(new_marg_atoms_cell{marg_id}), ...
        'new marginals mis-specified');
end

% compute the reassembly

% remove atoms with zero probability
nonzero_list = coup_probs > 0;
new_coup_probs = coup_probs(nonzero_list);

% this matrix will contain two types of columns: those columns that have
% not been coupled will store the coordinate of atoms; those columns that
% have been coupled will store the knot indices in the new marginals
new_coup_atom_indices = coup_atoms(nonzero_list, :);

% begin the loop to compute a reassembly
for marg_id = 1:marg_num
    [~, sorted_marg_atom_indices] = sort(new_marg_atoms_cell{marg_id});
    sorted_marg_probs = ...
        new_marg_probs_cell{marg_id}(sorted_marg_atom_indices);
    
    % sort the atoms into ascending order in the i-th dimension
    [~, asc_order] = sort(new_coup_atom_indices(:, marg_id), 'ascend');
    sorted_atom_indices = new_coup_atom_indices(asc_order, :);
    sorted_probs = new_coup_probs(asc_order);
    
    % the atoms and probabilities in the glued distribution
    atom_max_num = size(sorted_atom_indices, 1) ...
        + length(sorted_atom_indices);
    glued_atoms_indices = zeros(atom_max_num, marg_num + 1);
    glued_probs = zeros(atom_max_num, 1);
    
    % initialize the counters
    joint_counter = 1;
    marg_counter = 1;
    glued_counter = 0;
    
    % begin the loop to couple the marginals
    while true
        % take out one atom
        prob_min = min(sorted_probs(joint_counter), ...
            sorted_marg_probs(marg_counter));
        
        % record the atom and its probability
        glued_counter = glued_counter + 1;
        glued_atoms_indices(glued_counter, :) = ...
            [sorted_atom_indices(joint_counter, :), ...
            sorted_marg_atom_indices(marg_counter)];
        glued_probs(glued_counter) = prob_min;
        
        % decrease the probability from the remaining probabilities
        sorted_probs(joint_counter) = sorted_probs(joint_counter) ...
            - prob_min;
        sorted_marg_probs(marg_counter) = ...
            sorted_marg_probs(marg_counter) - prob_min;
        
        % advance the counters
        if sorted_probs(joint_counter) <= 1e-14
            joint_counter = joint_counter + 1;
        end
        
        if sorted_marg_probs(marg_counter) <= 1e-14
            marg_counter = marg_counter + 1;
        end
        
        if joint_counter > size(sorted_atom_indices, 1) ...
                || marg_counter > length(sorted_marg_atom_indices)
            break;
        end
    end
    
    % update the distribution
    glued_atoms_indices = glued_atoms_indices(1:glued_counter, :);
    new_coup_atom_indices = glued_atoms_indices(:, 1:marg_num);
    new_coup_atom_indices(:, marg_id) = glued_atoms_indices(:, end);
    new_coup_probs = glued_probs(1:glued_counter);
end

[new_coup_atom_indices, ~, umap] = unique(new_coup_atom_indices, 'rows', 'stable');
new_coup_probs = accumarray(umap, new_coup_probs, [size(new_coup_atom_indices, 1), 1]);

end

