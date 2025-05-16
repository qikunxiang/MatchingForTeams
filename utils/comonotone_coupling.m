function [coup_atom_indices, coup_probs] = comonotone_coupling( ...
    marg_atoms_cell, marg_probs_cell)
% Compute the comonotone coupling of finitely-supported marginal
% distributions
% Inputs:
%   marg_atoms_cell: cell array containing the atoms of the marginal
%   distributions
%   marg_probs_cell: cell array containing the probabilities of the
%   marginal distributions
% Outputs:
%   coup_atom_indices: matrix containing atoms in the comonotone coupling
%   described by the indices of the atoms from the marginals
%   coup_probs: the probabilities corresponding to the atoms in the
%   coupling

marg_num = length(marg_atoms_cell);
assert(length(marg_probs_cell), 'marginals mis-specified');

atom_num_list = zeros(marg_num, 1);

for marg_id = 1:marg_num
    % check the inputs
    assert(length(marg_atoms_cell{marg_id}) ...
        == length(marg_probs_cell{marg_id}), ...
        'marginals mis-specified');

    % normalize the probabilities (in case there is a small numerical
    % error)
    marg_probs_cell{marg_id} = marg_probs_cell{marg_id} ...
        / sum(marg_probs_cell{marg_id});

    nonzero_list = marg_probs_cell{marg_id} > 0;
    marg_probs_cell{marg_id} = marg_probs_cell{marg_id}(nonzero_list);
    marg_atoms_cell{marg_id} = marg_atoms_cell{marg_id}(nonzero_list);

    atom_num_list(marg_id) = length(marg_probs_cell{marg_id});
end

atom_max_num = max(atom_num_list);

% place atoms and probabilities in columns of matrices
atom_indices_mat = zeros(atom_max_num, marg_num);
probs_mat = zeros(atom_max_num, marg_num);

for marg_id = 1:marg_num
    [~, sorted_order] = sort(marg_atoms_cell{marg_id}, 'ascend');
    atom_indices_mat(1:atom_num_list(marg_id), marg_id) = sorted_order;
    probs_mat(1:atom_num_list(marg_id), marg_id) = ...
        marg_probs_cell{marg_id}(sorted_order);
end

% initialization before the loop
joint_atom_indices = zeros(sum(atom_num_list), marg_num);
joint_probs = zeros(sum(atom_num_list), 1);

% the counter for the number of atoms in the coupled measure
counter = 0;

% the index of the first remaining atom in each marginal
index_marg = ones(marg_num, 1);

% loop until all atoms have been coupled
while all(index_marg <= atom_num_list)
    % convert the row indices to linear indices in the
    % probability matrix
    lin_ind = sub2ind([atom_max_num, marg_num], index_marg, (1:marg_num)');

    % select the atom with the smallest probability
    prob_min = min(probs_mat(lin_ind));

    % couple the atoms
    counter = counter + 1;
    joint_atom_indices(counter, :) = atom_indices_mat(lin_ind);
    joint_probs(counter) = prob_min;

    % remove the probabilities that have already been coupled
    probs_mat(lin_ind) = probs_mat(lin_ind) - prob_min;

    % remove the atoms that have been fully coupled
    remove_list = probs_mat(lin_ind) <= 1e-14;
    index_marg(remove_list) = index_marg(remove_list) + 1;
end

% prepare output
coup_atom_indices = joint_atom_indices(1:counter, :);
coup_probs = joint_probs(1:counter);

end