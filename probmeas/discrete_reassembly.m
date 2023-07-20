function [reassembly_atoms, reassembly_probs] = discrete_reassembly( ...
    coup_atoms, coup_probs, marg_coup_atoms_cell, marg_coup_probs_cell, ...
    reassembled_marg_list)
% Compute a reassembly of a discrete measure with discrete marginals
% Inputs: 
%       coup_atoms: the atoms in the original discrete measure; each entry
%       represents the index of an atom in the respective marginal
%       coup_probs: the corresponding probabilities in the original
%       discrete measure
%       marg_coup_atoms_cell: cell array containing the atoms in the
%       coupling between the original marginal and the new marginal; each
%       entry represents the index of an atom
%       marg_coup_probs_cell: cell array containing the probabilities of
%       atoms in the coupling between the original marginal and the new
%       marginal
%       reassembled_marg_list: list of marginal indices to reassemble, the
%       rest will be skipped (default is 1:marg_num)
% Outputs: 
%       new_atoms: the atoms in the computed reassembly
%       new_probs: the probabilities of atoms in the computed reassembly

marg_num = size(coup_atoms, 2);

if ~exist('reassembled_marg_list', 'var') || isempty(reassembled_marg_list)
    reassembled_marg_list = 1:marg_num;
end

% check the inputs
assert(length(marg_coup_probs_cell) == length(reassembled_marg_list), ...
    'marginals mis-specified');
assert(size(coup_atoms, 1) == length(coup_probs), ...
    'original discrete measure mis-specified');
assert(min(reassembled_marg_list) >= 1 ...
    && max(reassembled_marg_list) <= marg_num, ...
    'reassembled marginal list mis-specified');

% remove atoms with zero probability
nonzero_list = coup_probs > 0;
coup_probs = coup_probs(nonzero_list);
coup_atoms = coup_atoms(nonzero_list, :);

% initialize the reassembly measure
reassembly_atoms = coup_atoms;
reassembly_probs = coup_probs;

input_counter = 1;

% begin the loop to compute a reassembly, where a single marginal is
% updated in each iteration
for marg_id = reassembled_marg_list
    % sort the atoms in the joint measure into descending order
    [reassembly_probs, sort_ind] = sort(reassembly_probs, 'descend');
    reassembly_atoms = reassembly_atoms(sort_ind, :);
    
    % get the i-th coupling of marginals
    marg_coup_atoms_i = marg_coup_atoms_cell{input_counter};
    marg_coup_probs_i = marg_coup_probs_cell{input_counter};
    input_counter = input_counter + 1;
    
    % remove atoms with zero probability from the coupling
    nonzero_list = marg_coup_probs_i > 0;
    marg_coup_probs_i = marg_coup_probs_i(nonzero_list);
    marg_coup_atoms_i = marg_coup_atoms_i(nonzero_list, :);
    
    % sort the atoms in the coupling in descending order; this will 
    % guarantee that the first atom found later will be the one with the 
    % largest probability
    [marg_coup_probs_i, sort_ind] = sort(marg_coup_probs_i, 'descend');
    marg_coup_atoms_i = marg_coup_atoms_i(sort_ind, :);
    
    % the atoms and probabilities in the glued measure
    atom_num_max = length(marg_coup_probs_i) + length(reassembly_probs);
    glued_atoms = zeros(atom_num_max, marg_num + 1);
    glued_probs = zeros(atom_num_max, 1);
    
    % initialize the counters
    counter_joint = 1;
    counter_glued = 0;
    
    % begin the loop to couple the marginals
    while counter_joint <= size(reassembly_probs, 1)
        % take out one atom from the coupling of the marginals
        marg_coup_ind = find(marg_coup_atoms_i(:, 1) ...
            == reassembly_atoms(counter_joint, marg_id), 1, 'first');

        if isempty(marg_coup_ind)
            % if this atom is not present in the coupling of the marginals,
            % the probability of this atom must be very small; this can
            % happen when the atom is missing due to numerical inaccuracies
            missing_atom = reassembly_atoms(counter_joint, marg_id);
            missing_prob = sum(reassembly_probs( ...
                reassembly_atoms(:, marg_id) == missing_atom));

            assert(missing_prob < 1e-8, ...
                'an atom is missing in the marginal coupling');

            % add the missing probability to an artificial atom coupled
            % with the first atom in the new marginal
            marg_coup_atoms_i = [marg_coup_atoms_i; missing_atom, 1]; ...
                %#ok<AGROW> 
            marg_coup_probs_i = [marg_coup_probs_i; missing_prob]; ...
                %#ok<AGROW> 
            
            marg_coup_ind = length(marg_coup_probs_i);
            prob_min = reassembly_probs(counter_joint);
        else
            if marg_coup_probs_i(marg_coup_ind) <= 0
                % in the case where the marginal probability of an atom in
                % the joint distribution is slightly larger than the
                % probability in the coupling of marginals (due to
                % numerical inaccuracies), skip this atom
                counter_joint = counter_joint + 1;
            
                [marg_coup_probs_i, sort_ind] = ...
                    sort(marg_coup_probs_i, 'descend');
                marg_coup_atoms_i = marg_coup_atoms_i(sort_ind, :);

                continue;
            end

            % take the smaller of the two probabilities
            prob_min = min(reassembly_probs(counter_joint), ...
                marg_coup_probs_i(marg_coup_ind));
        end
        
        % record the atom and its probability in the glued measure
        counter_glued = counter_glued + 1;
        glued_atoms(counter_glued, :) = ...
            [reassembly_atoms(counter_joint, :), ...
            marg_coup_atoms_i(marg_coup_ind, 2)];
        glued_probs(counter_glued) = prob_min;
        
        % decrease the probability from the remaining probabilities
        reassembly_probs(counter_joint) = ...
            reassembly_probs(counter_joint) - prob_min;
        marg_coup_probs_i(marg_coup_ind) = ...
            marg_coup_probs_i(marg_coup_ind) - prob_min;
        
        % advance the counter
        if reassembly_probs(counter_joint) <= 0
            counter_joint = counter_joint + 1;
        end
        
        % sort the atoms in the coupling again
        [marg_coup_probs_i, sort_ind] = sort(marg_coup_probs_i, 'descend');
        marg_coup_atoms_i = marg_coup_atoms_i(sort_ind, :);
    end
    
    % update the reassembly
    glued_atoms = glued_atoms(1:counter_glued, :);
    reassembly_atoms = glued_atoms(:, [1:marg_id - 1, marg_num + 1, ...
        marg_id + 1:marg_num]);
    reassembly_probs = glued_probs(1:counter_glued);
end

% normalize the updated coupling to remove potential numerical inaccuracies
reassembly_probs = reassembly_probs / sum(reassembly_probs);

end

