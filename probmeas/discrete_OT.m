function [coup_atoms, coup_probs, OT_cost] ...
    = discrete_OT(probs1, probs2, cost_mat)
% Compute a discrete optimal transport between two discrete measures by
% solving a linear programming problem
% Inputs:
%       probs1: the probabilities of atoms in the first discrete measure
%       probs2: the probabilities of atoms in the second discrete measure
%       cost_mat: cost matrix where each row corresponds to an atom in the
%       first measure and each column corresponds to an atom in the second
%       measure
% Outputs: 
%       coup_atoms: two-column matrix containing the coupling of the atoms
%       from the two discrete measures
%       coup_probs: the probabilities of atoms in the coupled measure
%       OT_cost: the computed optimal transport cost

n1 = length(probs1);
n2 = length(probs2);

[I2, I1] = meshgrid(1:n2, 1:n1);

C1_i = repelem((1:n1)', n2, 1);
C1_j = reshape(reshape(1:n1 * n2, n1, n2)', [], 1);
C1 = sparse(C1_i, C1_j, 1, n1, n1 * n2);

C2_i = repelem((1:n2)', n1, 1);
C2_j = (1:n1 * n2)';
C2 = sparse(C2_i, C2_j, 1, n2, n1 * n2);


model = struct;
model.obj = cost_mat(:);
model.A = [C1; C2];
model.rhs = [probs1; probs2];
model.sense = '=';
model.lb = zeros(n1 * n2, 1);

params = struct;
params.OutputFlag = 0;
params.FeasibilityTol = 1e-9;

result = gurobi(model, params);

OT_cost = result.objval;

supp = find(result.x > 0);
coup_atoms = [I1(supp), I2(supp)];
coup_probs = result.x(supp);

end

