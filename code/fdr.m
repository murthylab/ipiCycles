function qMatrix = fdr(pMatrix, minDiff)
% qMatrix = fdr(pMatrix, minDiff)
assert(minDiff>0, 'minDiff needs to be >0, is %1.2f', minDiff)
pVector = pMatrix(triu(pMatrix,minDiff)>0);  % vectorize upper triangular values with min diff
% [~, qVector] = mafdr(pVector, 'BHFDR', false);     % do fdr
[h, crit_p, adj_ci_cvrg, qVector]=fdr_bh(pVector,0.05,'pdep');
qMatrix = triu(ones(size(pMatrix)),minDiff); % prepare new matrix for qVec
qMatrix(qMatrix==1) = qVector; % put new q-values in matrix
