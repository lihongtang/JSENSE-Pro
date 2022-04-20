function [U, S, V] = svd_econ(C)
% SVD_ECON Performs svd(C,'econ') but checks for bad convergence and 
% then tries again using svds
%--------------------------------------------------------------------------

% Estimate subspace
try
    [U, S, V] = svd(C,'econ');
    if strcmp(lastwarn,'SVD did not converge at index = 1.')
        fprintf('\nCaught convergence problem in svd trying again using svds.\n');
        lastwarn('');
        [U, S, V] = svds(C,min(size(C)));
    end
catch err
    fprintf('\nCaught convergence problem in svd trying again using svds.\n');
        lastwarn('');
    [U, S, V] = svds(C,min(size(C)));
end
return;