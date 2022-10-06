function pcg_flag_check( flg, itr, relres )
%PCG_FLAG_CHECK Displays appropriate error message based on pcg error flag

% flg = 0
% pcg converged to the desired tolerance tol within maxit iterations.

    flg_codes = {'pcg iterated maxit times but did not converge.',...
                 'Preconditioner M was ill-conditioned.',...
                 'pcg stagnated. (Two consecutive iterates were the same.)',...
                 'One of the scalar quantities calculated during pcg became too small or too large to continue computing.'};
    if flg > 0
        switch nargin
            case 2
                fprintf('PCG Flag (itr: %i):  ',itr);
            case 3
                fprintf('PCG Flag (itr: %i relres: %g):  ',itr, relres);
            otherwise
                %pass
        end
        fprintf('%s\n',flg_codes{flg});
    end
    
    if nargin >1
        fprintf('Number of PCG iterations: %i\n',itr);
    end

    return;

end

