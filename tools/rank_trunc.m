function DX = rank_trunc(X,maxval,opt)
% RANK_TRUNC Singular value thresholding operator (see reference paper below)
%
% This code is based on the algorithm given in the following paper:
%   Candes, Xiaodong, Ma, Wright, "Robust Principal component Analysis".
%   2011. ACM, Journal of.
%
% Author: Bryan A. Clifford
% Date Created: 05/22/2013
%--------------------------------------------------------------------------

    [U S V] = svd_econ(X);
    DX = U*shrink(S,maxval)*V';
    
    if opt.watch
        s = diag(S);
        maxidx = find(s>maxval,1,'last');
        
        semilogy(opt.ax,diag(S),'color','b');
        title({'\bfNormalized Singular Values\rm';...
               'of U - Z/b'});
        xlabel('singular value index');
        ylabel('singular value');
        
        if isempty(maxidx);
            fprintf('*** Bad threshold value! (maxval = %e) ***\n',maxval);
        else
            line([maxidx maxidx],get(opt.ax,'ylim'),'color','k');
        end
        
        drawnow;
    end

    return

end