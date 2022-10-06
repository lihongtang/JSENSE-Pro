function SX = shrink(X,maxval)
% SHRINK Shrinkage operator (see reference paper below)
%
% This code is based on the algorithm given in the following paper:
%   Candes, Xiaodong, Ma, Wright, "Robust Principal component Analysis".
%   2011. ACM, Journal of.
%
% Author: Bryan A. Clifford
% Date Created: 05/22/2013
%--------------------------------------------------------------------------

    SX = abs(X) - maxval;
    SX = sign(X).*SX.*(SX > 0);

    return

end
