function dt = gen_dt(tvec) 
%% ==================================================================
%GEM_DT generate dt from t-vector
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-10-24
%
%   [INPUTS]
%   ---- Required ----
%   tvec                    vector of time
%
%   [OUTPUTS]
%   dt                      dt
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/10/24
%
%--------------------------------------------------------------------------

%%
    assert(isvector(tvec),'Input must be a vector.')
    assert(var(diff(tvec))<1e-3,'Non-uniform sampling detected.')

    dt = mean(diff(tvec));
    
end


