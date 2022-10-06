function data_flipped = flip_dims(data,dims) 
%% ==================================================================
%FLIP_DIMS flip in selected dimensions
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-09-04
%
%   [INPUTS]
%   ---- Required ----
%   data                    data
%   dims                    selected dimensions
%
%
%   [OUTPUTS]
%   data_flipped            flipped data
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/09/04
%
%--------------------------------------------------------------------------

%% 
    data_flipped = data;
    for n = 1:length(dims)
        data_flipped = flip(data_flipped,dims(n));
    end
    
end
