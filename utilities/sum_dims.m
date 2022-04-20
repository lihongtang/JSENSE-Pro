function data_sum = sum_dims(data,dims) 
%% ==================================================================
%SUM_DIMS summation in selected dimensions
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-10-08
%
%   [INPUTS]
%   ---- Required ----
%   data                    data
%
%   ---- Optional ----
%   dims                    selected dimensions [1:ndims(data)]
%
%
%   [OUTPUTS]
%   data_sum                summed data
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/10/08
%       Modified by Yibo Zhao @ UIUC, 2018/11/20:
%           Sum over all dimensions when dims is not specified
%
%--------------------------------------------------------------------------

%% 
    if ~exist('dims','var') || isempty(dims)
        dims = 1:ndims(data);
    end

%%
    data_sum = data;
    for n = 1:length(dims)
        data_sum = sum(data_sum,dims(n));
    end
    
end
