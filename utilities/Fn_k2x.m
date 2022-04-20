function [ data2image ] = Fn_k2x( data2image,dims,dont_shift )
%% ==================================================================
%FN_K2X do ifft on selected dimensions
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-07-23
%
%   [INPUTS]
%   ---- Required ----
%   data2image              signal in freq/k domain
%   dims                    dimensions along which to do the fft
%
%   ---- Optional ----
%   dont_shift              if true, don't do any ifft shifts [false]
%
%   [OUTPUTS]
%   data2image              signal in time/spatial domain
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/07/23
%
%   Note:
%       F3_k2x(data) is equivalent to Fn_k2x(data,[1,2,3])
%       F_k2x(data)  is equivalent to Fn_k2x(data,[1,2])
%       F1_k2x(data) is equivalent to Fn_k2x(data,1)
%
%   See also F1_K2X, F_K2X, F3_K2X, F1_X2K, F_X2K, F3_X2K, FN_X2K
%--------------------------------------------------------------------------

%% --- parse the inputs ---
    if ~exist('dont_shift','var')||isempty(dont_shift)
        dont_shift = false;
    end
    
%% --- ifft ---
    
    scale_fctr = sqrt(numel_dims(data2image,dims));
    
    if dont_shift
        for ind_dim = 1:length(dims)
            data2image = ifft(data2image,[],dims(ind_dim));
        end
    else
        for ind_dim = 1:length(dims)
            data2image = ifftshift(ifft(fftshift(data2image,dims(ind_dim)),[],dims(ind_dim)),dims(ind_dim));
        end
    end
    
    data2image = data2image*scale_fctr;
    
    return;

end
