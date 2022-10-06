function [ image2data ] = Fn_x2k( image2data,dims,dont_shift,points)
%% ==================================================================
%FN_X2K do fft on selected dimensions
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-07-23
%
%   [INPUTS]
%   ---- Required ----
%   image2data              signal in time/spatial domain
%   dims                    dimensions along which to do the fft
%
%   ---- Optional ----
%   dont_shift              if true, don't do any fft shifts [false]
%
%   [OUTPUTS]
%   image2data              signal in freq/k domain
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/07/23
%
%   Note:
%       F3_x2k(image) is equivalent to Fn_x2k(image,[1,2,3])
%       F_x2k(image)  is equivalent to Fn_x2k(image,[1,2])
%       F1_x2k(image) is equivalent to Fn_x2k(image,1)
%
%   See also F1_X2K, F_X2K, F3_X2K, F1_K2X, F_K2X, F3_K2X, FN_K2X
%--------------------------------------------------------------------------

%% --- parse the inputs ---
    if ~exist('dont_shift','var')||isempty(dont_shift)
        dont_shift = false;
    end
    
%% --- fft ---
     if ~exist('points','var')||isempty(points)
        for ind_dim = 1:length(dims)
            points(ind_dim) = size(image2data,dims(ind_dim));
        end
    end
    scale_fctr = 1/sqrt(numel_dims(image2data,dims));
    
    if dont_shift
        for ind_dim = 1:length(dims)
            image2data = fft(image2data,points(ind_dim),dims(ind_dim));
        end
    else
        for ind_dim = 1:length(dims)
            image2data = ifftshift(fft(fftshift(image2data,dims(ind_dim)),points(ind_dim),dims(ind_dim)),dims(ind_dim));
        end
    end
    
    image2data = image2data*scale_fctr;
    
    return;
    
end
