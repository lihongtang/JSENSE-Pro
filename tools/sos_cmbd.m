function sxt = sos_cmbd(sxtc,dims)
%% ==================================================================
%SOS_CMBD combine individual coil data by sum-of-squares, or compute
%spectral l2 integral
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-05-10
%
%   [INPUTS]
%   ---- Required ----
%   sxtc                    individual-coil EPSI data  [y,x,z,t,coil,(frame,ave...)]
%
%   ---- Optional ----
%   dim                     the dimension of coils
%                           [ default: min(5,ndims(sxtc)) ]
%
%   [OUTPUTS]
%   sxt                     coil-combined EPSI data    [y,x,z,t,1,(frame,ave...)]
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/05/10
%       Modified by Yibo Zhao @ UIUC, 2018/11/04:
%           Acccept a vector of dimensions
%
%--------------------------------------------------------------------------

    if(nargin <2)
        if(ndims(sxtc)<4)
            disp(['Data have only ',num2str(ndims(sxtc)),' dimensions; combine along the *last* dimension.']);
            dims = ndims(sxtc);
        elseif(ndims(sxtc)==4)
            dims = 4;
        elseif(ndims(sxtc)>5)
            dims = 5;
        else
            dims = 5;
        end
    end
    
    sxt = sxtc;
    for n = 1:length(dims)
        sxt = sqrt(sum(abs(sxt).^2,dims(n)));
    end
    
end
    
    
