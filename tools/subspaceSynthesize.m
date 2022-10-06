function [sxt_syn] = subspaceSynthesize(cxr,V_subspace,dim_r_data,dim_t_V)
%% ======================================================================
%SUBSPACESYNTHESIZE synthesis operation of a basis for a data: d = V*c
% =======================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2019-05-31
%
%   [INPUTS]
%   ---- Required ----
%   cxr                     spatial coefficients
%   V_subspace              temporal basis
%
%   ---- Optional ----
%   dim_r_data              the dimension of rank in data [4]
%   dim_t_basis             the dimension of t in basis [1]
%
%   [OUTPUTS]
%   sxt_syn                synthetic (x-t) domain data
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2019/05/31
%       Modified by Yibo Zhao @ UIUC, 2019/10/30:
%           Take care of rank-1 case.
%
%   See also SUBSPACEANALYZE, PROJONTOSUBSPACE
%
%--------------------------------------------------------------------------

%% parse inputs
    if ~exist('dim_r_data','var')||isempty(dim_r_data)
        dim_r_data = 4;
    end

    if ~exist('dim_t_V','var')||isempty(dim_t_V)
        dim_t_V = 1;
    end
    
    % permute rank to the last dimension
    perm_order             = 1:max(ndims(cxr),dim_r_data);
    perm_order(end)        = dim_r_data;
    perm_order(dim_r_data) = max(ndims(cxr),dim_r_data);
    cxr_perm               = permute(cxr,perm_order); % [y,x,z,coil,avg,...,rank]

    % take care of the basis
    switch dim_t_V
        case 1
            Vrt = V_subspace.'; % [rxt]
        case 2
            Vrt = V_subspace;   % [rxt]
        otherwise
            error('The fourth input must be either 1 or 2.');
    end
    
%% analysis
    sxt_syn = reshape(cxr_perm,[],size(cxr,dim_r_data))*Vrt; % [Pxr] * [rxt] = [Pxt]
    sxt_syn = reshape(sxt_syn,[size_dims(cxr_perm,1:max(ndims(cxr),dim_r_data)-1),size(Vrt,2)]); % [y,x,z,coil,avg,...,t]
    sxt_syn = permute(sxt_syn,perm_order);
    
end

