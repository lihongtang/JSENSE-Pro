function ppm_vec = gen_ppm_vec(Nt,dt,f0,ref_ppm) 
%% ==================================================================
%GEN_PPM_VEC generate a vector of ppm (parts per million) values
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-10-22
%
%   [INPUTS]
%   ---- Required ----
%   Nt                      number of temporal points
%   dt                      temporal interval (in second)
%
%   ---- Optional ----
%   f0                      systen frequency (in Hz) [123.2*1e6]
%   ref_ppm                 reference frequency (in ppm) [4.7]
%
%   [OUTPUTS]
%   ppm_vec                 ppm vector
%
%   Example:
%       get_voxel_spectrum3d(F3_t2f(sxt_mrsi),gen_ppm_vec(size(sxt_mrsi,4),dt_mrsi));
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/10/22
%       Modified by Yibo Zhao @ UIUC, 2018/10/23
%           Add an optional input: ref_ppm
%
%   Note: if the unit of Hz is wanted, use
%       gen_ppm_vec(Nt,dt,1e6,0);
%
%--------------------------------------------------------------------------

%%
    if ~exist('f0','var')||isempty(f0)
        f0 = 123.2*1e6;
    end

    if ~exist('ref_ppm','var')||isempty(ref_ppm)
        ref_ppm = 4.7;
    end

    ppm_vec = [-Nt/2:Nt/2-1]/Nt/dt; % in Hz

    ppm_vec = ppm_vec/f0*1e6;       % in ppm (parts per million)

    ppm_vec = ppm_vec + ref_ppm;    % with reference

end


