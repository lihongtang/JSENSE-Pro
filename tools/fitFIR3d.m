function [h12,dk12] = fitFIR3d(dk1,dk2,dims_gs,os,lambda,D,Omega)
%% ==================================================================
%FITFIR3D fits 3D circular FIR to compensate two k/f-space signals
%Model: h12 = argmin_h ||dk2 - Omega{dk1.*f(h)}|| + lambda*||D*h||
%       where f[m]=\sum_{n} h[n]*exp{i2pi*n*m/M}, m=0,...,M-1 for each dimension
%       dk12 = dk1.*f(h)
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2019-06-21
%
%   [INPUTS]
%   ---- Required ----
%   dk1                     dk1(ky,kx,kz)
%   dk2                     dk2(ky,kx,kz)
%   dims_gs                 bilateral order of generalized series
%
%   ---- Optional ----
%   os                      oversampling factor [1]
%   lambda                  regularization parameter [0]
%   D                       regularization weighting [1]
%   Omega                   sampling mask [true(size(sx2))]
%
%   [OUTPUTS]
%   h12                     generalized series coefficients (size: dims_gs)
%   dk12                    compensated function (ky,kx,kz)
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2019/06/21
%
%   See also FITGS1D, APPGS1D, FITGS2D, APPGS2D
%
%   Example:
%       GS_order = [9,9,1]; % 3D GS order, better to be odd numbers
%       [h12,dk12] = fitFIR3d(F3_x2k(I1),F3_x2k(I2),GS_order); % I1: reference, I2: target
%       I12 = F3_k2x(dk12); % I12: compensated I1
%
%--------------------------------------------------------------------------

%% ------ parse the input ------

    [Ny1,Nx1,Nz1] = size_dims(dk1,1:3);
    [Ny2,Nx2,Nz2] = size_dims(dk2,1:3);
    
    dims_gs = vec(dims_gs);
    assert(length(dims_gs)<=3,'GS order cannot be higher than 3 for 3D image.');
    dims_gs = cat(1,dims_gs,ones(3-length(dims_gs),1));
    
    if ~exist('os','var') || isempty(os)
        os = 1;
    end
    
    if ~exist('lambda','var') || isempty(lambda)
        lambda = 0;
    end
    
    if ~exist('Omega','var') || isempty(Omega)
        Omega = true([Ny2,Nx2,Nz2]);
    end
    
    if ~exist('D','var') || isempty(D)
        scale_fac = abs(dk1(cenInd(Ny1),cenInd(Nx1),cenInd(Nz1)));
        D = scale_fac*eye(prod(dims_gs));
    end
    
%% ------ fitting coefficients ------
    
    if nargout<=1
        linear_A = zeros([Ny2,Nx2,Nz2,dims_gs.'],'like',dk1);
    else
        linear_A = zeros([Ny1,Nx1,Nz1,dims_gs.'],'like',dk1);
    end
    
    for ii = 1:dims_gs(1)
        for jj = 1:dims_gs(2)
            for kk = 1:dims_gs(3)
                if isequal(os,1)
                    shift1 = ii-(dims_gs(1)+1)/2;
                    shift2 = jj-(dims_gs(2)+1)/2;
                    shift3 = kk-(dims_gs(3)+1)/2;
                    temp = circshift(dk1,[shift1,shift2,shift3]);
                else
                    shift1 = ii-(dims_gs(1)+1)/2;
                    shift2 = jj-(dims_gs(2)+1)/2;
                    shift3 = kk-(dims_gs(3)+1)/2;
                    temp = posShift3DinXspace(dk1,shift1/os,shift2/os,shift3/os);
                end
                if nargout<=1
                    temp = temp(cenInd(Ny1,Ny2),cenInd(Nx1,Nx2),cenInd(Nz1,Nz2));
                end
                linear_A(:,:,:,ii,jj,kk) = temp;
            end
        end
    end
    if nargout<=1
        linear_A = linear_A(:,:,:,:);
    else
        linear_A_full = reshape(linear_A(:,:,:,:),[Ny1*Nx1*Nz1,prod(dims_gs)]);
        linear_A = linear_A(cenInd(Ny1,Ny2),cenInd(Nx1,Nx2),cenInd(Nz1,Nz2),:);
    end
    linear_A = linear_A(repmat(Omega,[1,1,1,size(linear_A,4)]));
    linear_A = reshape(linear_A,[nnz(Omega),prod(dims_gs)]);
    
    % regularization
    linear_A = cat(1,linear_A,sqrt(lambda)*D);
    linear_b = cat(1,vec(dk2(Omega)),zeros(prod(dims_gs),1));
    
    % LS fitting
    h12      = pinv(linear_A)*linear_b;
    h12      = reshape(h12,dims_gs.');
   
    % apply fitted coefficients
    if nargout>1
        dk12 = reshape(linear_A_full*vec(h12),[Ny1,Nx1,Nz1]);
    end
end

%{
    % full forward matrix
    linear_A = full(convmtx2(sx1,N_gsx,N_gsy));
    
    % reshape and truncation
    linear_A = reshape(linear_A,[Ny1+N_gsx-1,Nx1+N_gsy-1,N_gsx*N_gsy]);
    linear_A = linear_A(cenInd(Ny1+N_gsx-1,Ny2),cenInd(Nx1+N_gsy-1,Nx2),:);
    linear_A = linear_A(repmat(Omega,[1,1,size(linear_A,3)]));
    
    % truncated forward matrix
    linear_A = reshape(linear_A,[nnz(Omega),N_gsx*N_gsy]);
%}



