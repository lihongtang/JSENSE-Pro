function [ AHAx ] = opAHA_3D( img_vec, inhomo, Nx, Ny, Nt, mask)
%opAHA Operator for minimization
%--------------------------------------------------------------------------
% est_xyt = min( |A*x - d|^2 + reg_para*R(D*x) )
% A = C*F*B
% C <--> Sampling Operator
% F <--> FFT (x to k)
% B <--> Phase From Field Inhomogeneity
% D = Difference Operator
% R(.) = Regularization function
% d = data
%--------------------------------------------------------------------------

    img = reshape(img_vec, Ny, Nx, Nt);
    AHAx = inhomo.*img;
    AHAx = F_x2k(AHAx,1);
    AHAx = opSampleH(opSample(AHAx,mask),mask);
    AHAx = F_k2x(AHAx,1);
    AHAx = conj(inhomo).*AHAx;
    AHAx = AHAx(:);
    
    return

end

