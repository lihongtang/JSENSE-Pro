function [ AHAx ] = opAHA_4D( img_vec, inhomo, Nx, Ny, Nz, Nt, mask)
%opAHA Operator for minimization
%--------------------------------------------------------------------------
% A = C*F*B
% C <--> Sampling Operator
% F <--> FFT (x to k)
% B <--> Phase From Field Inhomogeneity
% R(.) = Regularization function
% d = data
%--------------------------------------------------------------------------

    img = reshape(img_vec, Ny, Nx, Nz, Nt);
    AHAx = inhomo.*img;
    AHAx = F3_x2k(AHAx,1);
    AHAx = opSampleH(opSample(AHAx,mask),mask);
    AHAx = F3_k2x(AHAx,1);
    AHAx = conj(inhomo).*AHAx;
    AHAx = AHAx(:);
    
    return

end

