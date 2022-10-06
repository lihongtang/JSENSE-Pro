function y = opNormEq_uv_uOnly_v(x, U, inhomo, L, lambda, nu, Nx, Ny, Nz, Nt, mask)

    y = U*reshape(x,L,Nt);
    y = reshape(y, Ny, Nx, Nz, Nt);
    
    % Forward model:
    % A'*A
    % A = (sampling)(Fourier transform)(field inhomogeneity phase)
    y = inhomo.*y;
    y = F3_x2k(y,1);
    y(~mask) = 0;
    y = F3_k2x(y,1);
    y = conj(inhomo).*y;

    y = U'*reshape(y,Nx*Ny*Nz,Nt);
    
    % L2 regularization on U
    y = y(:) + .5*nu*x;
    
    return;
    
end