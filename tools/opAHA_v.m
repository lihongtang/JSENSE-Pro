function y = opAHA_v(x, U, inhomo, L, Nx, Ny, Nz, Nt, mask)

    y = U*reshape(x,L,Nt);
    
%     y = opAHA_4D(y(:), inhomo, Nx, Ny, Nz, Nt, mask);
    y = reshape(y, Ny, Nx, Nz, Nt);
    y = inhomo.*y;
    y = F3_x2k(y,1);
    y(~mask) = 0;
    y = F3_k2x(y,1);
    y = conj(inhomo).*y;

    y = U'*reshape(y,Nx*Ny*Nz,Nt);
    y = y(:);

end