function y = opAHA_u(x, V, inhomo, L, Nx, Ny, Nz, Nt, mask)

    y = reshape(x,Nx*Ny*Nz,L)*V;
    
%     y = opAHA_4D(y(:), inhomo, Nx, Ny, Nz, Nt, mask);
    y = reshape(y, Ny, Nx, Nz, Nt);
    y = inhomo.*y;
    y = F3_x2k(y,1);
    y(~mask) = 0;
    y = F3_k2x(y,1);
    y = conj(inhomo).*y;

    y = reshape(y,Nx*Ny*Nz,Nt)*V';
    y = y(:);

end