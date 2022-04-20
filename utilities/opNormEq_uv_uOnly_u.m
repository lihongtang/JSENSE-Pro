function y = opNormEq_uv_uOnly_u(x, V, inhomo, L, lambda, mu, Nx, Ny, Nz, Nt, mask, dirs, w)

    if lambda > 0
        
        % Forward model:
        % A'*A
        % A = (sampling)(Fourier transform)(field inhomogeneity phase)
        y = reshape(x,Nx*Ny*Nz,L)*V;
        y = inhomo.*reshape(y, Ny, Nx, Nz, Nt);
        y = F3_x2k(y,1);
        y(~mask) = 0;
        y = F3_k2x(y,1);
        y = conj(inhomo).*y;
        y = reshape(y,Nx*Ny*Nz,Nt)*V';

        % Weighted L2 regularization on UV
        % D'*W*W*D
        U = reshape(x, Ny, Nx, Nz, L);
        DWWDx = zeros(Ny, Nx, Nz, L);
        for n=1:length(dirs)
            DWWDx = DWWDx + diff_along_dir((w(:,:,:,:,n).^2).*diff_along_dir(U,dirs(n)),-dirs(n));
        end
        
        y = y(:) + .5*lambda*DWWDx(:);
        
    else
        % Forward model:
        % A'*A
        % A = (sampling)(Fourier transform)(field inhomogeneity phase)
        y = reshape(x,Nx*Ny*Nz,L)*V;
        y = reshape(y, Ny, Nx, Nz, Nt);
        
        y = inhomo.*y;
        y = F3_x2k(y,1);
        y(~mask) = 0;
        y = F3_k2x(y,1);
        y = conj(inhomo).*y;
        
        y = reshape(y,Nx*Ny*Nz,Nt)*V';
        
    end
    
    % L2 regularization on U
    y = y(:) + .5*mu*x;
    
    return;
    
end