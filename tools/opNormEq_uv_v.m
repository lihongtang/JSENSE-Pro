function y = opNormEq_uv_v(x, U, inhomo, L, lambda, nu, Nx, Ny, Nz, Nt, mask, dirs, w)

    y = U*reshape(x,L,Nt);
    y = reshape(y, Ny, Nx, Nz, Nt);
    
    if lambda > 0
        % Forward model:
        % A'*A
        % A = (sampling)(Fourier transform)(field inhomogeneity phase)
        AHAx = inhomo.*y;
        AHAx = F3_x2k(AHAx,1);
        AHAx(~mask) = 0;
        AHAx = F3_k2x(AHAx,1);
        AHAx = conj(inhomo).*AHAx;

        % Weighted L2 regularization on UV
        % D'*W*W*D
        DWWDx = zeros(Ny, Nx, Nz, Nt);
        for n=1:length(dirs)
            DWWDx = DWWDx + diff_along_dir((w(:,:,:,:,n).^2).*diff_along_dir(y,dirs(n)),-dirs(n));
        end
        
        y = U'*reshape(AHAx + .5*lambda*DWWDx,Nx*Ny*Nz,Nt);
        
    else
        % Forward model:
        % A'*A
        % A = (sampling)(Fourier transform)(field inhomogeneity phase)
        y = inhomo.*y;
        y = F3_x2k(y,1);
        y(~mask) = 0;
        y = F3_k2x(y,1);
        y = conj(inhomo).*y;
        
        y = U'*reshape(y,Nx*Ny*Nz,Nt);
        
    end
    
    % L2 regularization on U
    y = y(:) + .5*nu*x;
    
    return;
    
end