function y = DWWD_v( x, U, L, Nx, Ny, Nz, Nt, dirs, w)

    y = U*reshape(x,L,Nt);
    
%     y = opDWWD_4D( y(:), Nx, Ny, Nz, Nt, dirs, w);
    y = reshape(y,[Ny, Nx, Nz, Nt]);
    DWWDx = zeros(Ny, Nx, Nz, Nt);
    for n=1:length(dirs)
        DWWDx = DWWDx + diff_along_dir((w(:,:,:,:,n).^2).*diff_along_dir(y,dirs(n)),-dirs(n));
    end
    
    y = U'*reshape(DWWDx,Nx*Ny*Nz,Nt);
    y = y(:);
    
end