function y = DWWD_u( x, V, L, Nx, Ny, Nz, Nt, dirs, w)

    y = reshape(x,Nx*Ny*Nz,L)*V;
    
%     y = opDWWD_4D( y(:), Nx, Ny, Nz, Nt, dirs, w);
    y = reshape(y,[Ny, Nx, Nz, Nt]);
    DWWDx = zeros(Ny, Nx, Nz, Nt);
    for n=1:length(dirs)
        DWWDx = DWWDx + diff_along_dir((w(:,:,:,:,n).^2).*diff_along_dir(y,dirs(n)),-dirs(n));
    end
    
    y = reshape(DWWDx,Nx*Ny*Nz,Nt)*V';
    y = y(:);
    
end