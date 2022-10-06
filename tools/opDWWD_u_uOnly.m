function DWWDx = DWWD_u_uOnly( x, L, Nx, Ny, Nz, Nt, dirs, w)

    y = reshape(x,[Ny, Nx, Nz, L]);
    DWWDx = zeros(Ny, Nx, Nz, L);
    for n=1:length(dirs)
        DWWDx = DWWDx + diff_along_dir((w(:,:,:,1:L,n).^2).*diff_along_dir(y,dirs(n)),-dirs(n));
    end
    DWWDx = DWWDx(:);
    
end