function [ DWWDx ] = opDWWD_4D( x, Nx, Ny, Nz, Nt, dirs, w)
%opDWWD Operator for minimization
%--------------------------------------------------------------------------
% est_x = min( |A*x - d|^2 + reg_para*R(W*D*x) )
% A = Forward Operator
% D = Difference Operator
% W = Diagonal Weighting Operator
% R(.) = Regularization function
% d = data
%--------------------------------------------------------------------------

    x = reshape(x,[Ny, Nx, Nz, Nt]);
    DWWDx = zeros(Ny, Nx, Nz, Nt);
    for n=1:length(dirs)
        DWWDx = DWWDx + diff_along_dir((w(:,:,:,:,n).^2).*diff_along_dir(x,dirs(n)),-dirs(n));
    end
    DWWDx = DWWDx(:);

    return

end

