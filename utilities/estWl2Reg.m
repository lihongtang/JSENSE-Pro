function lambda = estWl2Reg(var_est,fm_hz,tvec,w,dirs,mask)

    dims = size(mask);
    if length(dims) == 3
        Ny = dims(1);
        Nx = dims(2);
        Nz = 1;
        Nt = dims(3);
        mask = permute(mask,[1,2,4,3]); % pretend to be 3d
        w = permute(w,[1,2,5,3,4]);
    elseif length(dims) == 4
        Ny = dims(1);
        Nx = dims(2);
        Nz = dims(3);
        Nt = dims(4);
    else
        error('bad data size');
    end
    assert(isvector(tvec));
    assert(length(tvec) == Nt);
    assert(length(dirs) == size(w,5));

    % generate fake noise with the estimated variance
    N = nnz(mask);
    noise = ( randn(N,1) + 1i*randn(N,1) )*sqrt(var_est/2);
    
    % minimum norm reconstruction
    phi = reshape(exp(2i*pi*fm_hz(:)*tvec),[Ny,Nx,Nz,Nt]);
    xt_mn = zeros(size(mask));
    xt_mn(mask) = noise;
    xt_mn = F3_k2x(xt_mn).*phi;
    clear phi;
    
    % regularization term
    e = real(.5*xt_mn(:)'*opDWWD_4D(xt_mn(:), Nx, Ny, Nz, Nt, dirs, sqrt(w)));
    
    % estimated lambda
    lambda = N*var_est/e;

end