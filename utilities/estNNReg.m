function lambda = estNNReg(d,var_est,fm_hz,tvec,mask)

    dims = size(mask);
    if length(dims) == 3
        Ny = dims(1);
        Nx = dims(2);
        Nz = 1;
        Nt = dims(3);
        mask = permute(mask,[1,2,4,3]); % pretend to be 3d
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

    % use semi-circular law to estimate value of largest noise singular value
    N = nnz(mask);
    noise_sv = sqrt(var_est)*(sqrt(N/Nt) + sqrt(Nt));
    
    % minimum norm reconstruction
    phi = reshape(exp(-2i*pi*fm_hz(:)*tvec),[Ny,Nx,Nz,Nt]);
    xt_mn = zeros(size(mask));
    xt_mn(mask) = d(:);
    xt_mn = F3_k2x(xt_mn).*phi;
    clear phi;
    
    % regularization term
    s = svd(reshape(xt_mn,Ny*Nx*Nz,Nt),'econ');
    L = find(s < noise_sv,1,'first');
    e = sum(s(1:L));
    
    % estimated lambda
    lambda = N*var_est/e;

end