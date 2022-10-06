function res = mtimes(a,b)

if a.adjoint
    % x=H'*y (x=res,b=y), x: object, y:multi-coil coil sensitivity map
    % multi-coil combination in the image domain
    b = reshape(Fn_k2x(reshape(b,a.dim_spatial).*a.mask,1:3),[],a.dim_spatial(4));
    for ch=1:size(a.b1,2)
        temp_basis = cell2mat(struct2cell((a.b1(ch))));
        coil_position = (sum(a.rank(1:(ch-1)),2)+1):(sum(a.rank(1:(ch)),2));
        res(coil_position,1) = temp_basis'*(conj(a.img).*b(:,ch));
        clearvars temp_basis; clearvars coil_position;
    end
else
    % y=H*x (y=res,b=x), x: object, y: multi-coil coil sensitivity map
    % a.b1:multi-coil basis from object;b is the desired parameters
    x_array = zeros([prod(a.dim_spatial(1:3)),a.dim_spatial(4)]);
    for ch=1:size(a.b1,2)
        temp_basis = cell2mat(struct2cell(a.b1(ch)));
        coil_position = (sum(a.rank(1:(ch-1)),2)+1):(sum(a.rank(1:(ch)),2));
        x_array(:,ch) = a.img.*(temp_basis*b(coil_position,1));
        clearvars temp_basis;        clearvars coil_position;
    end
    res = Fn_x2k(reshape(x_array,a.dim_spatial),1:3).*a.mask;
end