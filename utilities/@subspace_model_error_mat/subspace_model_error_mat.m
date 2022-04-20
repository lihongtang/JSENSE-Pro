
function  res =subspace_model_error_mat(img,b1,dim_spatial,undersampleMask)

%res = Subspace_mat(b1)
%
%
%	Implementation of multi-coil subspace projection matrix of sensitivity data
%	
%	input:
%			
%           b1 : coil sensitivity basis (Ny*Nx*Nz,Nc)
%
%	output: the operator
%
%	(c) Ricardo Otazo 2008

res.adjoint = 0;
res.b1 = b1;
res.img = img(:);
res.mask = undersampleMask;
res.dim_spatial =dim_spatial;
for i_coil = 1:size(b1,2)
    rank(i_coil) = size(b1(i_coil).subspace,2);
end
res.rank = rank;
res = class(res,'subspace_model_error_mat');


