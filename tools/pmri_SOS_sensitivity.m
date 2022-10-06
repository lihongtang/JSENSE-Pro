function [output_mask,output] = pmri_SOS_sensitivity(input,mask)
%
%	pmri_SOS_sensitivity	estimating sensitivity maps using the division
%	of images from each channel by the sum-of-squares of the set
%
%	[output_mask,output]=pmri_poly_sensitivity(input,order,['mask', mask,....])
%
%	INPUT:
%	input: 2D/3D input image [n_PE, n_PE,r_RO,coil].
%		n_PE: # of phase encoding
%		n_FE: # of frequency encoding

%	mask: 2D/3D input image mask [n_PE, n_PE].
%		n_PE: # of phase encoding
%		n_FE: # of frequency encoding
%		"0" indicates the correponding entries are included.
%		"1" indicates the correponding entries are excluded.
%	'flag_display': value of either 0 or 1
%		It indicates of debugging information is on or off.
%
%	OUTPUT:
%	output_mask: 2D estimated coil map with a spatial mask [n_PE, n_PE].
%		n_PE: # of phase encoding steps
%		n_FE: # of frequency encoding steps
%	output_mask: 2D estimated coil map without a spatial mask [n_PE, n_PE].
%		n_PE: # of phase encoding steps
%		n_FE: # of frequency encoding steps
%
%---------------------------------------------------------------------------------------
%	Lihong@SJTU	
%	lihong@may. 26, 2021
%



dim = size(input);
if exist('mask','var')
else
     mask = ones(dim(1:(end-1)));
end
output=reshape(zeros(dim),[],dim(end));
output_mask=reshape(zeros(dim),[],dim(end));
input = reshape(input,[],dim(end));
mask = reshape(mask,[],1);
sos_all_c = sqrt(sum(abs(input).^2,2));
for i_c = 1:dim(end)
    output(:,i_c) = input(:,i_c)./sos_all_c;
        output_mask(:,i_c) = output(:,i_c).*mask;

end
output = reshape(output,dim);
output_mask = reshape(output_mask,dim);





