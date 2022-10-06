function sen_norm = senNormRange(sen)
% normlize sensitivity to sum(abs(coil).^2)=1
% sen: [Nx,Ny,Nz,coil]
sen_norm = sen./repmat(sqrt(sum(abs(sen).^2,4)),[1,1,1,size(sen,4)]);
end