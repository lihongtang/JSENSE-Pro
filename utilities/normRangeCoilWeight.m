function dk_n = normRangeCoilWeight(dk,coil_weight)
dims = size(dk);
dk = reshape(dk,[],size(coil_weight,2));
dk_n = zeros(size(dk));

    for i_coil=1:size(coil_weight,2)
        temp = dk(:,i_coil);
       dk_n(:,i_coil) = temp./max(abs(temp(:))) .*coil_weight(1,i_coil);
    end
    dk_n = reshape(dk_n,dims);
end