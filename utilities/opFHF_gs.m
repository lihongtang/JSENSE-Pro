function y = opFHF_gs( x, inhomo, mask)
 
    y = inhomo.*reshape(x,size(inhomo));
    y = F3_x2k(y,1);
    y = opSampleH(opSample(y,mask),mask);
    y = F3_k2x(y,1);
    y = reshape(conj(inhomo).*y,[],1);
    
    return

end

