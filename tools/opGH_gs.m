function y = opAH_gs( x, ref, gsMask)

    y = conj(ref).*reshape(x,size(ref));
    y = F3_x2k(y,1);
    y = opSample(y,gsMask);
    
    return

end

