function y = opG_gs( x, ref, gsMask)

    y = opSampleH(x,gsMask);
    y = F3_k2x(y,1);  
    y = reshape(ref.*y,[],1);
    
    return

end

