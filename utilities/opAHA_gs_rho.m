function y = opAHA_gs(x,b,inhomo,ref,mask,gsMask)
    y = opG_gs(x,ref,gsMask);
    y = opFHF_gs(y,inhomo,mask) + b*y/2;
    y = opGH_gs(y,ref,gsMask);
    return
end