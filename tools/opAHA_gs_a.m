function y = opAHA_gs(x,b,inhomo,ref,mask,gsMask)
    y = opG_gs(x,ref,gsMask);
    y = opFHF_gs(y,inhomo,mask);
    y = opGH_gs(y,ref,gsMask);
    y = y + b*x/2;
    return
end