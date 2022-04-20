function e = l2err(x,y)
% computes the error of two vectors in the sense of L2 norm
    if islogical(x)
        x = double(x);
    end
    
    if islogical(y)
        y = double(y);
    end
    
    e = norm(x(:) - y(:))/norm(x(:));
    return;

end