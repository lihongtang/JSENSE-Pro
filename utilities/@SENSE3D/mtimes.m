function res = mtimes(a,x)
% This method applies the SENSE operator

maps = a.sen;
[sx,sy,sz,nc] = size(maps);

if a.adjoint    
    x = reshape(x,sx,sy,sz,nc);
%    res = sum(conj(maps).*x./sum(abs(maps).^2,4),4);    
res = sum(conj(maps).*x,4); 

%  res = sum(conj(maps).*x,4);
else
    x = reshape(x,[sx,sy,sz]);
   res = maps.*repmat(x,[1,1,1,nc]);   
  
       
end

