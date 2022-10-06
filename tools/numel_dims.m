function N = numel_dims(data,dims) 
% Compute number of elements in selected dimensions
%   Created by Yibo Zhao @ UIUC, 07/23/2018
%   Modified by Yibo Zhao @ UIUC, 09/17/2018:
%     Use prod(size_dims(...)) to accelerate.

    N = prod(size_dims(data,dims));
    
end

% previous method
%{
    cmd = sprintf('N = numel(data(');
    for n = 1:ndims(data)
        if n ~= ndims(data)
            if any(n == dims)
                cmd = sprintf('%s:,',cmd);
            else
                cmd = sprintf('%s1,',cmd);
            end
        else
            if any(n == dims)
                cmd = sprintf('%s:',cmd);
            else
                cmd = sprintf('%s1',cmd);
            end
        end
    end
    cmd = sprintf('%s));',cmd);
    eval(cmd);
    
    e.g.  for N = numel_dims(data,[1,2,3]), data is 4-D, actual 
    computation is N = numel(data(:,:,:,1))

%#ok<*STOUT>
%}
