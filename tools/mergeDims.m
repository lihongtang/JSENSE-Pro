function x = mergeDims(x,dims_cell)
% MERGEDISM merge dimensions of a matrix
%   Created by Yibo Zhao @ UIUC, 2019/10/14
%
% Example:
%   Say A is a 5D matrix of dimension [32,32,12,74,16], then
%   A_reshape = mergeDims(A,{1:3,4,5}); will return a 3D matrix of dimension [32x32x12,74,16]
%   A_csrt    = mergeDims(A,{[1:3,5],4}); will return a 2D Casorati matrix of dimension [32x32x12x16,74]
%
% Note: this function can be viewed as a generalized PERMUTE function
%
%% parse inputs
    dims_vec = cell2vec(dims_cell);
    if length(dims_vec)<ndims(x)
%         dims_res = setdiff(1:ndims(x),dims_vec);
%         dims_cell = cat(1,dims_cell,{dims_res});
%         dims_vec = cell2vec(dims_cell);
        error('dims_cell [%s] does not match the data [%s].',...
            gen_str_from_vec(dims_vec,','),gen_str_from_vec(1:ndims(x),','));
    end
    assert(length(unique(dims_vec))==length(dims_vec),'dims_cell cannot contain repeated values.');
    
%% preparation: size of data
    size_var = zeros([1,length(dims_cell)]);
    for ii = 1:length(dims_cell)
        size_var(ii) = numel_dims(x,dims_cell{ii});
    end

%% permute and reshape
    x = permute(x,dims_vec);
    x = reshape(x,size_var);
    
end

%% nested function: convert cell to matrix
function vec_out = cell2vec(cell_in)
    vec_out = [];
    for ii = 1:length(cell_in)
        vec_out = cat(1,vec_out,vec(cell_in{ii}));
    end
end