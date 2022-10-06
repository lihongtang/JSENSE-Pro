function [U,Z,npenalty_nn] = local_low_rank_sen(temp_sen,reg_para_nn,b,Z,trunc_opt,dim_spatial)
temp_sen_test = reshape(temp_sen(1:prod(dim_spatial(1:3)),:),dim_spatial);
temp_Z_test = reshape(Z(1:prod(dim_spatial(1:3)),:),dim_spatial);

C = mat2cell(temp_sen_test,[96,96,96,96],[96,96,96,96],[ones(1,dim_spatial(3))],dim_spatial(4));
C_Z = mat2cell(temp_Z_test,[96,96,96,96],[96,96,96,96],[ones(1,dim_spatial(3))],dim_spatial(4));

U_cell = cell(size(C));
Z_cell = cell(size(C));
dim_patch = size(C{1,1,1});
npenalty_nn = 0;
for i_patch = 1:size(C,1)
    for j_patch = 1:size(C,2)
        for k_patch = 1:size(C,3)
            temp_sen_patch = C{i_patch,j_patch,k_patch};
            temp_sen_patch = reshape(temp_sen_patch,[],size(temp_sen_patch,4));
            temp_Z_patch = C_Z{i_patch,j_patch,k_patch};
            temp_Z_patch = reshape(temp_Z_patch,[],size(temp_Z_patch,4));            
            U_cell{i_patch,j_patch,k_patch} = rank_trunc( temp_sen_patch- temp_Z_patch/b, reg_para_nn/b, trunc_opt);            
            Z_cell{i_patch,j_patch,k_patch} = temp_Z_patch + b*(U_cell{i_patch,j_patch,k_patch} -temp_sen_patch );
            U_cell{i_patch,j_patch,k_patch} = reshape(U_cell{i_patch,j_patch,k_patch},dim_patch);
            Z_cell{i_patch,j_patch,k_patch} = reshape(Z_cell{i_patch,j_patch,k_patch},dim_patch);
        npenalty_nn     = npenalty_nn+reg_para_nn*sum(abs(svd(temp_sen_patch)));

        end
    end
end

U = cell2mat(U_cell);
U = reshape(U,[prod(dim_spatial(1:3)),dim_spatial(4)]);

Z = cell2mat(Z_cell);
Z = reshape(Z,[prod(dim_spatial(1:3)),dim_spatial(4)]);
end