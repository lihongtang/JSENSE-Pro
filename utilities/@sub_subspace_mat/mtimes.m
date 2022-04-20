function res = mtimes(a,b)

if a.adjoint
    % x=H'*y (x=res,b=y), x: object, y:multi-coil coil sensitivity map
    % multi-coil combination in the image domain
    b = b.*a.brainMaskRef;
    for ch=1:size(a.b1,2)
        temp_basis = cell2mat(struct2cell((a.b1(ch))));
        coil_position = (sum(a.rank(1:(ch-1)),2)+1):(sum(a.rank(1:(ch)),2));
        res(coil_position,1) = temp_basis'*b(:,ch);
        clearvars temp_basis; clearvars coil_position;
    end
else
    % y=H*x (y=res,b=x), x: object, y: multi-coil coil sensitivity map
    % a.b1:multi-coil basis from object;b is the desired parameters
    
    for ch=1:size(a.b1,2)
        temp_basis = cell2mat(struct2cell(a.b1(ch)));
        coil_position = (sum(a.rank(1:(ch-1)),2)+1):(sum(a.rank(1:(ch)),2));
        res(:,ch) = temp_basis*b(coil_position,1).*a.brainMaskRef(:,ch);
        clearvars temp_basis;        clearvars coil_position;
        
        
    end
    
end