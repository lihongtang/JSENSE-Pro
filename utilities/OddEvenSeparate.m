function [OddData,EvenData] = OddEvenSeparate(IntData,dim,fill_zero)
%% ==================================================================
%INTERLEAVEODDEVENDATA separate odd and even echoes of 4D MRSI data
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-05-16
%
%   [INPUTS]
%   ---- Required ----
%   IntData                 Interleaved EPSI data    [y,x,z,t_oddeven,(frame,coil,ave...)]
%
%   ---- Optional ----
%   dim                     the dimension along which two data are to be
%                           separated [4]
%   fill_zero               flag to fill zeros for odd and even data [false]
%
%   [OUTPUTS]
%   OddData                 EPSI data of odd echoes  [y,x,z,t_odd,(frame,coil,ave...)]
%   EvenData                EPSI data of even echoes [y,x,z,t_even,(frame,coil,ave...)]
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/05/16
%       Modified by Yibo Zhao @ UIUC, 2018/07/11:
%           Support zero-filling for logical matrix.
%       Modified by Yibo Zhao @ UIUC, 2018/09/29
%           In the comments, add a new method using subsref and subsasgn, 
%           which is more readable (perhaps), but slower.
%
%   Example:
%       IntData            = 1:10; % 1x10 matrix
%       [OddData,EvenData] = OddEvenSeparate(IntData,2,false); % do fill zeros
%       [OddData,EvenData] = OddEvenSeparate(IntData,2,true); % don't fill zeros
%
%   See also ODDEVENINTERLEAVE
%--------------------------------------------------------------------------

%% ------ parse the input ------
    if(nargin < 2)||isempty(dim);
        dim = 4;
    end
    
    if(nargin < 3)||isempty(fill_zero);
        fill_zero = false;
    end
    
    data_size  = size(IntData);
    
    try
        assert(mod(data_size(dim),2)==0);
    catch
        error('ERROR: the size along the dimension to separate must be even.');
    end
    
%% ------ separation ------
    if fill_zero
        if islogical(IntData)
            OddData        = false(data_size,'like',IntData);
            EvenData       = false(data_size,'like',IntData);
        else
            OddData        = zeros(data_size,'like',IntData);
            EvenData       = zeros(data_size,'like',IntData);
        end
        
        % --- odd data ---
        cmd_odd         = '(';
        for n=1:dim-1
            cmd_odd     = sprintf('%s:,',cmd_odd);
        end
        cmd_odd         = sprintf('%s1:2:end',cmd_odd);
        if(ndims(IntData)>dim)
            for n=dim+1:ndims(IntData)
                cmd_odd = sprintf('%s,:',cmd_odd);
            end
        end
        cmd_odd         = sprintf('%s)',cmd_odd);
        cmd_odd         = sprintf('OddData%s = IntData%s;',cmd_odd,cmd_odd);
        eval(cmd_odd); 
        % cmd_odd = 'OddData(:,:,:,1:2:end,:,...,:) = IntData(:,:,:,1:2:end,:,...,:);';
        
        % --- even data ---
        cmd_even         = '(';
        for n=1:dim-1
            cmd_even     = sprintf('%s:,',cmd_even);
        end
        cmd_even         = sprintf('%s2:2:end',cmd_even);
        if(ndims(IntData)>dim)
            for n=dim+1:ndims(IntData)
                cmd_even = sprintf('%s,:',cmd_even);
            end
        end
        cmd_even         = sprintf('%s)',cmd_even);
        cmd_even         = sprintf('EvenData%s = IntData%s;',cmd_even,cmd_even);
        eval(cmd_even);
        % cmd_even = 'EvenData(:,:,:,2:2:end,:,...,:) = IntData(:,:,:,2:2:end,:,...,:);';
        
    else
        % --- odd data ---
        cmd_odd         = 'OddData = IntData(';
        for n=1:dim-1
            cmd_odd     = sprintf('%s:,',cmd_odd);
        end
        cmd_odd         = sprintf('%s1:2:end',cmd_odd);
        if(ndims(IntData)>dim)
            for n=dim+1:ndims(IntData)
                cmd_odd = sprintf('%s,:',cmd_odd);
            end
        end
        cmd_odd         = sprintf('%s);',cmd_odd);
        eval(cmd_odd); 
        % cmd_odd = 'OddData = IntData(:,:,:,1:2:end,:,...,:);';
        
        % --- even data ---
        cmd_even         = 'EvenData = IntData(';
        for n=1:dim-1
            cmd_even     = sprintf('%s:,',cmd_even);
        end
        cmd_even         = sprintf('%s2:2:end',cmd_even);
        if(ndims(IntData)>dim)
            for n=dim+1:ndims(IntData)
                cmd_even = sprintf('%s,:',cmd_even);
            end
        end
        cmd_even         = sprintf('%s);',cmd_even);
        eval(cmd_even); 
        % cmd_even = 'EvenData = IntData(:,:,:,2:2:end,:,...,:);';
    end
    
    % another possible method (use subsref and subsasgn, but about a little bit slower than current method)
    %{
    if fill_zero
        if islogical(IntData)
            OddData        = false(data_size,'like',IntData);
            EvenData       = false(data_size,'like',IntData);
        else
            OddData        = zeros(data_size,'like',IntData);
            EvenData       = zeros(data_size,'like',IntData);
        end
        
        S.type             = '()';
        S.subs             = cell(ndims(IntData),1);
        for ii = 1:ndims(IntData)
            S.subs{ii}     = ':';
        end
    
        S.subs{dim}        = 1:2:data_size(dim);
        temp_odd           = subsref(IntData,S);
        OddData            = subsasgn(OddData,S,temp_odd);
    
        S.subs{dim}        = 2:2:data_size(dim);
        temp_even          = subsref(IntData,S);
        EvenData           = subsasgn(EvenData,S,temp_even);
    else
        S.type             = '()';
        S.subs             = cell(ndims(IntData),1);
        for ii = 1:ndims(IntData)
            S.subs{ii}     = ':';
        end
    
        S.subs{dim}        = 1:2:data_size(dim);
        OddData            = subsref(IntData,S);
    
        S.subs{dim}        = 2:2:data_size(dim);
        EvenData           = subsref(IntData,S);
    end
    %}
   
end
    
    
