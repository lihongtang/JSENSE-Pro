function IntData = OddEvenInterleave(OddData,EvenData,dim)
%% ==================================================================
%ODDEVENINTERLEAVE interleave odd and even echoes of 4D MRSI data
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-04-30
%
%   [INPUTS]
%   ---- Required ----
%   OddData                 EPSI data of odd echoes  [y,x,z,t_odd,(frame,coil,ave...)]
%   EvenData                EPSI data of even echoes [y,x,z,t_even,(frame,coil,ave...)]
%
%   ---- Optional ----
%   dim                     the dimension along which two data are to be
%                           interleaved [4]
%
%   [OUTPUTS]
%   IntData                 Interleaved EPSI data    [y,x,z,t_oddeven,(frame,coil,ave...)]
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/04/30
%       Modified by Yibo Zhao @ UIUC, 2018/05/16
%           used a trick to break the limit on the number of dimensions, 
%           and enabled interleaving data along any dimension.
%       Modified by Yibo Zhao @ UIUC, 2018/09/29
%           in the comments, add a new method using subsasgn, which is more
%           readable (perhaps), but slower.
%
%   Example:
%       OddData  = 1:2:10; % 1x5 matrix
%       EvenData = 2:2:10; % 1x5 matrix
%       IntData  = OddEvenInterleave(OddData,EvenData,2); % 1x10 matrix
%
%   See also ODDEVENSEPARATE
%--------------------------------------------------------------------------

%% ------ parse the input ------

    if(nargin <3);
        dim = 4;
    end
    
    data_size_odd    = size(OddData);
    data_size_even   = size(EvenData);
    
    try
        assert( all(data_size_odd == data_size_even) );
    catch
        error('ERROR: Odd and even data must have the same size.');
    end
    
    data_size        = data_size_odd;
    data_size(dim)   = data_size_odd(dim) + data_size_even(dim);
    
    IntData          = zeros(data_size,'like',OddData);
    
%% --- a trick (updated on 2018/05/16) ---
    % --- odd data ---
    cmd_odd          = 'IntData(';
    for n=1:dim-1
        cmd_odd      = sprintf('%s:,',cmd_odd);
    end
    cmd_odd          = sprintf('%s1:2:end',cmd_odd);
    if(ndims(OddData)>dim)
        for n=dim+1:ndims(OddData)
            cmd_odd  = sprintf('%s,:',cmd_odd);
        end
    end
    cmd_odd          = sprintf('%s) = OddData;',cmd_odd);
    eval(cmd_odd); 
    % Explain: cmd_odd = 'IntData(:,:,:,1:2:end,:,...,:) = OddData';
    
    % --- even data ---
    cmd_even         = 'IntData(';
    for n=1:dim-1
        cmd_even     = sprintf('%s:,',cmd_even);
    end
    cmd_even         = sprintf('%s2:2:end',cmd_even);
    if(ndims(OddData)>dim)
        for n=dim+1:ndims(OddData)
            cmd_even = sprintf('%s,:',cmd_even);
        end
    end
    cmd_even         = sprintf('%s) = EvenData;',cmd_even);
    eval(cmd_even); 
    % Explain: cmd_even = 'IntData(:,:,:,2:2:end,:,...,:) = EvenData';
    
    
    % former method (kind of hard-code)
    %{
    if ndims(OddData)>8
        warning('Too many dimensions. Typical number dimensions should be less than 8.');
    end
    
    data_size = data_size_odd;
    data_size(4) = data_size_odd(4) + data_size_even(4);
    
    IntData = zeros(data_size,'like',OddData);
    IntData(:,:,:,1:2:end,:,:,:,:) = OddData;
    IntData(:,:,:,2:2:end,:,:,:,:) = EvenData;
    %}
    
    % another possible method (use subsasgn, but about two times slower than current method)
    %{
    S.type         = '()';
    S.subs         = cell(ndims(OddData),1);
    for ii = 1:ndims(OddData)
        S.subs{ii} = ':';
    end
    
    S.subs{dim}    = 1:2:data_size(dim);
    IntData        = subsasgn(IntData,S,OddData);
    
    S.subs{dim}    = 2:2:data_size(dim);
    IntData        = subsasgn(IntData,S,EvenData);
    %}
    
end
    
    
