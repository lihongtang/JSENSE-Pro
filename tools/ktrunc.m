function dktc_out = ktrunc(dktc_in,newdims,doNorm,doPad)
%% ========================================================================
% KTRUNC zero-pads or truncates a k-space data. If doNorm, the intensity in
% image domain is unchanged. If doPad and newdim<size(dktc_in), zero-pad the
% output to be the same size of dktc_in.
% =========================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-05-17
%
%   [INPUTS]
%   ---- Required ----
%   dktc_in                       k-space data to be resized
%   newdim                        new dimension
%
%   ---- Optional ----
%   doNorm                        do normalization to keep image domain intensity [false]
%   doPad                         do zero-padding if newdim is smaller [false]
%
%   [OUTPUTS]
%   dktc_out                      resized k-space data
%
%   Change log:
%       Created by Yibo Zhao @ UIUC, 2018/05/17
%       Modified by Yibo Zhao @ UIUC, 2018/07/11
%           Support trucation/padding of logical matrices
%       Modified by Yibo Zhao @ UIUC, 2018/07/25
%           Fix an error when newdims is longer than the number of dimensions of dktc_in
%
%   Example:
%     temp  = phantom(256);
%     d_in  = F3_k2x(temp);
%     d_out = ktrunc(d_in,[200,300],true,false);
%     figure;montagesc(abs(F3_k2x(d_in)));caxis([0,1]);
%     figure;montagesc(abs(F3_k2x(d_out)));caxis([0,1]);
%--------------------------------------------------------------------------

    %% --- parse inputs ---
    if nargin < 4
        doPad = false;
        if nargin < 3
            doNorm = false;
        end
    end
    
    dims = size(dktc_in);
    numNewdim = length(newdims);
    
    if ~(length(dims)==numNewdim)
            dims(length(dims)+1:numNewdim) = 1;
    end
    
    if(doPad == 1)
        try
            assert(all(dims(1:numNewdim) >= newdims))
        catch
            error('ERROR: when using doPad, the data should be truncated.');
        end
    end
    
    if ~exist('cenInd','file')
        cenInd = @(Nd,Cd) (floor(Nd/2)+1+ceil(-Cd/2):floor(Nd/2)+ceil(Cd/2)); %#ok<NASGU>
    end
    
    %% --- pre-allocate ---
    if(doPad == 1)
        if islogical(dktc_in)
            dktc_out = false(dims);
        else
            dktc_out = zeros(dims,'like',dktc_in);
        end
    else
        if islogical(dktc_in)
            dktc_out = false([newdims,dims(numNewdim+1:end)]);
        else
            dktc_out = zeros([newdims,dims(numNewdim+1:end)],'like',dktc_in);
        end
    end
    
    %% --- a trick ---
    % Explanation of 3 modes:
    % 1. zero-pad:        dktc_out(cenInd(N_new,N_old),:,...) = dktc_in(:,:,...);
    % 2. truncate:        dktc_out(:,:,...) = dktc_in(cenInd(N_old,N_new),:,...);
    % 3. truncate, doPad: dktc_out(cenInd(N_old,N_new),:,...) = dktc_in(cenInd(N_old,N_new),:,...);
    
    % --- 1. things in the first parentheses ---
    cmd1         = '(';
    for n=1:numNewdim-1
        if dims(n) >= newdims(n) % truncate
            if(doPad == 1)
                cmd1 = sprintf('%scenInd(dims(%s),newdims(%s)),',cmd1,int2str(n),int2str(n));
            else
                cmd1 = sprintf('%s:,',cmd1);
            end
        else % zero-pad
            cmd1 = sprintf('%scenInd(newdims(%s),dims(%s)),',cmd1,int2str(n),int2str(n));
        end
    end
    if dims(numNewdim) >= newdims(numNewdim) % truncate
        if(doPad == 1)
            cmd1    = sprintf('%scenInd(dims(%s),newdims(%s))',cmd1,int2str(numNewdim),int2str(numNewdim));
        else
            cmd1    = sprintf('%s:',cmd1);
        end
    else % zero-pad
        cmd1         = sprintf('%scenInd(newdims(%s),dims(%s))',cmd1,int2str(numNewdim),int2str(numNewdim));
    end
    for n=numNewdim+1:ndims(dktc_in)
        cmd1     = sprintf('%s,:',cmd1);
    end
    cmd1         = sprintf('%s)',cmd1);
    
    % --- 2. things in the second parentheses ---
    cmd2         = '(';
    for n=1:numNewdim-1
        if dims(n) >= newdims(n) % truncate
            cmd2 = sprintf('%scenInd(dims(%s),newdims(%s)),',cmd2,int2str(n),int2str(n));
        else % zero-pad
            cmd2 = sprintf('%s:,',cmd2);
        end
    end
    if dims(numNewdim) >= newdims(numNewdim) % truncate
        cmd2         = sprintf('%scenInd(dims(%s),newdims(%s))',cmd2,int2str(numNewdim),int2str(numNewdim));
    else % zero-pad
        cmd2         = sprintf('%s:',cmd2);
    end
    for n=numNewdim+1:ndims(dktc_in)
        cmd2     = sprintf('%s,:',cmd2);
    end
    cmd2         = sprintf('%s)',cmd2);
    
    % --- 3. put 'em together ---
    cmd         = sprintf('dktc_out%s = dktc_in%s;',cmd1,cmd2);
    eval(cmd);
    
    %% --- finalize ---
    if(doNorm == 1)&&(doPad == 0) % if zero-filling, no need to normalize
        dktc_out = dktc_out*sqrt(prod(newdims)/prod(dims));
    end
    
    
end
