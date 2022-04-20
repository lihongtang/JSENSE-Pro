function mask_new = changeMaskMannually(img,mask_orig,mode,verbose)
%% ==================================================================
%changeMaskMannually lets you mannually change a mask
% ===================================================================
%   Author: Rong Guo @ UIUC
%   Created: 2018-10-09
%
%   [INPUTS]
%   ---- Required ----
%   img                     image to generate the mask
%
%   ---- Optional ----
%   mask_orig               original mask [false(size(img))]
%   mode                    what to do to the original mask ['union']
%   verbose                 degree of verbosity
%
%   [OUTPUTS]
%   mask_new                new mask
%
%   Change log:
%       Created by   Rong Guo @ UIUC, 2018/10/09
%           original name: genMask3DManually
%       Modified by Yibo Zhao @ UIUC, 2018/10/20
%           a) make this function independent
%           b) add original mask and and modification modes as inputs
%
%   For example, generating a mask can be done by repeatedly using following three lines:
%       mask = false(size_dims(image,[1,2,3]));        % initialize or clear to zero
%       mask = changeMaskMannually(image,mask,'|',1);  % add an area
%       mask = changeMaskMannually(image,mask,'&~',1); % subtract an area
%   until you get a satisfying mask.
%
%   See also genMask3DManually
%--------------------------------------------------------------------------

%% --- parse the inputs ---
    if ~exist('mask_orig','var') || isempty(mask_orig)
        mask_orig = false(size(img));
    end
    
    if ~exist('mode','var') || isempty(mode)
        mode = 'union';
    end
    
    if ~exist('verbose','var') || isempty(verbose)
        verbose = 0;
    end
    
%% --- mannually generate the marginal mask ---
    figure;fh = colorImageOverlay_local(abs(img),mask_orig,0.5);
    
    try
        h         = imfreehand;
        mask_marg = h.createMask(fh); % marginal mask
    catch
        warning('Forced quit! no change is made to the original mask.');
        mask_new  = mask_orig;
        return;
    end
%% --- change the original mask ---
    mask_new  = false(size_dims(img,1:3));
    
    for iy = 1:size(mask_marg,1)
        for ix = 1:size(mask_marg,2)
            [nx,ny,nz] = montagesub2imgsub_local(ix,iy,permute(size_dims(mask_new,[1,2,3]),[2,1,3]));
            switch mode
                case {'union','|','||','add'}        % add area to mask
                    mask_new(ny,nx,nz) = mask_orig(ny,nx,nz)||mask_marg(iy,ix);
                case {'intersect','&','&&','select'} % select area in mask
                    mask_new(ny,nx,nz) = mask_orig(ny,nx,nz)&&(mask_marg(iy,ix));
                case {'setdiff','&~','&&~','minus'}  % subtract area from mask
                    mask_new(ny,nx,nz) = mask_orig(ny,nx,nz)&&(~mask_marg(iy,ix));
                otherwise
                    error('More modes coming soon.');
            end
        end
    end
    
    if verbose
        figure;colorImageOverlay_local(abs(img),mask_new,0.5);
    end
    
end

% local function: montagesub2imgsub_local
function [x,y,z] = montagesub2imgsub_local(c,r,dim,Ncc)
    
    Ny = dim(1);
    Nx = dim(2);
    Nz = dim(3);
    
    if ~exist('Ncc','var')
        Ncc = ceil(Nz/floor(sqrt(Nz))); % number of big columns
    end
    cc  = floor((c-1)/Nx);
    rr  = floor((r-1)/Ny);
    z = 1 + cc + Ncc*rr;
    x = c - Nx*cc;
    y = r - Ny*rr;
    
    return;
    
end

% local function: colorImageOverlay_local
function h = colorImageOverlay_local(img_udr, img_ovr, alpha, opts )
    
    if ~exist('opts','var')
        opts = struct;
    end
    
    normRange = @(x)( (x-min(x(:)))/max(abs(x(:))) );
    
    if (size(img_udr,4) ~= 1)
        error('img_udr must not be a multilayer image (size(img_udr,4) must = 1)');
    end
    
    assert(all(size(img_ovr)==size(img_udr)));
    
    img_udr = repmat(img_udr,[1,1,1,3]);
    img_udr = normRange(img_udr);
    
    montagesc(img_udr,opts);
    hold on;
    h = montagesc(img_ovr,opts);
    colormap([[0,0,0];jet(256)]);
    set(h, 'alphadata', alpha);
    
    if nargout == 0
        clear h;
    end

end

