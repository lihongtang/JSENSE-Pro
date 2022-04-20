function h = montagesc( I, opts )
%MONTAGESC A combination of imagesc and montage
%   Creates a single image montage of the stacks in I.
%
% ---- INPUTS ----
% I                 : a 3D matrix of the image stacks
% opts.showgrid     : if true (default) gridlines are shown between stacks in
%                     the figure
% opts.gridColor    : if given, sets the color of the grid (default white)
% opts.gridStyle    : if given, sets the style of the grid (default '-')
% opts.size         : grid size (default null)
% opts.xframes      : flag to place 'X' on extra frames (default true)
%
%   Modified by Yibo Zhao @ UIUC, 2018/07/31:
%       When the image is complex, warns the user and converts it to abs
%   Modified by Lihong tang @ SJTU, 2021/06/10:
%       adding scale option
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------

    if nargin < 2
        opts = struct;
    end
    if isfield(opts,'showgrid')
        showgrid = opts.showgrid;
    else
        showgrid = true;
    end
    if isfield(opts,'gridColor')
        gridColor = opts.gridColor;
    else
        gridColor = 'w';
    end
    if isfield(opts,'gridStyle')
        gridStyle = opts.gridStyle;
    else
        gridStyle = '-';
    end
    if isfield(opts,'xframes')
        xframes = opts.xframes;
    else
        xframes = true;
    end
    if isfield(opts,'scale')
        scale = opts.scale;
    else
        scale = [];
    end

    % Dimensions of subplot grid
    [M, N, P, Q] = size(I);
    if ~( (Q == 1) || (Q == 3) )
        error('Input image must be mono color or true color.');
    end
    
    I = permute(I,[1,2,4,3]);
    
    if isfield(opts,'size')
        a = opts.size(1);
        b = opts.size(2);
    else
        a = floor(sqrt(P));
        b = ceil(P/a);
    end
    extra_frames = a*b > P;
    if showgrid && extra_frames
        bad_frame_inds = zeros(2,2,P-a*b);
    else
        bad_frame_inds = [];
    end

    
    % Make image
    img = -zeros([a*M,b*N,Q]);
    p = 0;
    for ii = 0:a-1
        for jj = 0:b-1
            p = p+1;
            if p > P
                if showgrid
                    bad_frame_inds(:,:,a*b-p+1) = [ii*M + 1, ii*M + M;...
                                                   jj*N + 1, jj*N + N];
                else
                    break;
                end
            else
                img( ii*M + (1:M), jj*N + (1:N), 1:Q) = I(:,:,:,p);
            end
        end
    end
    
    % check
    if ~isreal(img)
        % @@ Yibo: this warning is necessary, e.g.
        % img_truth = phantom(64)-0.5;                                  % another example can be a field map
        % kspace_noisy = F_x2k(img_truth) + eps*rand(size(img_truth));  % just a little noise in k-space
        % figure;montagesc(F_k2x(kspace_noisy));                        % the figure is totally wrong!
        warning('Image data can not be complex! Converted to magnitude image.');
        
        img = abs(img);
    end
    
    % Make figure
    if ~isempty(scale)
    if nargout > 0
        h = imagesc(img,scale);
    else
        imagesc(img,scale);
    end
    else
            if nargout > 0
        h = imagesc(img);
    else
        imagesc(img);
    end
    end
    if Q == 1
        colormap(gray(256));
    end
    ah = gca;
    grid_opts.color = gridColor;
    grid_opts.style = gridStyle;
    grid_opts.linewidth = 2;
    if showgrid
        addgrid2img(ah,[a*M,b*N],[a,b],grid_opts);
        if xframes
            if extra_frames
                hold(ah,'all');
                for n=1:size(bad_frame_inds,3)
                    line(bad_frame_inds(2,[1,2],n) + [1,-1],...
                         bad_frame_inds(1,[1,2],n) + [1,-1],...
                         'linestyle',gridStyle,...
                         'color',gridColor,...
                         'linewidth',grid_opts.linewidth);
                    line(bad_frame_inds(2,[2,1],n) - [1,-1],...
                         bad_frame_inds(1,[1,2],n) + [1,-1],...
                         'linestyle',gridStyle,...
                         'color',gridColor,...
                         'linewidth',grid_opts.linewidth);
                end
                hold(ah,'off');
            end
        end
    end
    axes(ah);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addgrid2img( ah, img_dim, grid_dim, opts )
%ADDGRID2IMG Overlays gridlines to an image
%
% ---- INPUTS ----
% ah            :   axis handle to image axis
% img_dim       :   dimensions (size) of the image to overlay the grid on
% grid_dim      :   dimensions of the grid
% opts.color    :   color of grid lines
% opts.style    :   style of grid lines
% opts.lwidth   :   line width of grid lines
%-------------------------------------------------------------------------------

    if nargin < 4
        opts = struct;
    end
    if isfield(opts,'color')
        lcolor = opts.color;
    else
        lcolor = 'g';
    end
    if isfield(opts,'style')
        lstyle = opts.style;
    else
        lstyle = '-';
    end
    if isfield(opts,'lwidth')
        lwidth = opts.lwidth;
    else
        lwidth = 2;
    end
    
    hold(ah,'all');
    dyx = img_dim./(grid_dim);
    for py=1:grid_dim(1)-1
        line([0,img_dim(2)]+.5,[dyx(1)*py, dyx(1)*py]+.5,'linestyle',lstyle,'color',lcolor,'linewidth',lwidth)
    end
    for px=1:grid_dim(2)-1
        line([dyx(2)*px, dyx(2)*px]+.5,[0,img_dim(1)]+.5,'linestyle',lstyle,'color',lcolor,'linewidth',lwidth)
    end
end






