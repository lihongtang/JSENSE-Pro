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

