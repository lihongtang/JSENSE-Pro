function MakeFigPretty( fh, pos, lw,FontSize )
%MAKEFIGPRETTY Makes a figure pretty
%   fh      = figure handle
%   pos     = figure position vector. pos = [ loc_x loc_y width height ] in
%             inches.
%   lw      = width of line objects in the figure
%--------------------------------------------------------------------------

    if nargin < 2
        pos = [.25 .25 10 8];
    end

    axh = findobj(fh,'type','axes','-not','tag','legend');
    
    if ~isempty(pos)
        set(fh,'Units','Inches', 'Position', pos);
    end
     if ~exist('FontSize')
        FontSize=30;
    end
    for n=length(axh):-1:1 % oldest first
        title_h = get(axh(n),'Title');
%                 sgtitle_h = get(axh(n),'sgTitle');

        xlabel_h = get(axh(n),'Xlabel');
        ylabel_h = get(axh(n),'Ylabel');
        
        set([title_h,xlabel_h,ylabel_h],'FontName','Calibri');
        set(title_h,'FontSize',FontSize,'FontWeight','Bold');
%                 set([sgtitle_h,xlabel_h,ylabel_h],'FontName','Calibri');
%         set(sgtitle_h,'FontSize',FontSize,'FontWeight','Bold');
        set(ylabel_h,'FontSize',FontSize);
        set(xlabel_h,'FontSize',FontSize);
        set(axh(n),'FontSize',FontSize,'LineWidth',3);
    end
    
    if ~exist('lw')||isempty(lw)
        lw = 3;
        
    end
    lines = findall(fh,'type','line');
    set(lines(:),'LineWidth',lw);
     set(fh,'color','w');
    set(gca,'Fontname','Calibri');
    set(gca,'color','none');


end

