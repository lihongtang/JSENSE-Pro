function saveFigPng(fh,name,opt)
% SAVEFIGPNG Saves the figure with handle 'fh' as a .fig file and an .png file 
%
% INPUTS
% fh    : figure handle
% name  : filename without extension
%
% Last edited 04/2018 Bryan A Clifford
%-------------------------------------------------------------------------------

    if exist('opt','var')
        if ~isfield(opt,'loose')
            opt.loose = false;
        end
        if ~isfield(opt,'do_save')
            opt.do_save = true;
        end
        if ~isfield(opt,'do_save_fig')
            opt.do_save_fig = true;
        end
        if ~isfield(opt,'frame')
            opt.frame = true;
        end
    else
        opt.loose = false;
        opt.frame = true;
        opt.do_save = true;
        opt.do_save_fig = true;
    end

    if opt.do_save

        % save the figure
        if opt.do_save_fig
            hgsave(fh, [name '.fig']);
        end

        % print the figure as a png file
        set(fh, 'PaperPositionMode', 'auto');
        if opt.frame
            imwrite(getfield(getframe(fh),'cdata'), [name '.png']);
        else
            if opt.loose
                print(fh, '-dpng', '-r300','-painters', '-loose', [name '.png']);
            else
                print(fh, '-dpng', '-r300','-painters', [name '.png']);
            end
        end
    end

end
