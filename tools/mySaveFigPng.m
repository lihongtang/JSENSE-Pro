function mySaveFigPng(fh,save_fig_dir,filename,opt)
% MYSAVEFIGPNG save a figure in .fig and .png formats
% [Input]
%   fh           : figure handle to save
%   save_fig_dir : directory to save figures
%   filename     : saving file name
%
% Modified by Yibo Zhao, 2018/06/26: don't use frame
% Modified by Yibo Zhao, 2018/07/17: an opt for frame
% Modified by Yibo Zhao, 2018/09/03: use fullfile function
%
% See also SAVEFIGPNG

if isempty(save_fig_dir)
    save_fig_dir = fullfile(filesep,'home','yibo','Temp_transfer',filesep);
end

if ~exist('opt','var') || isempty(opt)
    opt = struct('frame',false);
end

if ~(exist(save_fig_dir,'dir')==7)
%     disp('Warning: saving figure directory not found; creating it now.');
%     warning('Saving figure directory not found; creating it now.');
    mkdir(save_fig_dir);
end

saveFigPng(fh,fullfile(save_fig_dir,filename),opt);

% Future: arbitrary number of input strings

