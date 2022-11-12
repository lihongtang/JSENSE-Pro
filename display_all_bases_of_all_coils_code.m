%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating all the basis functions for all the coils 
%based on the 15 training data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
fullpath = mfilename('fullpath'); 
[path,name]=fileparts(fullpath);
addpath(genpath(path));
cd(path);
% load data
load(fullfile('./data','bases.mat'));
results_data_subj = fullfile(path,'results');
if ~exist(results_data_subj)
    mkdir(results_data_subj)
end
my_pos = [0.25,0.25,8.5,4.2];
for coil = 1:20
slice = 1;
%% mag
opt.scale = [0,0.5e-3];
figure_name = sprintf('mag_coil%d_slice%d.png',coil,slice);
figure;montagesc(abs(bases(:,:,:,coil)),opt);
MakeFigPretty(gcf,my_pos);
axis off;
im = gcf;
print(im,'-dpng','-r300',fullfile(results_data_subj,figure_name));
close all;

%% phase
opt.scale = [-pi,pi];
figure_name = sprintf('phase_coil%d_slice%d.png',coil,slice);
figure;montagesc(angle(bases(:,:,:,coil)),opt);colormap('jet');
% figure;montagesc(reshape(show_sen.*brainMask(:,:,slice),[dim_spatial(1),dim_spatial(2)*15]),opt);colormap('jet')
MakeFigPretty(gcf,my_pos);
axis off;
im = gcf;
print(im,'-dpng','-r300',fullfile(results_data_subj,figure_name));
close all;
sprintf("finish coil %d",coil)
end
