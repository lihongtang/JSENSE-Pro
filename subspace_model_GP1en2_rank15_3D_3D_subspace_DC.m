%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampling, varibale density R = 4;
% initialize: P-LORAKS
%% step 1, update coil sensitivity functions with fixed image
%% step 2, update image with fixed coil sensitivity functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
fullpath = mfilename('fullpath'); 
[path,name]=fileparts(fullpath);
addpath(genpath(path));
cd(path);
% load data
load(fullfile('./data','data.mat'));
load(fullfile('./data','GT.mat'));
results_data_subj = fullfile(path,'results');
if ~exist(results_data_subj)
    mkdir(results_data_subj)
end
show_slice = 1;

coil = 7;
%% display initialization
opt.scaleError = 10;
imgInitPLORAKS = Fn_k2x(opt.recon,1:2);
tic
[img_update,~] = adapt_array_3d_no_phase_corr_size1(squeeze(imgInitPLORAKS));% we can remove the subject independent phase
toc
% figure;montagesc(img_update)
opt.imgCSSENSE_origin = img_update;
toc
% opt.scale = [0,0.5];
% figure;montagesc(abs(img_update),opt);
% MakeFigPretty(gcf);
% axis off;
% im = gcf;
% print(im,'-dpng','-r300',fullfile(results_data,sprintf('SENSE_combine.png')));

% tic
% [img_update3D,~] = adapt_array_3d_no_phase_corr(imgInitPLORAKS,1);% we can remove the subject independent phase
% toc

%% load subspace
dim_spatial = size(imgInitPLORAKS); 
fileNameSubspace = fullfile(path,'data',sprintf('Basis.mat'));
load(fileNameSubspace);
T = 15;
for i_coil=1:size(B,2)
    temp = B(1,i_coil).subspace;
    N_rank = size(temp,2);
    temp = reshape(temp,[dim_spatial(1:3),N_rank]);
    temp = rot90(temp);
    B(1,i_coil).subspace = reshape(temp,[],N_rank);
end
%% iteration recon
z_cur = opt.imgCSSENSE_origin;
tol = 1e-3;


opt.recon_cg = [];
%%
iteration_algorithm = 1;
iter_whole = 100;
l2_error = zeros(iter_whole,1);

myMetricAll_CS = [];
myMetricAll_CS_poly = [];
myMetricAll_CS_poly_subspace =[];
myMetricAll_CS_poly_phase_conserve = [];
myMetricAll_CS_poly_subspace_phase_conserve =[];
myMetricAll_CS_poly_alter = [];
opt.debug=1;
opt.show_slice = show_slice;
myMetricAll_CS_subspace = [];
show_coil = 7;
rank_all = repmat(15,[iter_whole,1]);
opt.scaleError = 5;
%% investigate novel parts imapact
% sparsity constraint (lambda)
reg_all = [1e-3];
% Gaussian prior of coefficients (beta1)
beta_reg1 = [1e-2]
% nuclear-norm penalty
beta_reg2 = [1e-17]

coil_weight = max(reshape(abs(undersampledData),[],dim_spatial(4)),[],1);
opt.coil_weight = coil_weight;
opt.dim_spatial = dim_spatial;
while(1)
    if iteration_algorithm>iter_whole
        disp('Meet the max iteration:',num2str(iteration_algorithm)); 
        disp(sprintf('Stop Iteration:%d',iteration_algorithm-1)); 
        break;
    end
    

    %% subspace projection
    % generate subspace
    my_rank = rank_all(iteration_algorithm);
    for i_coil=1:size(B,2)
        
        temp = B(1,i_coil).subspace;
        temp_basis = temp(:,1:my_rank);
        temp_constant = zeros(size(temp_basis))+0;
        temp_mask = repmat(reshape(brainMask_dilate,[],1),1,my_rank);
        temp = reshape(temp_constant.*(~temp_mask)+temp_basis.*temp_mask,[dim_spatial(1:3),my_rank]);
        B_temp(1,i_coil).subspace = reshape(temp,[],my_rank);
        
    end
    z = z_cur;
    results_data = fullfile(results_data_subj,sprintf('SparseSENSET%d-poly-rank%d-GP1en2-NN-3D-3DSub-AdaptSize1-DC/iter%d',T,my_rank,iteration_algorithm));
    
    if ~isdir(results_data)
        mkdir(results_data);
    end
    % subsapce model
 %% in k-space
    opt.method = 'p';
    fileNameSubspace = fullfile(path,'data',sprintf('train_c_est.mat'));
load(fileNameSubspace);
   opt.c_est_var = c_est_var(:);%% c_est_var is the std of c
    opt.c_est_mean = c_est_mean(:);
    opt.pro_beta = beta_reg1(1);
    opt.reg_para_nn = beta_reg2(1); 
    opt.b = 1e-13;
     opt.results_data = results_data;
     opt.brainMaskRef = brainMask_dilate;
     opt.dim_spatial = dim_spatial;
     New_dim = [prod(dim_spatial(1:3)),dim_spatial(4)];
     opt.debug = 1;
     opt.order = 5;
     opt.brainMaskRef = brainMask;
     opt.brainMask = brainMask;
     opt.conf_mask = brainMask;
     [opt] = step5UpdateSenModelError_subspace(B_temp,undersampledData,opt);
     opt.scale = [0,0.5];
     projection_ESPIRiT = opt.senSubspaceProjectionL2;
     projection_ESPIRiT = projection_ESPIRiT.*repmat(brainMask_dilate,[1,1,1,dim_spatial(4)]);
     projection_ESPIRiT(isnan(projection_ESPIRiT)) = eps;
     projection_ESPIRiT(projection_ESPIRiT==0) = eps;
  
    %% gen SENSE recon
     opt.ssim_location = [10,35]; % ssim text location
     opt.groundTruth = rGold;
     opt.undersampleMask = undersampleMask;
     opt.results_data = results_data;
     opt.sen_name = 'ESPIRiT proposed';
     opt.i_subj = 1;
     opt.debug=1;
     opt.R_sen = 1;
     opt.R = 1;
     opt.iteration_algorithm = 1;
     opt.reg = reg_all(1);
     opt.sen_methd = 1;
     opt.scaleError = 5;
     se = strel('disk',0);
     opt.brainMask = imdilate(brainMask,se);
     opt.myMetricAll_CS = myMetricAll_CS_subspace;
     opt.show_slice = show_slice;
     opt.scale = [];
     opt.temp_groundTruth_sen = sen_GT_SOS;
     opt = recon_with_sen_sparseSENSE(projection_ESPIRiT,undersampledData,opt);
     myMetricAll_CS_subspace = opt.myMetricAll_CS;
         z_cur = opt.imgCSSENSE_origin;
           % sythesize k-space
    SENSE3DOp = SENSE3D(projection_ESPIRiT);
    dkUpdate = Fn_x2k(SENSE3DOp*opt.imgCSSENSE_origin,1:3);
    % update k-space
    dkUpdate = dkUpdate.*(1-undersampleMaskC)+undersampledData;
    l2_error(iteration_algorithm) = norm(dkUpdate(:)-tempKsapce(:),2);
    save(fullfile(opt.results_data,'l2_error_kspace.mat'),'l2_error')
    recon_cg = opt.recon_cg;
    save(fullfile(opt.results_data,'recon_cg_sen.mat'),'recon_cg');
    tic
    [img_update,~] = adapt_array_3d_no_phase_corr_size1(squeeze(Fn_k2x(dkUpdate,1:3)));% we can remove the subject independent phase
    toc
    opt.imgCSSENSE_origin = img_update;
    t = (norm(z_cur(:)-z(:))/norm(z(:)));
    % out iteration error
    fprintf('OutIteration: %03d    Rel Error: %10.4e\n',iteration_algorithm,t);
    %% stop or not
    % check for convergence
    if t < tol
    disp('Convergence tolerance met: change in solution is small');
    disp(sprintf('Stop Iteration:%d',iteration_algorithm)); 
    break;
    end
    iteration_algorithm = iteration_algorithm+1;
end

