%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampling, varibale density R = 16;
% initialize: P-LORAKS
%% step 1, update coil sensitivity functions with fixed image
%% step 2, update image with fixed coil sensitivity functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
addpath(genpath(pwd));
% load data
load(fullfile('./data','data.mat'));
results_data_subj = fullfile(pwd,'results');
if ~exist(results_data_subj)
    mkdir(results_data_subj)
end
%% generate GT
opt.tempKsapce = tempKsapce;
dim_spatial = size(tempKsapce);
rGold = my_rss(Fn_k2x(tempKsapce,1:3));
tic
[~,sen_GT_SOS] = pmri_SOS_sensitivity(Fn_k2x(tempKsapce,1:3),brainMask);
toc

save(fullfile(results_data_subj,'GT'),'rGold','sen_GT_SOS');
se = strel('disk',5);
brainMask_dilate = imdilate(brainMask,se);

%% display initialization
opt.scaleError = 10;
imgInitPLORAKS = Fn_k2x(opt.recon,1:2);
tic
[img_update,~] = adapt_array_3d_no_phase_corr(imgInitPLORAKS,4);% we can remove the subject independent phase
toc

tic
[~,sensitivityUpdate] = pmri_SOS_sensitivity(imgInitPLORAKS,brainMask);
toc
sensitivityUpdate = sensitivityUpdate.*exp(-1i*angle(img_update));
% sen_ESPIRiT = sen_ESPIRiT_3D(opt.recon);
opt.scale = [];
figure;montagesc(squeeze(sensitivityUpdate(:,:,show_slice,:)),opt);
title('Sensitivity update Amp');colorbar;
colorbar;
MakeFigPretty(gcf,my_pos,[],my_FontSize)
axis off
saveFigPng(gcf,fullfile(results_data,sprintf('Sensitivity update Amp')));

opt.scale = [];
figure;montagesc(angle(squeeze(sensitivityUpdate(:,:,show_slice,:))),opt);
title('Sensitivity update Phase');colorbar;colormap('jet');
colorbar;
MakeFigPretty(gcf,my_pos,[],my_FontSize)
axis off
saveFigPng(gcf,fullfile(results_data,sprintf('Sensitivity update Phase')));

opt.scale = [0,1];
figure;montagesc(img_update(:,:,show_slice,:).*brainMask(:,:,show_slice),opt);
title('P-LORAKS');
colorbar;
MakeFigPretty(gcf,my_pos,[],my_FontSize)
axis off
saveFigPng(gcf,fullfile(results_data,sprintf('update image Amp')));


%% load subspace
fileNameSubspace = fullfile(pwd,'data',sprintf('Basis_all_T2W_complex_ref_unalignment_15_poly_remove_img_phase.mat'));
load(fileNameSubspace);
%% show first basis
temp = B(1,7).subspace;
temp = reshape(temp(:,1),dim_spatial(1:3));
figure;montagesc(temp.*brainMask);

%% iteration recon
opt.imgCSSENSE_origin =img_update;
z_cur = opt.imgCSSENSE_origin;
tol = 1e-3;
opt.recon_cg = [];
iteration_algorithm = 1;
iter_whole = 20;
l2_error = zeros(iter_whole,1);
opt.debug=1;
opt.show_slice = show_slice;
myMetricAll_CS_subspace = [];
rank_all = repmat(15,[iter_whole,1]);
opt.scaleError = 5;
reg_all = [1e-2,1e-2,1e-2,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
beta_reg = [1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6];
coil = 7;
coil_weight = max(reshape(abs(undersampledData),[],dim_spatial(4)),[],1);
opt.coil_weight = coil_weight;
opt.dim_spatial = dim_spatial;
while(1)
    % make document for iteration debug results
    z = z_cur;
    results_data = fullfile(results_data_subj,sprintf('iter%d',iteration_algorithm));
    opt.results_data = results_data;
    
    if ~isdir(results_data)
        mkdir(results_data);
    end
    

%% step 1, update coil sensitivity functions with fixed image
    % generate subspace
    my_rank = rank_all(iteration_algorithm);
    for i_coil=1:size(B,2)
        
        temp = B(1,i_coil).subspace;
        temp_basis = temp(:,1:my_rank);
        temp_constant = zeros(size(temp_basis))+0;
        temp_mask = repmat(reshape(brainMask_dilate,[],1),1,my_rank);
        B_temp(1,i_coil).subspace = temp_constant.*(~temp_mask)+temp_basis.*temp_mask;
        
    end
    
    % subsapce model
 %% in k-space
    opt.method = 'p+nn';
    load(fullfile( pwd,'data','/train_c_est_update_in_k_T15.mat'));
    opt.c_est_var = c_est_var(:);%% c_est_var is the std of c
    opt.c_est_mean = c_est_mean(:);
    opt.pro_beta = beta_reg(iteration_algorithm);
     opt.results_data = results_data;
     opt.brainMaskRef = brainMask_dilate;
     opt.dim_spatial = dim_spatial;
     New_dim = [prod(dim_spatial(1:3)),dim_spatial(4)];
     opt.debug = 1;
     opt.order = 5;
     opt.brainMaskRef = brainMask;
     opt.brainMask = brainMask;
     opt.conf_mask = brainMask;
     opt.reg_para_nn = 1e-5; 
     [opt] = step5UpdateSenModelError_subspace(B_temp,undersampledData,opt);
     opt.scale = [0,0.5];
     projection_ESPIRiT = opt.senSubspaceProjectionL2;
     projection_ESPIRiT = senNormRange(projection_ESPIRiT).*repmat(brainMask_dilate,[1,1,1,dim_spatial(4)]);
     projection_ESPIRiT(isnan(projection_ESPIRiT)) = eps;
     projection_ESPIRiT(projection_ESPIRiT==0) = eps;
     figure;montagesc(projection_ESPIRiT(:,:,:,7),opt);
       
    %% check the results of step 1
     opt.scale = [0,0.5];
    figure;montagesc(abs(squeeze(projection_ESPIRiT(:,:,show_slice,:))),opt);
    title('Subspace model magnitude');
    opt.scale = [];
    figure;montagesc(angle(squeeze(projection_ESPIRiT(:,:,show_slice,:))),opt);colormap('jet');
    title('Subspace model phase');
    opt.ErrorScale = 5;
    opt.scale = [0,0.5];
    figure;montagesc(squeeze(sen_GT_SOS(:,:,show_slice,:)-projection_ESPIRiT(:,:,show_slice,:)).*opt.ErrorScale,opt);colorbar;
    title('GT Error');
    %% step 2, update image with fixed coil sensitivity functions
    % gen sparseSENSE recon
     opt.ssim_location = [10,35]; % ssim text location
     opt.groundTruth = rGold;
     opt.undersampleMask = undersampleMask;
     opt.results_data = fullfile(results_data,sprintf('Pro-Subspace-SENSET15-CS'));
     opt.sen_name = 'ESPIRiT proposed';
     opt.i_subj = 1;
     opt.debug=1;
     opt.R_sen = 1;
     opt.R = 1;
     opt.iteration_algorithm = 1;
     opt.reg = reg_all(iteration_algorithm);
     opt.sen_methd = 1;
     opt.scaleError = 5;
     se = strel('disk',0);
     opt.brainMask = imdilate(brainMask,se);
     opt.myMetricAll_CS = myMetricAll_CS_subspace;
     opt.show_slice = show_slice;
     opt.scale = [];
     opt.temp_groundTruth_sen = sen_GT_SOS;
     opt = recon_with_sen_sparseSENSE(projection_ESPIRiT,double(undersampledData),opt);
     myMetricAll_CS_subspace = opt.myMetricAll_CS;
     figure;montagesc(abs(opt.imgCSSENSE_origin(:,:,1)))
     figure;montagesc(angle(opt.imgCSSENSE_origin(:,:,1)));colormap('jet');
   
         z_cur = opt.imgCSSENSE_origin;
    t = (norm(z_cur(:)-z(:))/norm(z(:)));
        % check for convergence
    if t < tol
        disp('Convergence tolerance met: change in solution is small');
        break;
    end
    % display the status
    if ~rem(iteration_algorithm,1)
        disp(['iter ' num2str(iteration_algorithm) ', relative change in solution: ' num2str(t)]);
    end
   

           % sythesize k-space
    SENSE3DOp = SENSE3D(projection_ESPIRiT);
    dkUpdate = Fn_x2k(SENSE3DOp*opt.imgCSSENSE_origin,1:3);
    % update k-space
    dkUpdate = dkUpdate.*(1-undersampleMaskC)+undersampledData;
    opt.scale = [0,0.5];
    figure;montagesc(my_rss(Fn_k2x(dkUpdate,1:3)),opt);colorbar;
    l2_error(iteration_algorithm) = norm(dkUpdate(:)-tempKsapce(:),2);
    save(fullfile(opt.results_data,'l2_error_kspace.mat'),'l2_error')
    tic
[img_update,~] = adapt_array_3d_no_phase_corr(Fn_k2x(dkUpdate,1:3),4);% we can remove the subject independent phase
toc
% figure;montagesc(img_update)
%% update image with data consitency
opt.imgCSSENSE_origin = img_update;
 tic
    [~,sensitivityUpdate] = pmri_SOS_sensitivity(Fn_k2x(dkUpdate,1:3),brainMask);
    toc
    sensitivityUpdate = sensitivityUpdate.*exp(-1i*angle(img_update));
    save(fullfile(opt.results_data,'sen_show.mat'),'sensitivityUpdate');
  
    %% iterative
    fprintf('finish iteration %d ',iteration_algorithm)
    iteration_algorithm = iteration_algorithm+1;
    
    %% stop or not
    if iteration_algorithm>iter_whole
        disp('Meet the max iteration:',num2str(iteration_algorithm));
        
        break;
    end
end


my_metric_plot = struct2table(cell2mat(myMetricAll_CS_subspace));
reg_log = 1:iteration_algorithm-1;
figure;
subplot(2,2,1);plot(1:iteration_algorithm-1,my_metric_plot.ssim);
hold on;plot(reg_log,my_metric_plot.ssim,'k*');
hold on;plot(reg_log,my_metric_plot.mssim,'r*')
xlabel('log(\mu)');
ylabel('SSIM');
MakeFigPretty(gcf)

subplot(2,2,2);
plot(reg_log,abs(my_metric_plot.psnr));
xlabel('log(\mu)');
ylabel('PSNR');
MakeFigPretty(gcf,[],5)

subplot(2,2,3);
plot(reg_log,abs(my_metric_plot.l2error));
xlabel('log(\mu)');
ylabel('L2Error');
MakeFigPretty(gcf,[],5)


subplot(2,2,4);
plot(reg_log,abs(my_metric_plot.nmse));
xlabel('log(\mu)');
ylabel('NMSE');
MakeFigPretty(gcf)
saveFigPng(gcf,fullfile(results_data,sprintf('metric vs.mu %s iter %d',opt.i_subj,opt.iteration_algorithm)));

