function opt = recon_with_sen_sparseSENSE(usedSensitivity,temp_undersample,opt,large_iter)

results_data = opt.results_data;
if ~isdir(results_data)
    mkdir(results_data);
end

if ~exist('large_iter','var')||isempty(large_iter)
    large_iter = 1;
end

scaleError = opt.scaleError;
brainMask = opt.brainMask;
dim_spatial = opt.dim_spatial;
% coil_weight = opt.coil_weight;
myMetricAll_CS = opt.myMetricAll_CS;
groundTruth = opt.groundTruth;
show_slice = opt.show_slice;

if ~isfield(opt,'scale')
     opt.scale = [0,0.8];
end
if ~isfield(opt,'scale_sen')
     opt.scale_sen = [0,0.5];
end

% sensitivityfullsampleZF = opt.sensitivityfullsampleZF;
data_R = round(1./(nnz(temp_undersample)./size(temp_undersample(:),2)));
if opt.sen_methd==1
    usedSensitivity = usedSensitivity;
    
    
else
    
    if opt.sen_methd==2
        constant_matrix = 1./sqrt(dim_spatial(4))*ones(size(usedSensitivity)).*~repmat(brainMask,[1,1,1,dim_spatial(4)]);
        usedSensitivity = usedSensitivity.*repmat(brainMask,[1,1,1,dim_spatial(4)])+constant_matrix;
    else
        
        constant_matrix = opt.sensitivityCS.*~repmat(brainMask,[1,1,1,dim_spatial(4)]);
        usedSensitivity = usedSensitivity.*repmat(brainMask,[1,1,1,dim_spatial(4)])+constant_matrix;
    end
    
end
if opt.debug==1
% sen
opt_temp = opt;
opt_temp.scale = opt.scale_sen;
show_img = [];
show_img.sen = squeeze(usedSensitivity(:,:,show_slice,:));
figure;montagesc(show_img.sen,opt_temp);
title(sprintf('Sensitivity (iter %d)\n (R = %d)',opt.iteration_algorithm,opt.R_sen));
colorbar
MakeFigPretty(gcf)
axis off
saveFigPng(gcf,fullfile(results_data,sprintf('sensitivity %s (R = %d) iter = %d ',opt.sen_name,opt.R_sen,opt.iteration_algorithm)));

% sen error

opt_temp = opt;
opt_temp.scale = opt.scale_sen;
temp_error = (opt.temp_groundTruth_sen-usedSensitivity).*scaleError;
show_img.temp_error = squeeze(temp_error(:,:,show_slice,:)).*brainMask(:,:,show_slice);
figure;montagesc(show_img.temp_error,opt_temp);
title(sprintf('Error (x %d) \n (R = %d)',scaleError,opt.R_sen));
colorbar
MakeFigPretty(gcf)
axis off
saveFigPng(gcf,fullfile(results_data,sprintf('sensitivity error %s (R = %d) iter = %d',opt.sen_name,opt.R_sen,opt.iteration_algorithm)));
end
opt_temp = [];

opt.usedSensitivity = usedSensitivity;
usedKsapce = double(temp_undersample);
[~,opt] = step6_sparseSENSE_SSIM_for_optimal_regularization(usedKsapce,usedSensitivity,groundTruth,opt.undersampleMask,opt);
% DC for each recon
opt.imgCSSENSE_origin = opt.imgCSSENSE;
opt.imgCSSENSE = double(abs(opt.imgCSSENSE));
% [~,opt.multicoil_img,~] = dataReplacement(opt.imgCSSENSE_origin,temp_undersample,usedSensitivity,coil_weight);
close all;
% figure(1);montagesc(rot90(my_rss(opt.multicoil_img(:,:,show_slice,:))))
% figure(2);montagesc(rot90(abs(opt.imgCSSENSE_origin(:,:,show_slice,:))))
% opt.imgCSSENSE = my_rss(opt.multicoil_img);
%% recon
myMetricInitCS_sen = myAssessmentImageQuality(opt.imgCSSENSE(:,:,show_slice),groundTruth(:,:,show_slice),opt,brainMask(:,:,show_slice));
myMetricAll_CS{opt.iteration_algorithm} = myMetricInitCS_sen;
opt.myMetricAll_CS = myMetricAll_CS;
% opt.scale = [0,0.8];
if opt.debug==1

figure;montagesc(squeeze(abs(opt.imgCSSENSE(:,:,show_slice))),opt);
title(sprintf('CS-SENSE (iter %d)\n (R = %d)',opt.iteration_algorithm,opt.R));
colorbar;
% show the metrics
text( opt.ssim_location(1), opt.ssim_location(2),...
    sprintf('\\color{red}\\bfSSIM = %.3f\nPSNR = %.2f\nNMSE = %.3f',...
    myMetricInitCS_sen.mssim,myMetricInitCS_sen.psnr,myMetricInitCS_sen.nmse),...
    'FontSize',17 );
MakeFigPretty(gcf)
axis off
saveFigPng(gcf,fullfile(results_data,sprintf('CS-SENSE %s (R = %d) reg = %d iter = %d SSIM = %d',opt.sen_name,opt.R,opt.opt_reg,opt.iteration_algorithm,opt.myMetric.mssim )));
end
show_img.recon = squeeze(abs(opt.imgCSSENSE(:,:,show_slice)));
%% Error
% opt.scale = [];
temp_error = (opt.imgCSSENSE(:,:,show_slice)-groundTruth(:,:,show_slice)).*scaleError;
show_img.recon_error = temp_error;
if opt.debug==1
figure;montagesc(squeeze(abs(temp_error)),opt);
title(sprintf('Error (x %d)\n (R = %d)',scaleError,opt.R));
colorbar
text( opt.ssim_location(1), opt.ssim_location(2),...
    sprintf('\\color{red}\\bfError x %d',...
    scaleError),...
    'FontSize',17 );

axis off
MakeFigPretty(gcf)
saveFigPng(gcf,fullfile(results_data,sprintf('CS-SENSE Error %s(R = %d) reg = %d iter = %d SSIM = %d',opt.sen_name,opt.R,opt.opt_reg,opt.iteration_algorithm,myMetricInitCS_sen.mssim)));
end
% test_recon = opt.imgCSSENSE(:,:,show_slice);
% test_gt = groundTruth(:,:,show_slice);
% test = normRange(cat(3,test_recon,test_gt));
% if opt.debug==1
% 
% figure;montagesc(rot90(squeeze(normRange(abs(test(:,:,2))))),opt);
% title(sprintf('CS-SENSE (iter %d)\n (R = %d)',opt.iteration_algorithm,opt.R));
% colorbar;
% % show the metrics
% text( opt.ssim_location(1), opt.ssim_location(2),...
%     sprintf('\\color{red}\\bfSSIM = %.3f\nPSNR = %.2f\nNMSE = %.3f',...
%     myMetricInitCS_sen.mssim,myMetricInitCS_sen.psnr,myMetricInitCS_sen.nmse),...
%     'FontSize',17 );
% MakeFigPretty(gcf)
% axis off
% saveFigPng(gcf,fullfile(results_data,sprintf('CS-SENSE %s (R = %d) reg = %d iter = %d SSIM = %d',opt.sen_name,opt.R,opt.opt_reg,opt.iteration_algorithm,opt.myMetric.mssim )));
% end
% %% no normrage
% %% Error
% opt.scale = [0,0.5];
% temp_error = abs(abs(test(:,:,1))-abs(test(:,:,2))).*scaleError;
% show_img.recon_error = temp_error;
% if opt.debug==1
% figure;montagesc(rot90(squeeze(abs(temp_error))),opt);
% title(sprintf('Error (x %d)\n (R = %d)',scaleError,opt.R));
% colorbar
% text( opt.ssim_location(1), opt.ssim_location(2),...
%     sprintf('\\color{red}\\bfError x %d',...
%     scaleError),...
%     'FontSize',17 );
% 
% axis off
% MakeFigPretty(gcf)
% saveFigPng(gcf,fullfile(results_data,sprintf('CS-SENSE Error %s(R = %d) reg = %d iter = %d SSIM = %d',opt.sen_name,opt.R,opt.opt_reg,opt.iteration_algorithm,myMetricInitCS_sen.mssim)));
% end

save(fullfile(results_data,sprintf('Show_data_sen_recon.mat')),'myMetricInitCS_sen','show_img');
opt_save.imgCSSENSE = opt.imgCSSENSE_origin;
opt_save.usedSensitivity = opt.usedSensitivity; 
opt_save.temp_groundTruth_sen = opt.temp_groundTruth_sen; 
if opt.iteration_algorithm >= large_iter
save(fullfile(results_data,sprintf('results.mat')),'opt_save','-v7.3');
end
close all;
opt.show_img = show_img;
