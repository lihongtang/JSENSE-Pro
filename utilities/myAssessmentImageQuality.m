function my_metric = myAssessmentImageQuality(V_recon,V_ground,opt,mask4SSIM)
%% ==================================================================
%CoilSensitivityMapsRecon3D
% reference: fastMRI: an open dataset and benchmarks for accelerated MRI
% ===================================================================
%   Author: lihong @SJTU
%   Created: 2021-06-02
%
%   [INPUTS]
%   ---- Required ----
%   V_recon                  reconstruction results
%   V_ground                 ground truth
%   ---- Optional ----
%   opt:
%     - method             metrics method

%   [OUTPUTS]
%   Metric_value
%
%   Change log:

%--------------------------------------------------------------------------

%% ------ parse the input ------

if ~exist('mask4SSIM','var')
    mask4SSIM = ones(size(V_recon));
end
% use abs to calclulate structure SSIM and mean SSIM
%% ssim
temp = zeros(size(V_ground,3),1);
temp2 = zeros(size(V_ground,3),1);

for i_slice = 1:size(V_ground,3)   
%     temp_ground = normRange(squeeze(abs(V_ground(:,:,i_slice).*mask4SSIM(:,:,i_slice)))).*255;
%     temp_recon = normRange(squeeze(abs(V_recon(:,:,i_slice).*mask4SSIM(:,:,i_slice)))).*255; 
%     temp_ground = normRange(squeeze(abs(V_ground(:,:,i_slice).*mask4SSIM(:,:,i_slice))));
%     temp_recon = normRange(squeeze(abs(V_recon(:,:,i_slice).*mask4SSIM(:,:,i_slice)))); 
    temp_ground = squeeze(abs(V_ground(:,:,i_slice).*mask4SSIM(:,:,i_slice)));
    temp_recon = squeeze(abs(V_recon(:,:,i_slice).*mask4SSIM(:,:,i_slice)));  
%     [temp(i_slice),~,~,temp2(i_slice),~] = ssim_para(temp_recon,temp_ground,[0.01,0.03],ones(7),max(temp_ground(:)));
[temp(i_slice), ~] = ssim_index(temp_recon, temp_ground, [0.01,0.03], ones(7), max(temp_ground(:)));
temp2(i_slice) = temp(i_slice);
    %  [test1,test2] = ssim(temp_recon,temp_ground);
%     figure;imshow([temp_recon temp_ground],[])
 close all
end

my_metric.ssim = mean(temp);
my_metric.mssim = mean(temp2);
% use the complex to obtain other metrics
V_recon = V_recon(:);
V_ground = V_ground(:);

my_metric.nmse = nmse(V_recon,V_ground);

my_metric.psnr = psnr(V_recon,V_ground);

my_metric.l1error = l1error(V_recon,V_ground);
my_metric.l2error = norm(V_recon-V_ground,2);
end

function Metric_value = nmse(V_recon,V_ground)

Metric_value = norm(V_recon-V_ground,2).^2./(norm(V_ground,2).^2);
end

function Metric_value = psnr(V_recon,V_ground)
% use abs for calculate PSNR
V_recon = abs(V_recon);
V_ground = abs(V_ground);

Metric_value = 10*log10(max(V_ground).^2./(mse(V_recon,V_ground)));
end

function Metric_value = l1error(V_recon,V_ground)

Metric_value = norm(V_recon-V_ground,1);
end

