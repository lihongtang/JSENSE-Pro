function opt = step5UpdateSenModelError_subspace(B,temp_undersample,opt)%% projection onto subsapce L2
% load(fileNameSubspace);

if ~isfield(opt,'method')
    method = 'L2';
else
    method = opt.method;
end
if ~isfield(opt,'debug')
    debug = 1;
else
    debug = opt.debug;
end
if ~isfield(opt,'results_data')
    error('Please set results_data path')
else
    results_data = opt.results_data;
end

if ~isfield(opt,'brainMaskRef')
    error('Please set reference brain mask');
else
    brainMaskRef = opt.brainMaskRef;
    figure;montagesc(brainMaskRef);
    
end

if ~isfield(opt,'dim_spatial')
    error('Please set dim_spatial')
else
    dim_spatial = opt.dim_spatial;
end
img = opt.imgCSSENSE_origin;
switch(lower(method))
    case 'l2'
        %         clear classes
        undersampleMask = logical(abs(temp_undersample));
        % parameters for reconstruction
        % figure;montagesc(brainMaskRef);
        % brainMaskRef_trust = repmat(reshape(ones(size(brainMaskRef)),[],1),[1,size(testSen,2)]);
        param.E = subspace_model_error_mat(img,B,dim_spatial,undersampleMask);
        
        
        % L1Weight was chosen in an empirical way, but it seems to work very well
        % if the images are normalized from 0-1
        % param.W = TV_Temp();param.L1Weight=0.1;
        param.W = 1;param.L1Weight=0;
        param.TV = TVOP();param.TVWeight=0;
        param.y = temp_undersample;
        param.nite = 8;
        param.display=1;
        
        % linear reconstruction
        if ~isfield(opt,'recon_cg')||isempty(opt.recon_cg)
            
            recon_coe = param.E'*param.y;
        else
            recon_coe = opt.recon_cg;
        end
        % %% check the operator
        % checkAdjoint(@(x)param.E*x,...
        %     @(x)param.E'*x,...
        %     [size(recon_coe,1),1],1e-10,1);
        
        
        AHA = @(x)param.E'*(param.E*x);
        AHb = param.E'*param.y;
        maxit = 200;
        tol = 1e-6;
        recon_cg = mypcg(AHA,AHb,tol,maxit,[],[],recon_coe);
        
        % k-t SPARSE-SENSE reconstruction
        fprintf('\n  sum of subsapce projection \n')
        
        %% combination the coefficients and basis
        senSubspaceProjectionL2 =  zeros(size(B(1).subspace,1),size(B,2));
        for i_coil = 1:size(B,2)
            rank(i_coil) = size(B(i_coil).subspace,2);
            start_p =sum(rank(1:(i_coil-1)),2)+1;
            end_p = sum(rank(1:(i_coil)),2);
            coil_position = start_p:end_p;
            
            senSubspaceProjectionL2(:,i_coil) = B(i_coil).subspace*recon_cg(coil_position);
        end
        
        
        %% display results
        % Objective
        if debug==1
            
            % projection
            opt.scale = [0,0.5];
            slice = 1;
            senSubspaceProjectionL2 = reshape(senSubspaceProjectionL2,dim_spatial(1:4));
            figure;montagesc(rot90(squeeze(senSubspaceProjectionL2(:,:,slice,:)).*repmat(brainMaskRef(:,:,slice),[1,1,dim_spatial(4),1])),opt);
            title(sprintf('Projection onto subspace\n (reference mask)'))
            colorbar;
            axis off
            MakeFigPretty(gcf)
            saveFigPng(gcf,fullfile(results_data,sprintf('Projection coil sensitivity map L2 slice %d',slice)));
            
        end
        close all;
        opt.senSubspaceProjectionL2 = senSubspaceProjectionL2;
        opt.recon_cg = recon_cg;
        opt.brainMaskRef = brainMaskRef;
    case 'p'
        %         clear classes
        undersampleMask = logical(abs(temp_undersample));
        % parameters for reconstruction
        % figure;montagesc(brainMaskRef);
        % brainMaskRef_trust = repmat(reshape(ones(size(brainMaskRef)),[],1),[1,size(testSen,2)]);
        param.E = subspace_model_error_mat(img,B,dim_spatial,undersampleMask);
        
        
        % L1Weight was chosen in an empirical way, but it seems to work very well
        % if the images are normalized from 0-1
        % param.W = TV_Temp();param.L1Weight=0.1;
        param.W = 1;param.L1Weight=0;
        param.TV = TVOP();param.TVWeight=0;
        param.y = temp_undersample;
        param.nite = 8;
        param.display=1;
        
        % linear reconstruction
        if ~isfield(opt,'recon_cg')
            
            
        elseif isempty(opt.recon_cg)
            recon_coe = param.E'*param.y;
        else
            recon_coe = opt.recon_cg;
        end
        % %% check the operator
        % checkAdjoint(@(x)param.E*x,...
        %     @(x)param.E'*x,...
        %     [size(recon_coe,1),1],1e-10,1);
        
        
        AHA = @(x)param.E'*(param.E*x)+opt.pro_beta.*(x./sqrt(opt.c_est_var));
        AHb = param.E'*param.y+opt.pro_beta.*(opt.c_est_mean./sqrt(opt.c_est_var));
        maxit = 200;
        tol = 1e-6;
        recon_cg = mypcg(AHA,AHb,tol,maxit,[],[],recon_coe);
        
        % k-t SPARSE-SENSE reconstruction
        fprintf('\n  sum of subsapce projection \n')
        
        % recon_coe_l2 = recon_coe;
        %% the non-linear CG algorithm will re-start two times to improve
        % convergence using the last result as starting point
        % for i_L1_iter = 1:1
        %     tic
        % [recon_coe_l2,cost_L2(i_L1_iter).loss] = sub_subspace_projectionL1NlCg(recon_coe_l2,param);
        % toc
        % end
        
        
        %% combination the coefficients and basis
        senSubspaceProjectionL2 =  zeros(size(B(1).subspace,1),size(B,2));
        for i_coil = 1:size(B,2)
            rank(i_coil) = size(B(i_coil).subspace,2);
            start_p =sum(rank(1:(i_coil-1)),2)+1;
            end_p = sum(rank(1:(i_coil)),2);
            coil_position = start_p:end_p;
            
            senSubspaceProjectionL2(:,i_coil) = B(i_coil).subspace*recon_cg(coil_position);
        end
        
        
        %% display results
        % Objective
        if debug==1
            
            % projection
            opt.scale = [0,0.5];
            slice = 1;
            senSubspaceProjectionL2 = reshape(senSubspaceProjectionL2,dim_spatial(1:4));
            figure;montagesc(rot90(squeeze(senSubspaceProjectionL2(:,:,slice,:)).*repmat(brainMaskRef(:,:,slice),[1,1,dim_spatial(4),1])),opt);
            title(sprintf('Projection onto subspace\n (reference mask)'))
            colorbar;
            axis off
            MakeFigPretty(gcf)
            saveFigPng(gcf,fullfile(results_data,sprintf('Projection coil sensitivity map L2 slice %d',slice)));
            
        end
        close all;
        opt.senSubspaceProjectionL2 = senSubspaceProjectionL2;
        opt.recon_cg = recon_cg;
        opt.brainMaskRef = brainMaskRef;
        
    case 'p+nn'
        if isfield(opt,'tol_cg')
    tol_cg = opt.tol_cg;
else
    tol_cg = 1e-5;
end

if isfield(opt,'maxit_cg')
    maxit_cg = opt.maxit_cg;
else
    maxit_cg = 500;
end

% Parse augmented lagrange minimization params
if isfield(opt,'a')
    a = opt.a;
else
    a = 1.2;
end

if isfield(opt,'b')
    b = opt.b;
else
    b = 1;
end

if isfield(opt,'tol_al')
    tol_al = opt.tol_al;
else
    tol_al = 1e-5;
end

if isfield(opt,'maxit_al')
    maxit_al = opt.maxit_al;
else
    maxit_al = 500;
end
reg_para_nn = opt.reg_para_nn;
if isfield(opt,'watch_svs')
    watch_svs = opt.watch_svs;
else
    watch_svs = false;
end
        undersampleMask = logical(abs(temp_undersample));
        my_rank = size(B(1).subspace,2);
        brainMaskRef_trust = repmat(reshape([reshape(brainMaskRef,[],1)],[],1),[1,size(B,2)]);
        dim_spatial_new = size(brainMaskRef_trust);
        
        param.E = subspace_model_error_mat(img,B,dim_spatial,undersampleMask);
        param.E_B = sub_subspace_mat(B,brainMaskRef_trust); % multicoil model (b1: coil sensitivities)
        
        param.W = 1;param.L1Weight=0;
%         param.TV = TVOP();param.TVWeight=0;
        param.y = temp_undersample;
        param.nite = 8;
        param.display=1;
        
        % linear reconstruction
        if ~isfield(opt,'recon_cg')
            
            
        elseif isempty(opt.recon_cg)
            recon_coe = param.E'*param.y;
        else
            recon_coe = opt.recon_cg;
        end
        % %% check the operator
        % checkAdjoint(@(x)param.E*x,...
        %     @(x)param.E'*x,...
        %     [size(recon_coe,1),1],1e-10,1);
        
        
        AHA = @(x)param.E'*(param.E*x)+opt.pro_beta.*(x./sqrt(opt.c_est_var));
        AHb = param.E'*param.y+opt.pro_beta.*(opt.c_est_mean./sqrt(opt.c_est_var));
        
        %% Solve for xt
        disp('Projection');
        % Set up for watching singular values
        trunc_opt.watch = false;
        if watch_svs
            figure;
            trunc_opt.watch = true;
            trunc_opt.ax = gca;
        end
        iter = 1;
        rel_err = inf;
        U    = zeros(dim_spatial_new);
        Z    = zeros(dim_spatial_new);
        my_A  = AHb;
        
        fprintf('Iteration: %03d    Rel Error: %10.4e\n',iter,rel_err);
        
        while rel_err > tol_al
            tic
            if iter > maxit_al
                fprintf('\nExceeded maximum number of iterations during nuclear norm minimization! ');
                fprintf('RelErr = %g\n',rel_err);
                break;
            end
            
            % Don't forget where you came from
            my_A_last = my_A;
            
            % Minimization via conjugate gradient descent
            [my_A,flag,relres,iter_pgc]  = ...
                pcg( @(x)(AHA(x) + .5*b*(param.E_B'*(param.E_B*x)) ), ...
                AHb + .5*(param.E_B'*(b.*U + Z)), tol_cg, maxit_cg, [], [], my_A );
            
            pcg_flag_check(flag,iter_pgc,relres);
            temp_sen = param.E_B*my_A;
            
            [U,Z] = local_low_rank_sen(temp_sen,reg_para_nn,b,Z,trunc_opt,dim_spatial);
            rel_err = norm(my_A - my_A_last)/norm(my_A_last);
            
            iter = iter + 1;
            b = a*b;
            toc
            fprintf('Iteration: %03d    Rel Error: %.4e\n',iter,rel_err);
        end
        fprintf('\n');
        clear AHd;
        % Compute normalized data consistency
%         resid           = temp_sen -  param.y;
%         opt.nconsist        = sum(abs(resid(:)).^2);
%         % Compute normalized nuclear norm regularization penalty
%         opt.npenalty_nn = local_low_rank_pen(temp_sen,reg_para_nn,dim_spatial);
%         % reconstruction
        fprintf('\n  sum of subsapce projection \n')
        recon_cg=my_A;
        
        %% combination the coefficients and basis
        senSubspaceProjectionL2 =  zeros(size(B(1).subspace,1),size(B,2));
        for i_coil = 1:size(B,2)
            rank(i_coil) = size(B(i_coil).subspace,2);
            start_p =sum(rank(1:(i_coil-1)),2)+1;
            end_p = sum(rank(1:(i_coil)),2);
            coil_position = start_p:end_p;
            
            senSubspaceProjectionL2(:,i_coil) = B(i_coil).subspace*recon_cg(coil_position);
        end
        
        
        %% display results
        % Objective
        if debug==1
            
            % projection
            opt.scale = [0,0.5];
            slice = 1;
            senSubspaceProjectionL2 = reshape(senSubspaceProjectionL2,dim_spatial(1:4));
            figure;montagesc(rot90(squeeze(senSubspaceProjectionL2(:,:,slice,:)).*repmat(brainMaskRef(:,:,slice),[1,1,dim_spatial(4),1])),opt);
            title(sprintf('Projection onto subspace\n (reference mask)'))
            colorbar;
            axis off
            MakeFigPretty(gcf)
            saveFigPng(gcf,fullfile(results_data,sprintf('Projection coil sensitivity map L2 slice %d',slice)));
            
        end
        close all;
        opt.senSubspaceProjectionL2 = senSubspaceProjectionL2;
        opt.recon_cg = recon_cg;
        opt.brainMaskRef = brainMaskRef;
        
    case 'l2+l1'
        %% L2+L1
        % parameters for reconstruction
        param.TVWeight = 0;
        % parameters for reconstruction
        brainMaskRef_trust = repmat(reshape(brainMaskRef,[],1),[1,size(testSen,2)]);
        param.E = sub_subspace_mat(B,brainMaskRef_trust); % multicoil model (b1: coil sensitivities)
        param.TV = TVOP();
        param.y = testSen.*brainMaskRef_trust;
        param.nite = 100;
        param.display=1;
        % test_beta = linspace(5e-5,5e-1,10);
        test_beta = 3e-5;
        cost_L2_L1 = [];
        my_metricSubspace = struct;
        % linear reconstruction
        recon_coe=param.E'*param.y;
        
        %% for multi param.L1Weight
        for i_beta = 1:length(test_beta)
            param.W = 1;param.L1Weight=test_beta(i_beta);
            % sum of subsapce projection
            fprintf('\n  sum of subsapce projection \n')
            
            recon_coe_L2_l1 = recon_coe;
            % the non-linear CG algorithm will re-start two times to improve
            % convergence using the last result as starting point
            for i_L1_iter = 1:10
                tic
                [recon_coe_L2_l1,cost_L2_L1(i_beta).loss] = sub_subspace_projectionL1NlCg(recon_coe_L2_l1,param);
                toc
            end
            
            %% combination the coefficients and basis
            sen_subspace_projection_L2_L1 =  zeros(size(B(1).subspace,1),size(B,2));
            for i_coil = 1:size(B,2)
                rank(i_coil) = size(B(i_coil).subspace,2);
                start_p =sum(rank(1:(i_coil-1)),2)+1;
                end_p = sum(rank(1:(i_coil)),2);
                coil_position = start_p:end_p;
                
                sen_subspace_projection_L2_L1(:,i_coil) = B(i_coil).subspace*recon_coe_L2_l1(coil_position);
            end
            
            % projection
            opt.scale = [0,0.8];
            slice = 1;
            sen_subspace_projection_L2_L1 = reshape(sen_subspace_projection_L2_L1,dim_spatial(1:4));
            figure;montagesc(rot90(squeeze(sen_subspace_projection_L2_L1(:,:,slice,:)).*repmat(brainMaskRef(:,:,slice),[1,1,dim_spatial(4),1])),opt);
            title(sprintf('Projection onto subspace\n (reference mask)'))
            colorbar;
            axis off
            MakeFigPretty(gcf)
            saveFigPng(gcf,fullfile(results_data,sprintf('Projection coil sensitivity map  L2_L1(weight %d) slice %d',param.L1Weight,slice)));
            
            %error
            test_sen_array = reshape(testSen,dim_spatial(1:4));
            slice = 1;
            opt.scale = [0,0.04];
            Error_projection_L2L1 = test_sen_array-sen_subspace_projection_L2_L1;
            figure;montagesc(rot90(squeeze(Error_projection_L2L1(:,:,slice,:)).*repmat(brainMaskRef(:,:,slice),[1,1,dim_spatial(4),1])),opt);
            title(sprintf('Error\n (reference mask)'))
            colorbar;
            axis off
            MakeFigPretty(gcf)
            saveFigPng(gcf,fullfile(results_data,sprintf('Error L2_L1(weight %d) %d',param.L1Weight,slice)));
            
            % Error vs iteration
            figure;plot(cost_L2_L1(i_beta).loss,'*');
            xlabel('Iteration');
            ylabel('Error');
            axis off
            title(sprintf('Error vs. iteration'))
            MakeFigPretty(gcf,[],5)
            saveFigPng(gcf,fullfile(results_data,sprintf('Error vs. iteration L2_L1(weight %d)',param.L1Weight)));
            
            
            %% calculate metrics
            % my_metric_CS_SENSE_BART(i_reg) = myAssessmentImageQuality(sen_subspace_projection_L2_L1,test_sen_array,opt);
        end
        close all;
        opt.senSubspaceProjectionL2 = sen_subspace_projection_L2_L1;
        
        % %% display the metric
        % my_metric_plot = struct2table(my_metric_CS_SENSE_BART);
        % figure;plot(test_beta(1:6),my_metric_plot.ssim);
        % xlabel('Beta');
        % ylabel('SSIM');
        % title(sprintf('SSIM vs.Beta'))
        % MakeFigPretty(gcf,[],5)
        % saveFigPng(gcf,fullfile(results_data,sprintf('SSIM vs.Beta')));
        %
        % figure;plot(test_beta(1:6),abs(my_metric_plot.psnr));
        % xlabel('Beta');
        % ylabel('PSNR');
        % title(sprintf('PSNR vs.Beta'))
        % MakeFigPretty(gcf,[],5)
        % saveFigPng(gcf,fullfile(results_data,sprintf('PSNR vs.Beta')));
        %
        % figure;plot(test_beta(1:6),abs(my_metric_plot.nmse));
        % xlabel('Beta');
        % ylabel('NMSE');
        % title(sprintf('NMSE vs.Beta'))
        % MakeFigPretty(gcf,[],5)
        % saveFigPng(gcf,fullfile(results_data,sprintf('NMSE vs.Beta')));
        %
        %
        % figure;plot(test_beta(1:6),abs(my_metric_plot.l1error));
        % xlabel('Beta');
        % ylabel('L1-Error');
        % title(sprintf('L1-Error vs.Beta'))
        % MakeFigPretty(gcf,[],5)
        % saveFigPng(gcf,fullfile(results_data,sprintf('L1-Error vs.Beta')));
        %
        %% display basis
        % opt.scale = [];
        % for i_order = 1:1:20
        % slice = 1;
        % for i_coil = 1:16
        %  temp = B(i_coil).subspace;
        %  temp = reshape(temp,[dim_spatial(1:3),size(temp,2)]);
        %   temp_basis(:,:,i_coil) = temp(:,:,slice,i_order);
        % end
        % figure;montagesc(rot90(temp_basis.*repmat(brainMaskRef(:,:,slice),[1,1,dim_spatial(4),1])),opt);
        % title(sprintf('Basis\n (order = %d)',i_order))
        % colorbar;
        % MakeFigPretty(gcf)
        % saveFigPng(gcf,fullfile(results_data,sprintf('Basis (order = %d) slice %d',i_order,slice)));
        % end
end
end


% %% L2+L1 SoftTresh for L1
% % parameters for reconstruction
% param.E = sub_subspace_mat(B); % multicoil model (b1: coil sensitivities)
% % L1Weight was chosen in an empirical way, but it seems to work very well
% % if the images are normalized from 0-1
% % param.W = TV_Temp();param.L1Weight=0.1;
% param.TV = TVOP();param.TVWeight=0;
% param.y = test_sen+my_noise;
% param.nite = 100;
% param.display=1;
% test_beta =0.005e-2;
% param.W = 1;
% % decreasing factor for lambda (lambda_ite=lambda*(1-ite*param.dec))
% param.dec = 0.01;
% param.nite=30;
% param.tol=5e-4;
% param.display=1;
%
% cost_L2_L1_soft_thresh = [];
% my_metric_L2L1_soft_thresh = [];
% % linear reconstruction
% % recon_coe=param.E'*param.y;
% recon_coe=zeros(size(param.E'*param.y));
% %% for multiparam.lambda
% for i_beta = 1:length(test_beta)
% param.lambda = test_beta(i_beta);
% % sum of subsapce projection
% fprintf('\n  sum of subsapce projection \n')
%
% recon_coe_L2_l1_soft_thresh = recon_coe;
% % the non-linear CG algorithm will re-start two times to improve
% % convergence using the last result as starting point
% for i_L1_iter = 1:3
% tic
% [recon_coe_L2_l1_soft_thresh,cost_L2_L1_soft_thresh(i_beta).loss] = sub_subspace_projectionL1SoftThresh(recon_coe_L2_l1_soft_thresh,param);
% toc
% end
%
%
% %% combination the coefficients and basis
% sen_subspace_projection_L2_L1_soft_thresh =  zeros(size(B(i_coil).subspace,1),size(B,2));
% for i_coil = 1:size(B,2)
% rank(i_coil) = size(B(i_coil).subspace,2);
%             coil_position = (sum(rank(1:(i_coil-1)),1)+1):(sum(rank(1:(i_coil)),1));
%
%     sen_subspace_projection_L2_L1_soft_thresh(:,i_coil) = B(i_coil).subspace*recon_coe_L2_l1_soft_thresh(coil_position);
% end
%
% % projection
% opt.scale = [];
% slice = 1;
% sen_subspace_projection_L2_L1_soft_thresh = reshape(sen_subspace_projection_L2_L1_soft_thresh,dim_maps(1:4));
% figure;montagesc(rot90(squeeze(sen_subspace_projection_L2_L1_soft_thresh(:,:,slice,:)).*repmat(brain_mask_all(:,:,slice,ind),[1,1,dim_maps(4),1])),opt);
% title(sprintf('Projection onto subspace\n (reference mask)'))
% colorbar;
% MakeFigPretty(gcf)
% saveFigPng(gcf,fullfile(results_data,sprintf('Projection coil sensitivity map  L2_L1(thresh %d) slice %d sof thresh',param.lambda,slice)));
%
% %error
% slice = 1;
% opt.scale = [0,0.04];
% Error_projection_L2L1_soft_thresh = test_sen_array-sen_subspace_projection_L2_L1_soft_thresh;
% figure;montagesc(rot90(squeeze(Error_projection_L2L1_soft_thresh(:,:,slice,:)).*repmat(brain_mask_all(:,:,slice,ind),[1,1,dim_maps(4),1])),opt);
% title(sprintf('Error\n (reference mask)'))
% colorbar;
% MakeFigPretty(gcf)
% saveFigPng(gcf,fullfile(results_data,sprintf('Error L2_L1(thresh %d) %d sof thresh',param.lambda,slice)));
%
% % Error vs iteration
% figure;plot(cost_L2_L1_soft_thresh(i_beta).loss);
% xlabel('Iteration');
% ylabel('Error');
% title(sprintf('Error vs. iteration'))
% MakeFigPretty(gcf,[],5)
% saveFigPng(gcf,fullfile(results_data,sprintf('Error vs. iteration L2_L1(thresh %d) sof thresh',param.lambda)));
%
%
% %% calculate metrics
% opt.metric = 'ssim';
% my_metric_L2L1_soft_thresh(i_beta).ssim = myAssessmentImageQuality(sen_subspace_projection_L2_L1_soft_thresh,test_sen_array,opt);
% opt.metric = 'nmse';
% my_metric_L2L1_soft_thresh(i_beta).nmse = myAssessmentImageQuality(sen_subspace_projection_L2_L1_soft_thresh,test_sen_array,opt);
% opt.metric = 'psnr';
% my_metric_L2L1_soft_thresh(i_beta).psnr = myAssessmentImageQuality(sen_subspace_projection_L2_L1_soft_thresh,test_sen_array,opt);
% opt.metric = 'l1error';
% my_metric_L2L1_soft_thresh(i_beta).l1error = myAssessmentImageQuality(sen_subspace_projection_L2_L1_soft_thresh,test_sen_array,opt);
%
% end
%
%
% %% display the metric
% my_metric_plot_soft_thresh = struct2table(my_metric_L2L1_soft_thresh);
% figure;plot(test_beta(1:10),my_metric_plot_soft_thresh.ssim);
% xlabel('Beta');
% ylabel('SSIM');
% title(sprintf('SSIM vs.Beta'))
% MakeFigPretty(gcf,[],5)
% saveFigPng(gcf,fullfile(results_data,sprintf('SSIM vs.Beta sof thresh')));
%
% figure;plot(test_beta(1:10),abs(my_metric_plot_soft_thresh.psnr));
% xlabel('Beta');
% ylabel('PSNR');
% title(sprintf('PSNR vs.Beta'))
% MakeFigPretty(gcf,[],5)
% saveFigPng(gcf,fullfile(results_data,sprintf('PSNR vs.Beta sof thresh')));
%
% figure;plot(test_beta(1:10),abs(my_metric_plot_soft_thresh.nmse));
% xlabel('Beta');
% ylabel('NMSE');
% title(sprintf('NMSE vs.Beta'))
% MakeFigPretty(gcf,[],5)
% saveFigPng(gcf,fullfile(results_data,sprintf('NMSE vs.Beta sof thresh')));
%
%
% figure;plot(test_beta(1:10),abs(my_metric_plot_soft_thresh.l1error));
% xlabel('Beta');
% ylabel('L1-Error');
% title(sprintf('L1-Error vs.Beta'))
% MakeFigPretty(gcf,[],5)
% saveFigPng(gcf,fullfile(results_data,sprintf('L1-Error vs.Beta sof thresh')));
