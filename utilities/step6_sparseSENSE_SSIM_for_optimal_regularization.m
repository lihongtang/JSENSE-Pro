function [image_out,opt] = step6_sparseSENSE_SSIM_for_optimal_regularization(kspace_slice,map_slice,groundtruth,undersampleMask,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using the BART CS-SENSE
% Discrete wavelet transform;
% Optimal regularization paramter with grid search with SSIM as metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
groundtruth = abs(groundtruth);
Nc = size(kspace_slice,4);
if ~exist('opt','var')
    opt = struct;
end
if isfield(opt,'debug')
    debug = opt.debug;
else
    debug = 0;
end

if isfield(opt,'reg')
    reg = opt.reg;
else
    reg = linspace(1e-4,1e-1,5);
end
if isfield(opt,'results_data')
    results_data = opt.results_data;
else
    results_data = pwd;
end

tempKsapce = opt.tempKsapce;
if length(reg) ==1
    opt.opt_reg = reg;
    eval(sprintf('image_out = SparseSENSE3D(map_slice,repmat(undersampleMask,[1,1,1,Nc]),kspace_slice,opt.imgCSSENSE_origin,%d);',reg))
    
    
    % CS-SENSE BART
    opt.myMetric= myAssessmentImageQuality(abs(image_out),groundtruth,opt);
    
    my_metric_CS_SENSE_BART= myAssessmentImageQuality(abs(image_out),groundtruth,opt);
    fprintf('reg 1\n')
    
else
    for i_reg = 1:length(reg)
        eval(sprintf('image_out = SparseSENSE3D(map_slice,repmat(undersampleMask,[1,1,1,Nc]),kspace_slice,opt.imgCSSENSE_origin,%d);',reg(i_reg)))
        
        % CS-SENSE BART
        my_metric_CS_SENSE_BART(i_reg)= myAssessmentImageQuality(abs(image_out),groundtruth,opt);
        fprintf('reg %d\n',i_reg)
        figure;montagesc(image_out);
        figure;montagesc(groundtruth);
        SENSE3DOp = SENSE3D(map_slice);
        update_k = Fn_x2k(SENSE3DOp*image_out,1:3);
        l2norm_error(i_reg) = norm(update_k(:)-tempKsapce(:),2);
        
        
    end
    % display the metric
    my_metric_plot = struct2table(my_metric_CS_SENSE_BART);
    reg_log = log10(reg);
    if debug == 1
        figure;
        subplot(2,2,1);plot(reg_log,my_metric_plot.ssim);
        hold on;plot(reg_log,my_metric_plot.ssim,'k*');
        hold on;plot(reg_log,my_metric_plot.mssim,'r*')
        xlabel('log(\mu)');
        ylabel('SSIM');
        MakeFigPretty(gcf)
        
        subplot(2,2,2);
        plot(reg_log,abs(my_metric_plot.psnr),'*');
        xlabel('log(\mu)');
        ylabel('PSNR');
        MakeFigPretty(gcf,[],5)
        
        subplot(2,2,3);
        plot(reg_log,abs(l2norm_error),'r*');
        xlabel('log(\mu)');
        ylabel('L2Error');
        MakeFigPretty(gcf,[],5)
        
        
        subplot(2,2,4);
        plot(reg_log,abs(my_metric_plot.nmse),'r*');
        xlabel('log(\mu)');
        ylabel('NMSE');
        MakeFigPretty(gcf)
        saveFigPng(gcf,fullfile(results_data,sprintf('metric vs mu %d iter %d',opt.i_subj,opt.iteration_algorithm)));
    end
    while(1)
        %     locs = opt_regularization_parameter(my_metric_plot);
        %     % use the mssim as th emetric
        %% l2 error as the best parameetrs;
        [~,locs] = findpeaks(-l2norm_error);
        [~,temp_locs] = max(-l2norm_error(locs));
        locs = locs(temp_locs);
        if isempty(locs)
            disp('Original reg!')
            reg
            reg_value = input('Do not find the optimal regularization parameters!\n Please input([min,max,num] or number)');
            if length(reg_value)>1
                reg = 10.^linspace(reg_value(1),reg_value(2),reg_value(3));
                reg_log = log10(reg);
            else
                locs = reg_value;
                break;
                
            end
            
            for i_reg = 1:length(reg)
                eval(sprintf('image_out = SparseSENSE3D(map_slice,repmat(undersampleMask,[1,1,1,Nc]),kspace_slice,opt.imgCSSENSE_origin,%d);',reg(i_reg)))
                
                my_metric_CS_SENSE_BART(i_reg) = myAssessmentImageQuality(abs(image_out),groundtruth,opt);
                SENSE3DOp = SENSE3D(map_slice);
                update_k = Fn_x2k(SENSE3DOp*image_out,1:3);
                l2norm_error(i_reg) = norm(update_k(:)-tempKsapce(:),2);
                fprintf('reg %d',i_reg)
            end
            % display the metric
            my_metric_plot = struct2table(my_metric_CS_SENSE_BART);
            
            if debug == 1
                figure;
                subplot(2,2,1);plot(reg_log,my_metric_plot.ssim);
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
                plot(reg_log,abs(l2norm_error));
                xlabel('log(\mu)');
                ylabel('L2Error');
                MakeFigPretty(gcf,[],5)
                
                
                subplot(2,2,4);
                plot(reg_log,abs(my_metric_plot.nmse));
                xlabel('log(\mu)');
                ylabel('NMSE');
                MakeFigPretty(gcf)
                saveFigPng(gcf,fullfile(results_data,sprintf('metric vs.mu %s iter %d',opt.i_subj,opt.iteration_algorithm)));
            end
        else
             
            break;
        end
    end
    eval(sprintf('image_out = SparseSENSE3D(map_slice,repmat(undersampleMask,[1,1,1,Nc]),kspace_slice,opt.imgCSSENSE_origin,%d);',reg(locs)))
    opt.myMetric= myAssessmentImageQuality(abs(image_out),groundtruth,opt);
    opt.opt_reg = reg(locs);
end


opt.imgCSSENSE = image_out;
opt.my_metric_CS_SENSE_BART = my_metric_CS_SENSE_BART;
end