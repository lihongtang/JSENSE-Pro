function [ var_est, info ] = estNoiseVar( noise, opts )
%ESTNOISEVAR Estimates the noise variance in a MRSI data set
% ---- INPUTS ----
% noise         2D matrix with noise samples. Each row should be a noise sample from
%               a different spatial location
% opts.disp     (optional default = 1) A flag specifying that the noise histogram
%               and covariance matrix should be displayed
% opts.siglvl   (optional default = 0.05) The significance level to use for the
%               confidence interval estimation
% opts.get_cv   (optional default = true) If true, estimate the covariance
%               matrix of the noise sample.
%             
% ---- OUTPUTS ----
% var_est             The estimated variance of the noise
% info.conf_intvl     The confidence interval
% info.ks_real_p      The estimated p value for the KS-Test of the real data
%                     to a normal distribution
% info.ks_imag_p      The estimated p value for the KS-Test of the imaginary
%                     data to a normal distribution
% info.ks_h           The result of the hypothosis test. 0 if both the real
%                     and imaginary parts of the noise passed the KS-Test.
% info.cvmat          The estimated spatial covariance matrix of the noise
% ------------------------------------------------------------------------------

    %% Parse the input arguments
    if ndims(noise) ~= 2
        error('ERROR: noise data has wrong dimensions!');
    end
    [M N] = size(noise);


    if nargin > 1
        if isfield(opts,'disp')
            doPlot = opts.disp;
        else
            doPlot = 1;
        end

        if isfield(opts,'siglvl')
            siglvl = opts.siglvl;
        else
            siglvl = 0.05;
        end
        
        if isfield(opts,'get_cv')
            get_cv = opts.get_cv;
        else
            get_cv = true;
        end

        clear opts;
    else
        doPlot = 1;
        siglvl = 0.05;
        get_cv = true;
    end


    %% Estimate noise standard deviation
    var_est = var(noise(:));
    d = N*M - 1;
    a = (d*var_est)/chi2inv(1 - siglvl/2, d);
    b = (d*var_est)/chi2inv(siglvl/2, d);
    fprintf('estimated noise variance: %g\n',var_est);
    fprintf('estimated %6.4f%% confidence interval: [%g, %g]\n',1-siglvl,a,b);


    %% Estimate Noise Spatial Covariance Matrix
    %---------------------------------------------------------------------------
    if get_cv
        info.cvmat = zeros([M M]);
        fprintf('estimating covariance matrix of noise...');
        info.cvmat = noise*noise'/(N-1);
        fprintf('done\n');
    end
    
    %% Get real and imaginary noise
    noise_r = sort(real(noise(:)).','ascend');
    noise_r = (noise_r - mean(noise_r))/std(noise_r);

    noise_i = sort(imag(noise(:)).','ascend');
    noise_i = (noise_i - mean(noise_i))/std(noise_i);


    %% Display Noise Covariance Matrix
    if doPlot && get_cv
        figure;

        subplot(2,1,1);
        imagesc(real(info.cvmat));
        title('\bfEstimated Covariance Matrix (Real Part)');
        axis square off;
        colormap(gray(256));
        colorbar;

        subplot(2,1,2);
        imagesc(imag(info.cvmat));
        title('\bfEstimated Covariance Matrix (Imaginary Part)');
        axis square off;
        colormap(gray(256));
        colorbar;
    end


    %% Estimate Noise Distribution
    if doPlot
        figure;
        title('\bfEstimated Noise Distribution');
        ylabel('noise pdf')
        xlabel('noise value')
        hold all;

        Nbins = 100;
        nmax  = max( [noise_r, noise_i] );
        nmin  = min( [noise_r, noise_i] );
        bin_width = (nmax - nmin)/Nbins;
        bin_edges = (0:Nbins+1)*bin_width + nmin - bin_width/2;
        bin_ctrs  = bin_edges + bin_width/2;
        
        hr = histc(noise_r,bin_edges);
        hi = histc(noise_i.',bin_edges);

        assert(sum(hr) == (M*N));
        assert(sum(hi) == (M*N));

        pdfr = hr/(M*N)/bin_width;
        pdfi = hi/(M*N)/bin_width;

        stairs(bin_ctrs,pdfr,'b')
        stairs(bin_ctrs,pdfi,'r')
        plot(bin_ctrs, exp(-(bin_ctrs.^2)/2)/sqrt(2*pi), 'k');

        legend('Real Part','Imaginary Part','Normal Distribution')
    end
    

    %% Compare the distribution of the noise to the normal distribution
    [r_pass info.ks_real_p] = kstest(noise_r,[],.05);
    [i_pass info.ks_imag_p] = kstest(noise_i,[],.05);
    if ~r_pass
        fprintf('real part of estimated noise PASSED ks-test      (p-value: %f)\n',info.ks_real_p);
    else
        fprintf('real part of estimated noise FAILED ks-test      (p-value: %f)\n',info.ks_real_p);
    end

    if ~i_pass
        fprintf('imaginary part of estimated noise PASSED ks-test (p-value: %f)\n',info.ks_imag_p);
    else
        fprintf('imaginary part of estimated noise FAILED ks-test (p-value: %f)\n',info.ks_imag_p);
    end
    
    
    %% Package data for output
    info.conf_intvl = [a, b];
    info.ks_pass = r_pass*i_pass;

end

