function w = genWeights2d(My,Mx,Nt,dirs,Iref,Iref_conf,percent_max_weight,display_weight_stats)
    
    ndirs  = length(dirs);
    c = zeros([My Mx ndirs]);
    for m = 1:size(Iref,3)
        cur_Iref    = Iref(:,:,m);

        cur_Iref    = cur_Iref - min(cur_Iref(:));
        cur_Iref    = cur_Iref./max(cur_Iref(:));

        cur_c       = zeros(My, Mx, ndirs);
        for n=1:ndirs
            cur_c(:,:,n) = abs(diff_along_dir(cur_Iref, dirs(n))).^2;
        end

        c = c + Iref_conf(m)*cur_c;
    end

    % Normalize weights
    wmax = quantile(1./c(:),1-percent_max_weight);
    if wmax == inf
        disp('Warning: wmax equals inf');
    end
    w = double(c <= 1/wmax);
    w( w == 0 ) = 1./c( w == 0 )/wmax;

    % Plot statistics about the weights
    if display_weight_stats
        figure

        subplot(2,1,1)
        hist(w(:),1000)
        title('\bfWeight Distribution')

        subplot(2,1,2)
        plot((1:numel(w))/numel(w),sort(w(:)))
        title('\bfSorted Weights'),xlabel('weight index percentile');

        figure, imagesc(mean(c,3)), colormap(gray(256)), title('\bfMean Squared Finite Difference Values');
        figure, imagesc(mean(w,3)), colormap(gray(256)), title('\bfMean Weight Values');
    end

    % Replicate weights for all time points
    w = reshape(w, [My*Mx ndirs]);
    w = repmat(w,  [Nt 1]);
    w = reshape(w, [My Mx Nt ndirs]);

end