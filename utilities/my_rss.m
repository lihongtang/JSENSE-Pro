function img_rss = my_rss(dx,coil)
if nargin < 2
    img_rss = sqrt(sum(abs(dx).^2,4));
else
    img_rss = squeeze(sqrt(sum(abs(dx).^2,coil)));
end