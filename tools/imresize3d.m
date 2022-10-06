function [ c3d_up, orig_coords, new_coords ] = imresize3d( c3d, newDim, meth, threshold )
% 3D Interpolation
% ---- INPUTS ----
% c3d       : original image
% newDim    : new dimensions
% meth      : interpolation method (default is spline)
%
% ---- OUTPUTS ----
% c3d_up        : interpolated image
% orig_coords   : coordinates for the original image (used for interpolation)
% new_coords    : coordinates for the new image (used for interpolation)
% meth          : interpolation method, default spline
% threshold     : threshold to generate binary output, default 0.5
%
% Author: Bryan A. Clifford
% Date Created: 10-17-2014
% Last Updated: Chao Ma, 08-08-2015, Add threhold
% ------------------------------------------------------------------------------

    %% ---- 3D interpolation ----
    if nargin < 3
        meth = 'spline';
    end
    if nargin < 4
        threshold = 0.5;
    end

    dims    = size(c3d);
    if length(dims) == 2
        dims(3) = 1;
    end
    if length(newDim) == 2
        newDim(3) = 1;
    end

    % [gy_orig, gx_orig, gz_orig] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
    ny = (0:dims(1)-1) + .5;
    nx = (0:dims(2)-1) + .5;
    nz = (0:dims(3)-1) + .5;
    [gy_orig, gx_orig, gz_orig] = ndgrid(ny, nx, nz);

    ny_new = ( ( 0:newDim(1)-1 ) + .5 )*dims(1)/newDim(1);
    nx_new = ( ( 0:newDim(2)-1 ) + .5 )*dims(2)/newDim(2);
    nz_new = ( ( 0:newDim(3)-1 ) + .5 )*dims(3)/newDim(3);
    [gy_new, gx_new, gz_new] = ndgrid(ny_new, nx_new, nz_new);

    orig_coords.x = gx_orig;
    orig_coords.y = gy_orig;
    orig_coords.z = gz_orig;
    new_coords.x = gx_new;
    new_coords.y = gy_new;
    new_coords.z = gz_new;

    if dims(3) > 1
        c3d_up = interp3(gx_orig, gy_orig, gz_orig, double(c3d), gx_new, gy_new, gz_new, meth);
    else
        c3d_up = interp2(gx_orig, gy_orig, double(c3d), gx_new, gy_new, meth);
    end

    if islogical(c3d)
        c3d_up = c3d_up >= threshold;
    end

    return;
end


