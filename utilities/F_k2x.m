function [ image ] = F_k2x( data, dont_shift )
%F_K2X Converts 2D spatial-spectral data from data space to image space.
%
% image         The spatial-spectral image. First 2 dims should be data space.
%               The 3rd dim can be time or frequency.
% 
% dont_shift    (Optional) If true, don't do any fft shifts.
%--------------------------------------------------------------------------

    if nargin < 2
        dont_shift = false;
    end

    [Ny, Nx, ~] = size(data);
    if dont_shift
        image = sqrt(Nx*Ny)*ifft(ifft(data,[],1),[],2);
    else
        image = sqrt(Nx*Ny)*ifftshift(ifft(fftshift(...
                                ifftshift(ifft(fftshift(...
                                    data,...
                                1),[],1),1),...
                            2),[],2),2);
    end
    return;

end
