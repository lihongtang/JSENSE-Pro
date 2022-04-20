function [ data ] = F_x2k( image, dont_shift )
%F_X2K Converts 2D spatial-spectral data from image space to data space.
%
% image         The spatial-spectral image. First 2 dims should be image space.
%               The 3rd dim can be time or frequency.
% 
% dont_shift    (Optional) If true, don't do any fft shifts.
%--------------------------------------------------------------------------

    if nargin < 2
        dont_shift = false;
    end
    
    [Ny, Nx, ~] = size(image);
    if dont_shift
        data = sqrt(1/(Nx*Ny))*fft(fft(image,[],1),[],2);
    else
        data = sqrt(1/(Nx*Ny))*ifftshift(fft(fftshift(...
                                    ifftshift(fft(fftshift(...
                                        image,...
                                     1),[],1),1),...
                                2),[],2),2);
    end
    return;

end

