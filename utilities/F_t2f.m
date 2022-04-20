function [ spectrum ] = F_t2f( signal )
%F_T2F Converts 2D spatial-temporal data to spatial-spectral data
%
% signal    The spatial-temporal data. First 2 dims should be space.
%           (Spatial can mean k-space or x-space)
%--------------------------------------------------------------------------

    spectrum = fftshift(fft( signal, [], 3 ), 3)/sqrt( size(signal,3) );
    return;

end

