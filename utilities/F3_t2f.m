function [ spectrum ] = F3_t2f( signal )
%F3_T2F Converts 3D spatial-temporal data to spatial-spectral data
%
% signal    The spatial-temporal data. First 3 dims should be space.
%           (Spatial can mean k-space or x-space)
%--------------------------------------------------------------------------

    spectrum = fftshift(fft( signal, [], 4 ), 4)/sqrt( size(signal,4) );
    return;

end

