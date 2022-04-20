function [ signal ] = F_f2t( spectrum )
%F_F2T Converts 2D spatial-spectral data to spatial-temporal data
%
% spectrum    The spatial-spectral data. First 2 dims should be space.
%             (Spatial can mean k-space or x-space)
%--------------------------------------------------------------------------

    signal = ifft(ifftshift( spectrum, 3 ), [], 3)*sqrt( size(spectrum,3) );
    return;

end

