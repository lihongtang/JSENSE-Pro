function [ signal ] = F3_f2t( spectrum )
%F3_F2T Converts 3D spatial-spectral data to spatial-temporal data
%
% spectrum    The spatial-spectral data. First 3 dims should be space.
%             (Spatial can mean k-space or x-space)
%--------------------------------------------------------------------------

    signal = ifft(ifftshift( spectrum, 4 ), [], 4)*sqrt( size(spectrum,4) );
    return;

end

