function [ data ] = opSampleH( sample_data, mask  )
%OPSAMPLEH Hermatian conjugate of sampling operator (see opSample.m)
%   Represents the sampling operation of data collection
%   
%   data    the input data to be sampled
%   mask    the k-space sampling mask
%--------------------------------------------------------------------------

    data = zeros(size(mask));
    data(mask) = sample_data;
    
    return;

end
