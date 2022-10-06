function [ sample_data ] = opSample(data, mask )
%OPSAMPLE Sampling operator
%   Represents the sampling operation of data collection
%   
%   data    the input data to be sampled
%   mask    the logical sampling mask
%--------------------------------------------------------------------------
    
    sample_data = data(mask);
    
    return

end

