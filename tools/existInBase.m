function X = existInBase(varName) 
% check whether a variable exist in the base
%   Created by Yibo Zhao @ UIUC, 09/12/2018
    
    X = true;
    
    try
        evalin('base',[varName,';']);
    catch
        X = false;
    end
    
end