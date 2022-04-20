function mySaveVars(fileName,varargin) 
%% ========================================================================
%MYSAVEVARS save data into a file, either creating a new file or appending
%to an existing file, while skipping non-existent variables.
% =========================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-09-12
%
%   [INPUTS]
%   ---- Required ----
%   fileName                name of the file to save data into
%   varName                 name of variables to save, either a string or a cell of strings
%
%   ---- Optional ----
%   verbose                 degree of verbosity [0]
%
%   [OUTPUTS]
%   N/A
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/09/12
%       Modified by Yibo Zhao @ UIUC, 2018/10/22
%           If one of the variables to save doesn't exist, assign [] to it
%       Modified by Yibo Zhao @ UIUC, 2018/10/25
%           Parse inputs more gracefully
%
%--------------------------------------------------------------------------

%% ------ parse the input ------
    [varName,verbose] = parseInputs(varargin{:});
    % Rules: 
    %   first scalar variable goes to verbose
    %   cell of strings and string variable go into structure varName
    %   other variables are ignored
    
%% ------ save variables ------
    % This code is kind of ugly, since in MATLAB the values of variables are not
    % transferred into the function, so we have to save data in the base...
    fileName = check_and_add_mat_extension(fileName);

    if verbose>0
        disp(['Saving data to ',fileName]);t_save = tic;
    end
    
    % If a variable doesn't exist in base, then assign empty to it, warn the user, and save.
    for nn = 1:length(varName) % for each variable
        if ~existInBase(varName{nn})
            warning('%s not found in base. Assigning empty matrix to it and save...',varName{nn});
            assignin('base',varName{nn},[]);
        end
    end
    
    % If the file already exists, use '-append'.
    if exist(fileName,'file')
        name_temp = 'THISISAVERYLONGVARIABLETOAVOIDDECLAREDVARIABLES'; % this is a very long variable to avoid declared variables
        assignin('base',name_temp,varName);
        evalin('base',['save(''',fileName,''',''-append'',',name_temp,'{:});']);
        evalin('base',['clear ',name_temp,';']);
        
    % If the file doesn't exist, use '-v7.3'.
    else
        name_temp = 'THISISAVERYLONGVARIABLETOAVOIDDECLAREDVARIABLES';
        assignin('base',name_temp,varName);
        evalin('base',['save(''',fileName,''',''-v7.3'',',name_temp,'{:});']);
        evalin('base',['clear ',name_temp,';']);
    end
    
    if verbose>0
        disp(['Finished! Elapsed time: ',num2str(toc(t_save)),'.']);
    end
    
end

function fileName = check_and_add_mat_extension(fileName)

    if ~exist(fileName,'file') % if doesn't exist
        if ~strcmpi(fileName(end-3:end),'.mat'); % if doesn't have a ',mat' extension
            fileName = strcat(fileName,'.mat'); % add file extension '.mat'
        end
    end
    
end

function [varName,verbose] = parseInputs(varargin)
    
    narginchk(1,inf);
    
    % default values
    verbose    = 0;

    % parse args    
    flag_verbose_parsed = false;
    
    for ii = 1:nargin
        
        % scalar should be verbose
        if isscalar(varargin{ii})&&(~ischar(varargin{ii}))&&(~iscellstr(varargin{ii}))&&(~flag_verbose_parsed) 
            verbose             = varargin{ii};
            flag_verbose_parsed = true;
        
        % matrix should be weights
        elseif ischar(varargin{ii})
            if ~exist('varName','var')
                varName             = {varargin{ii}};
            else
                varName             = [varName,varargin{ii}];
            end
            
        elseif iscellstr(varargin{ii})
            if ~exist('varName','var')
                varName             = varargin{ii};
            else
                varName             = [varName,varargin{ii}{:}];
            end
            
        % unexpected case
        else
            warning('Input #%d has illegal data type. Please check.',ii);
        end
    end
    
end

%#ok<*AGROW>
