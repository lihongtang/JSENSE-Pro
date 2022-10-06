function loadStruct = myLoadVars(fileName,varargin) 
%% ======================================================================
%MYLOADVARS a) load variable(s) to base only if they don't exist in base
%           b) load variables to structure
%           c) load one variable to a variable
%Note: be careful when you are to change the variables to load!
% =======================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-09-12
%
%   [INPUTS]
%   ---- Required ----
%   fileName                name of the file to load data from
%   dataName                name of variables to load, either a string or a cell of strings
%
%   ---- Optional ----
%   verbose                 degree of verbosity [1]
%
%   [OUTPUTS]
%   loadStruct              a structure containing the loaded data
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/09/12
%       Modified by Yibo Zhao @ UIUC, 2018/09/15:
%           Fix a bug when loading multiple variable with a struct as
%           output, using a new function "mergeStructs".
%       Modified by Yibo Zhao @ UIUC, 2018/10/23:
%           Separete the case of loading to base and loading to a structure.
%           If there is a single variable with single output, directly assign
%           the loaded data to the input (instead of a structure of it).
%       Modified by Yibo Zhao @ UIUC, 2018/10/24:
%           Warn the user when the data file doesn't exist, instead of issuing an error.
%       Modified by Yibo Zhao @ UIUC, 2018/10/25:
%           Parse inputs more gracefully, and check for file extensions.
%
%   Example:
%       myLoadVars(preparedDataFile,'D2Info'); % load to base
%       t_struct     = myLoadVars(preparedDataFile,'tD2_o','tD2_e'); % load to structure
%       sxt_epsi_odd = myLoadVars(preparedDataFile,'sxt_sense_odd'); % load to variable
%
%   See also MERGESTRUCTS
%
%--------------------------------------------------------------------------

%% ------ parse the input ------
    [dataName,verbose] = parseInputs(varargin{:});
    % Rules: 
    %   first scalar variable goes to verbose
    %   cell of strings and string variable go into structure varName
    %   other variables are ignored
    
%% ------ load variables ------
    fileName = check_and_add_mat_extension(fileName);
    
    if ~exist(fileName,'file') % if still doesn't exist
        if verbose>0, warning('The specified file %s doesn''t exist.\nCannot load any data.',fileName);end
        return;
    end
    
    if nargout==0 % if load to base
        for nn = 1:length(dataName) % for each data
            if ~existInBase(dataName{nn}) % if it is not loaded yet
                if verbose>0, fprintf(['Loading ',dataName{nn},'...']);tic; end
                
                loadStruct_temp = load(fileName,dataName{nn}); %#ok<NASGU>
                
                try
                    temp = eval(['loadStruct_temp.',dataName{nn}]);
                    assignin('base',dataName{nn},temp);
                catch
                    continue; % the user should already be warned
                end
                
                if verbose>0, fprintf('\b\b\b, done! (%ds)\n',round(toc)); end
            else
                if verbose>0, disp([dataName{nn},' already exists in the workspace.']); end
            end
        end
    else % if load to struct
        if length(dataName)==1
            if verbose>0, fprintf(['Loading ',dataName{1},'...']);tic; end
            
            loadStruct_temp = load(fileName,dataName{1});
            loadStruct = loadStruct_temp.(dataName{1});
            
            if verbose>0, fprintf('\b\b\b, done! (%ds)\n',round(toc)); end
        else
            loadStruct = struct;
            for nn = 1:length(dataName) % for each data
                if verbose>0, fprintf(['Loading ',dataName{nn},'...']);tic; end

                loadStruct_temp = load(fileName,dataName{nn});
                loadStruct      = mergeStructs(loadStruct,loadStruct_temp);

                if verbose>0, fprintf('\b\b\b, done! (%ds)\n',round(toc)); end
            end
        end
    end
    
end

function fileName = check_and_add_mat_extension(fileName)

    if ~exist(fileName,'file') % if doesn't exist
        if ~strcmpi(fileName(end-3:end),'.mat'); % if doesn't have a ',mat' extension
            fileName = strcat(fileName,'.mat'); % add file extension '.mat'
        end
    end
    
end

function [dataName,verbose] = parseInputs(varargin)
    
    narginchk(1,inf);
    
    % default values
    verbose    = 1;

    % parse args    
    flag_verbose_parsed = false;
    
    for ii = 1:nargin
        
        % scalar should be verbose
        if isscalar(varargin{ii})&&(~ischar(varargin{ii}))&&(~iscellstr(varargin{ii}))&&(~flag_verbose_parsed) 
            verbose             = varargin{ii};
            flag_verbose_parsed = true;
        
        % matrix should be weights
        elseif ischar(varargin{ii})
            if ~exist('dataName','var')
                dataName             = {varargin{ii}};
            else
                dataName             = [dataName,varargin{ii}];
            end
            
        elseif iscellstr(varargin{ii})
            if ~exist('dataName','var')
                dataName             = varargin{ii};
            else
                dataName             = [dataName,varargin{ii}{:}];
            end
            
        % unexpected case
        else
            warning('Input #%d has illegal data type. Please check.',ii);
        end
    end
    
end

%#ok<*AGROW>
