function S_new = SetStructDefault(S, fieldname, defaultval)
%   Function:       S_new = SETSTRUCTDEFAULT(S, fieldname, defaultval)
%   Author:         Qiang Ning
%   Created:        Jul 08, 2014
%   Last update:    Aug 19, 2014
%   Description:    Set default values for missing fields
%   Required external files:
%                   (none)
%   Update history: Qiang@Aug 19, 2014. Changed header file; allowed
%                   "fieldname" not to be cell-typed; added a feature that
%                   S might have one field but its val is empty.
%   Input:
%       S--the struct to be initialized.
%       fieldname--a single string or a cell of strings;fields that needs 
%           to be initialized.
%       defaultval--if the corresponding field does not exists, set it to
%           be the default value specified by "defaultval".
%   Output:
%       S_new: struct with default vals
%   Example:
%       If S=struct('a',1,'b',2),
%       Let S_new=SETSTRUCTDEFAULT(S, {'b','c'},{3,4}),
%       Then S_new=struct('a',1,'b',2,'c',4).
%       Let S_new2=SETSTRUCTDEFAULT(S, 'c', 4),
%       Then S_new2=struct('a',1,'b',2,'c',4).

%%  Check inputs
if iscell(fieldname)
    if ~iscell(defaultval)
        error('defaultval must be a cell.');
    end
    if length(fieldname) ~= length(defaultval)
        error('Dimension mismatch.');
    end
    for ii = 1 : length(fieldname)
        if ~ischar(fieldname{ii})
            error('fieldname must be string(s).');
        end
    end
    mode        = 1;
else
    if ~ischar(fieldname)
        error('fieldname must be a string.');
    end
    if ~isrow(fieldname)
        error('fieldname must be a row string.');
    end
    mode        = 0;
end

%%  Set default vals
if mode
    for kk = 1 : length(fieldname)
        if ~isfield(S, fieldname{kk}) ||...
                isempty(eval(sprintf('S.%s',fieldname{kk})))
            eval(sprintf('S.%s=defaultval{kk};',fieldname{kk}));
        end
    end
else
    if ~isfield(S, fieldname) ||...
                isempty(eval(sprintf('S.%s',fieldname)))
        eval(sprintf('S.%s=defaultval;',fieldname));
    end
end
S_new = S;