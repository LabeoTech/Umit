function result = exportMatToJson(matFile, jsonFile, options)
%EXPORTMATTOJSON Write a readable JSON copy of a report or log MAT file.
%
%   result = exportMatToJson(matFile, jsonFile)
%   result = exportMatToJson(matFile, jsonFile, 'MaxArrayElements', 200, 'Full', false)
%
%   Reads every variable of MATFILE and writes one UTF-8 JSON document, so a
%   report or log written by DataViewer (the alignment report, pipeLog,
%   dataHistory, the motion-correction log, ...) can be read without MATLAB.
%   The source file is never modified.
%
%   The converter is available only in the compiled application
%   (ISDEPLOYED is true): with MATLAB the MAT file can be read directly.
%   Otherwise it raises Umitoolbox:MatToJson:NotDeployed unless 'Force' is
%   true (developers and tests).
%
%   Document layout:
%       { "sourceFile": "<file name>", "exportedOn": "<ISO 8601 time>",
%         "summarisedVariables": [ ... ], "variables": { "<name>": <value>, ... } }
%
%   Encoding rules:
%       struct            object; a struct array is an array of objects
%       table             { "_type": "table", "variableNames": [...], "rows": [ {...}, ... ] }
%       cell              array
%       char / string     text (a char matrix or string array is an array of text)
%       logical, numeric  number or (nested) array; NaN and Inf are written as null
%       datetime          ISO 8601 text; duration, categorical, enum: text
%       function handle   its text; containers.Map: object of its keys
%       other objects     object of their readable properties plus "_class"
%                         (the class name alone when none is readable)
%   An array with more than 'MaxArrayElements' elements (numeric, logical,
%   char matrix, string, datetime, cell or struct array; default 200) is
%   replaced by a descriptor { "_type": "array", "class": ..., "size": [...],
%   "min": ..., "max": ..., "mean": ..., "note": "values left out ..." } and
%   the variable's name is listed in summarisedVariables. A table is never
%   summarised row-wise, but the cells inside it follow the same rule.
%   'Full', true writes every value.
%
%   Options:
%       MaxArrayElements - element limit before an array is summarised. Default 200.
%       Full             - write every value. Default false.
%       Force            - allow the conversion outside the compiled application.
%                          Default false.
%
%   Output (struct): JsonFile (path written), SummarisedVariables (string
%   column), VariableNames (string column).
%
%   Errors: Umitoolbox:MatToJson:NotDeployed, :MissingFile, :NotAMatFile,
%   :SameFile (the target is the source file), :CannotWrite.

arguments
    matFile (1, 1) string
    jsonFile (1, 1) string
    options.MaxArrayElements (1, 1) double {mustBePositive, mustBeInteger} = 200
    options.Full (1, 1) logical = false
    options.Force (1, 1) logical = false
end

if ~isdeployed && ~options.Force
    error('Umitoolbox:MatToJson:NotDeployed', ...
        ['JSON export is available in the compiled application only. With MATLAB, ' ...
         'open the file with load("%s").'], matFile);
end

if ~isfile(matFile)
    error('Umitoolbox:MatToJson:MissingFile', 'File not found: %s', matFile);
end

matFullPath = iFullPath(matFile);
jsonFullPath = iFullPath(jsonFile);
if strcmpi(matFullPath, jsonFullPath)
    error('Umitoolbox:MatToJson:SameFile', ...
        'The JSON file must not be the MAT file itself: %s', matFullPath);
end

try
    variableInfo = whos('-file', matFile);
    loaded = load(matFile);
catch ME
    error('Umitoolbox:MatToJson:NotAMatFile', ...
        '"%s" could not be read as a MAT file: %s', matFile, ME.message);
end

limit = Inf;
if ~options.Full
    limit = options.MaxArrayElements;
end

% Keys that are not valid MATLAB field names are written with this prefix and
% renamed in the text afterwards.
typeKey = 'JSONKEYTYPE';
classKey = 'JSONKEYCLASS';

summarised = strings(0, 1);
currentVariable = "";

names = string({variableInfo.name}).';
variables = struct();
for iName = 1:numel(names)
    currentVariable = names(iName);
    variables.(names(iName)) = encode(loaded.(names(iName)));
end
summarised = unique(summarised, 'stable');

[~, matName, matExt] = fileparts(matFullPath);
document = struct( ...
    'sourceFile', [matName matExt], ...
    'exportedOn', char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')), ...
    'summarisedVariables', {cellstr(summarised(:).')}, ...
    'variables', variables);

text = jsonencode(document, 'PrettyPrint', true, 'ConvertInfAndNaN', true);
text = strrep(text, ['"' typeKey '"'], '"_type"');
text = strrep(text, ['"' classKey '"'], '"_class"');

fid = fopen(jsonFullPath, 'w', 'n', 'UTF-8');
if fid < 0
    error('Umitoolbox:MatToJson:CannotWrite', 'Could not write: %s', jsonFullPath);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, unicode2native(text, 'UTF-8'), 'uint8');
clear cleanup

result = struct( ...
    'JsonFile', jsonFullPath, ...
    'SummarisedVariables', summarised, ...
    'VariableNames', names);

    % ---------------------------------------------------------------------
    % Encoding (nested: shares limit, summarised and currentVariable)
    % ---------------------------------------------------------------------

    function value = encode(x)
        if isa(x, 'timetable')
            x = timetable2table(x);
        end
        if isa(x, 'table')
            value = encodeTable(x);
        elseif isstruct(x)
            value = encodeStruct(x);
        elseif iscell(x)
            value = encodeCell(x);
        elseif isa(x, 'datetime')
            if numel(x) > limit
                value = descriptor(x);
            else
                x.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSS';
                value = textValue(string(x));
            end
        elseif isa(x, 'duration') || isa(x, 'calendarDuration') || isa(x, 'categorical')
            if numel(x) > limit
                value = descriptor(x);
            else
                value = textValue(string(x));
            end
        elseif isstring(x)
            if numel(x) > limit
                value = descriptor(x);
            else
                value = textValue(x);
            end
        elseif ischar(x)
            value = encodeChar(x);
        elseif islogical(x) || isnumeric(x)
            value = encodeNumeric(x);
        elseif isa(x, 'function_handle')
            value = string(func2str(x));
        elseif isa(x, 'containers.Map')
            value = encodeMap(x);
        elseif isenum(x)
            value = textValue(string(x));
        elseif isobject(x)
            value = encodeObject(x);
        else
            value = string(class(x));
        end
    end

    function value = encodeStruct(s)
        if isempty(s)
            value = {};
        elseif isscalar(s)
            value = struct();
            fields = fieldnames(s);
            for iField = 1:numel(fields)
                value.(fields{iField}) = encode(s.(fields{iField}));
            end
        elseif numel(s) > limit
            value = descriptor(s);
        else
            value = cell(1, numel(s));
            for iElement = 1:numel(s)
                value{iElement} = encodeStruct(s(iElement));
            end
        end
    end

    function value = encodeTable(t)
        variableNames = t.Properties.VariableNames;
        rows = cell(1, height(t));
        for iRow = 1:height(t)
            row = struct();
            for iVar = 1:numel(variableNames)
                column = t.(variableNames{iVar});
                if iscell(column)
                    cellValue = column{iRow};
                else
                    cellValue = column(iRow, :);
                end
                row.(variableNames{iVar}) = encode(cellValue);
            end
            rows{iRow} = row;
        end
        value = struct(typeKey, 'table', 'variableNames', {variableNames}, 'rows', {rows});
    end

    function value = encodeCell(c)
        if isempty(c)
            value = {};
        elseif numel(c) > limit
            value = descriptor(c);
        else
            value = cell(1, numel(c));
            for iElement = 1:numel(c)
                value{iElement} = encode(c{iElement});
            end
        end
    end

    function value = encodeChar(x)
        if isempty(x)
            value = '';
        elseif isrow(x)
            value = x;
        elseif numel(x) > limit
            value = descriptor(x);
        else
            value = cellstr(x).';
        end
    end

    function value = encodeNumeric(x)
        if numel(x) > limit
            value = descriptor(x);
        elseif ~isreal(x)
            value = descriptor(x);
        else
            if issparse(x)
                x = full(x);
            end
            if isempty(x)
                value = [];
            elseif islogical(x)
                value = x;
            else
                value = double(x);
            end
        end
    end

    function value = encodeMap(m)
        value = struct();
        keyList = keys(m);
        for iKey = 1:numel(keyList)
            key = keyList{iKey};
            fieldName = matlab.lang.makeValidName(char(string(key)));
            value.(fieldName) = encode(m(key));
        end
    end

    function value = encodeObject(obj)
        value = struct();
        try
            propertyNames = properties(obj);
            for iProp = 1:numel(propertyNames)
                name = propertyNames{iProp};
                try
                    value.(name) = encode(obj.(name));
                catch
                    value.(name) = '<unreadable>';
                end
            end
        catch
            propertyNames = {};
        end
        if isempty(propertyNames)
            value = string(class(obj));
        else
            value.(classKey) = class(obj);
        end
    end

    function value = textValue(s)
        s = s(:).';
        s(ismissing(s)) = "";
        if isempty(s)
            value = {};
        elseif isscalar(s)
            value = s;
        else
            value = cellstr(s);
        end
    end

    function d = descriptor(x)
        %DESCRIPTOR Small description of an array that is too large to write.
        d = struct(typeKey, 'array', 'class', class(x), 'size', size(x));
        if isnumeric(x) || islogical(x)
            values = double(x(:));
            values = values(isfinite(values));
            if ~isempty(values)
                d.min = min(values);
                d.max = max(values);
                d.mean = mean(values);
            end
        end
        d.note = ['values left out (more than ' num2str(limit) ...
            ' elements); use the Full option to write them'];
        summarised(end + 1, 1) = currentVariable;
    end
end

function p = iFullPath(file)
%IFULLPATH Absolute path of a file that may not exist yet.
p = char(file);
[folder, name, ext] = fileparts(p);
if isempty(folder)
    folder = pwd;
end
if isfolder(folder)
    [~, info] = fileattrib(folder);
    if isstruct(info) && isfield(info, 'Name')
        folder = info.Name;
    end
end
p = fullfile(folder, [name ext]);
end
