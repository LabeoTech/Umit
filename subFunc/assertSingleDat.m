function assertSingleDat(source, varargin)
%ASSERTSINGLEDAT Error unless a .dat file stores single-precision values.
%
%   assertSingleDat(filePath)
%   assertSingleDat(Info)
%   assertSingleDat(..., 'Context', text)
%
%   umIT supports only single-precision .dat files: the header can describe
%   other classes (uint8, uint16, int16, uint32, int32, double), but the
%   toolbox neither writes nor processes them. This is the one shared gate
%   for that policy; PipelineManager and the file rewriters call it.
%
%   Inputs:
%       source - Either the path of a .dat file (its class is read with
%                loadMetaData, so headered and legacy sidecar files are
%                both covered) or an Info struct with a dataClass field.
%                A path that is not a .dat file, or a file that does not
%                exist, is ignored: callers check existence themselves.
%
%   Name-Value options:
%       Context - Text that starts the error message, for example
%                 'Step 3 (myFunction) output "outData"'. Default: the file
%                 name, or 'The data'.
%
%   Errors:
%       Umitoolbox:Dat:nonSingle - the file or Info stores another class.
%       Anything loadMetaData raises for an unreadable file is not caught
%       here: this function enforces the class policy only.
%
%   See also: loadMetaData, saveData, datHeaderFromInfo

p = inputParser;
p.FunctionName = 'assertSingleDat';
addParameter(p, 'Context', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, varargin{:});
context = char(string(p.Results.Context));

if isstruct(source)
    if ~isscalar(source) || ~isfield(source, 'dataClass')
        error('Umitoolbox:assertSingleDat:invalidInput', ...
            'Info must be a scalar struct with a dataClass field.');
    end
    dataClass = char(string(source.dataClass));
    if isempty(context)
        context = 'The data';
    end
else
    if ~(ischar(source) || (isstring(source) && isscalar(source)))
        error('Umitoolbox:assertSingleDat:invalidInput', ...
            'The source must be a .dat file path or an Info struct.');
    end
    filePath = char(string(source));
    [~, fileName, extension] = fileparts(filePath);
    if ~strcmpi(extension, '.dat') || ~isfile(filePath)
        return
    end
    dataClass = char(string(loadMetaData(filePath).dataClass));
    if isempty(context)
        context = sprintf('"%s%s"', fileName, extension);
    end
end

if ~strcmpi(dataClass, 'single')
    error('Umitoolbox:Dat:nonSingle', ...
        ['%s stores %s values. umIT supports only single-precision .dat files: ' ...
        'convert the data to single before it is written (for example single(data)).'], ...
        context, dataClass);
end
end
