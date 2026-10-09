function filePath = writeCSVConditionFile(folderPath, fileName, spec)
%WRITECSVCONDITIONFILE Write a synthetic CSV condition file.
%
% Syntax:
%   filePath = writeCSVConditionFile(folderPath, fileName, labels)
%   filePath = writeCSVConditionFile(folderPath, fileName, tableStruct)
%
% Inputs:
%   folderPath - Destination folder.
%   fileName   - CSV file name.
%   spec       - Either:
%                  1) cellstr/string vector with one label per event, or
%                  2) struct with one field per CSV column and one element
%                     per event.
%
% Output:
%   filePath - Full path to the written CSV file.

    folderPath = convertStringsToChars(folderPath);
    fileName = convertStringsToChars(fileName);
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end

    if isstring(spec)
        spec = cellstr(spec(:));
    end

    if iscell(spec)
        T = table(spec(:), 'VariableNames', {'Condition'});
    elseif isstruct(spec)
        fn = fieldnames(spec);
        assert(~isempty(fn), 'The specification struct must have at least one field.');
        nRows = numel(spec.(fn{1}));
        vars = cell(1, numel(fn));
        for ii = 1:numel(fn)
            col = spec.(fn{ii});
            if isstring(col)
                col = cellstr(col(:));
            elseif ischar(col)
                col = cellstr(string(col(:)));
            else
                col = cellstr(string(col(:)));
            end
            assert(numel(col) == nRows, 'All CSV columns must have the same number of rows.');
            vars{ii} = col(:);
        end
        T = table(vars{:}, 'VariableNames', fn);
    else
        error('Unsupported CSV specification format.');
    end

    filePath = fullfile(folderPath, fileName);
    writetable(T, filePath);
end
