function txt = tableToTextConverter(tab)
% Convert a MATLAB table to a formatted text representation.
%
% This function takes a MATLAB table 'tab' as input and converts it into a
% text-based representation. It handles text and numeric data, while
% non-convertible items (e.g., structures, cell arrays) are replaced with a
% placeholder text.
%
% Inputs:
%   tab (table) - The input MATLAB table to be converted.
%
% Outputs:
%   txt (char) - A character array containing the text-based representation of
%   the input table.
%
% Notes:
%   - The input table 'tab' should be a valid MATLAB table.
%   - The function assumes that the table's data is well-structured and
%   compatible for conversion.
%
% Example:
%   inputTable = readtable('data.csv'); % Load a sample table from a CSV file
%   textTable = tableToTextConverter(inputTable);
%   disp(textTable); % Display the text-based representation of the table

% Check if all elements of the table are convertible to text format:
headers = tab.Properties.VariableNames;
tabC = table2cell(tab);
bConvertibleItems = cellfun(@(x) isnumeric(x) || ischar(x), tabC, 'UniformOutput', true);

% Replace non-convertible items with a placeholder text:
placeHolderText = 'Non-convertible';
tabC(~bConvertibleItems) = {placeHolderText};

% Convert all numeric data to text:
tabC = cellfun(@(x) num2str(x), tabC, 'UniformOutput', false);
tabC = [headers; tabC];

% Calculate maximum width for all columns
colWidth = cellfun(@(x) length(x), tabC);
colWidth = max(colWidth, [], 1);

% Center Headers:
for ii = 1:size(tabC, 2)
    tabC(1, ii) = pad(tabC(1, ii), colWidth(ii) + 2, 'both');
    tabC(2:end, ii) = pad(tabC(2:end, ii), colWidth(ii) + 1, 'right');
    tabC(2:end, ii) = pad(tabC(2:end, ii), colWidth(ii) + 2, 'left');
end

% Create text-formatted table:
colWidth = cellfun(@(x) length(x), tabC(1, :));
txt = '';
headerSep = '';

% Create headers:
for ii = 1:size(tabC, 2)
    headerSep = [headerSep sprintf('+%s', repmat('-', 1, colWidth(ii)))];%#ok
end
headerSep(end + 1) = '+';

txt = [txt sprintf('%s\n|%s|\n%s\n', headerSep, strjoin(tabC(1, :), '|'), headerSep)];

for ii = 2:size(tabC, 1)
    txt = [txt sprintf('|%s|\n', strjoin(tabC(ii, :), '|'))];%#ok
end
txt = [txt sprintf('%s\n', headerSep)];
end
