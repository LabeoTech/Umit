function out = ReadConfigFile()

Folder = mfilename('fullpath');
Folder = Folder(1:strfind(Folder,'\ReadConfig'));

CList = dir([Folder filesep 'cfg*.txt']);

if( isempty(CList) )
    warndlg(['No configuration file found. Please contact '...
        'info@labeotech.com to get the config file for your setup.'],...
        'Configuration File');
    ME = MException('Configuration:FileMissing', ...
        ['No configuration file found. Please contact '...
        'info@labeotech.com to get the config file for your setup.']);
    throw(ME);
elseif( size(CList,1) > 1 )
    str = {CList.name};
    [s, v] = listdlg('PromptString', 'Multiple config files were found. Select one:',...
        'SelectionMode','single','ListString',str);
    if( v < 1 )
        ME = MException('Configuration:FileSelection', ...
        'Multiple configuration files found. One must be selected.');
        throw(ME);
    end
    fid = fopen([Folder CList(s).name]);
    out.ConfigName = CList(s).name;
else
    fid = fopen([Folder CList.name]);
    out.ConfigName = CList.name;
end


while ~feof(fid)
    tline = fgetl(fid);
%    disp(tline)
    Pos = strfind(tline, ':');
    if( ~isempty(Pos) )
        tline(1:Pos) = regexprep(tline(1:Pos), ' ', '_');
        Param = tline(1:(Pos-1));
        Value = (tline((Pos+2):end));
        eval(['out.' Param ' = ''' Value ''';']);
    end
end