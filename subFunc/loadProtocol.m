function protocol = loadProtocol(filename)
% LOADPROTOCOL is a custom "load" function that automatically updates the
% "SaveDir" path when the hard drive letter is changed by the OS.
[Path , filename] = fileparts(filename);
if ~isfolder(Path)
    Path = fileparts(which(filename));
end    
a = load(fullfile(Path, filename));
protocol = a.protocol; clear a
hd_letter = Path(1);
if ~startsWith(protocol.SaveDir, hd_letter)
    % Update just the HD letter for now:
    protocol.SaveDir(1) = hd_letter;
end
end
