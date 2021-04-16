function getNamesOfOpenFigs

fH = get(groot);

names = sort({fH.Children.Name});
fid = fopen('OpenFigNames.txt', 'w');
fprintf(fid, '%s\n', 'List of Opened Figures:');
for i = 1:length(names)
       fprintf(fid, '%s\n', names{i});
end
fclose(fid);
winopen(pwd);
end
