function safeFclose(fid)
    if fid > 0 & fopen(fid) ~= -1 %#ok
        fclose(fid);
    end
end
