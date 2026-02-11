function safeFclose(fid)
    if fid > 0 & fopen(fid) ~= -1
        fclose(fid);
    end
end
