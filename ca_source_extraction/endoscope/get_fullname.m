%% function for get files full name by replacing relative locations like '.', '~'
function nam = get_fullname(nam)
tmp_dir = cd();
[dir_nm, file_nm, ext] = fileparts(nam);
if isempty(dir_nm)
    dir_nm = tmp_dir;
end
cd(dir_nm);
dir_nm = cd();
nam = [dir_nm,filesep, file_nm, ext];
cd(tmp_dir);
end