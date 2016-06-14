%% edit the file directory 
temp = get(edit_dir, 'string'); 
if exist(dir_nm, 'dir')
    if (temp(end))~=filesep
        temp(end+1) = filesep; 
    end
    dir_nm = temp; 
else
    set(edit_dir, 'string', dir_nm); 
end