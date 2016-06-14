%% change file folders
dir_nm = uigetdir(dir_nm, 'select file folder');
if dir_nm~=0  % save the current result 
    if dir_nm(end)~=filesep; dir_nm(end+1) = filesep; end
    
    save([CNMFE_DIR, filesep, '.dir.mat'],  'dir_nm');
    set(edit_dir, 'string', dir_nm); 
else  % fail to select new folder, use the old one 
    dir_nm = get(edit_dir, 'string'); 
end

