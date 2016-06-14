%% edit the file directory
temp = get(edit_dir, 'string');
if exist([dir_nm, temp], 'file')
    set(edit_dir, 'string', temp);
    [~, file_nm, file_type] = fileparts(temp); %#ok<ASGLU>
else 
    set(edit_file, 'string', [file_nm, file_type]); 
    return; 
end

%% convert tif to mat file and map it to memory  
gui_load_mat; 
set(edit_height, 'string', d1); 
set(edit_width, 'string', d2); 
set(edit_frame, 'string', numFrame); 

%% update file types and file names 
file_type = '.mat'; 
nam = nam_mat; 
set(edit_file, 'string', [file_nm, file_type]); 