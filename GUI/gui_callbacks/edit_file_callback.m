%% edit the file directory
temp = get(edit_dir, 'string');
if exist([dir_nm, temp], 'file')
    set(edit_dir, 'string', temp);
    [~, file_nm, file_type] = fileparts(temp); %#ok<ASGLU>
else 
    set(edit_file, 'string', [file_nm, file_type]); 
    return; 
end

%% convert tif to mat file and map it to memory, create Sources2D object 
gui_load_mat; 
set(edit_height, 'string', d1); 
set(edit_width, 'string', d2); 
set(edit_frame, 'string', numFrame); 

neuron_raw = Sources2D('d1',d1,'d2',d2, 'bas_nonneg', 1, ...
    'gSig', 4, 'gSiz', 15);   % dimensions of datasets
sframe = 1; 
eframe = numFrame; 
num2read = eframe-sframe+1; 
set(edit_begin, 'string', 1); 
set(edit_end, 'string', eframe); 
set(edit_total, 'string', num2read); 

%% update file types and file names 
file_type = '.mat'; 
nam = nam_mat; 
set(edit_file, 'string', [file_nm, file_type]); 