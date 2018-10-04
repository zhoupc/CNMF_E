%% select file 
if ~exist('zip_file_path', 'var') || ~exist(zip_file_path, 'file')
    [file_name, folder_path] = uigetfile('*.zip'); 
    zip_file_path = fullfile(folder_path, file_name); 
end 

%% unzip the zip file
[s_dir, s_name, ~] = fileparts(zip_file_path);
output_folder = fullfile(s_dir, s_name);
unzip(zip_file_path, output_folder); 

%% select the mat file 
temp = dir(output_folder); 
for m=1:length(temp)
   [~, ~, file_type] = fileparts(temp(m).name); 
   if strcmpi(file_type, '.mat')
      matfile_name = fullfile(output_folder, temp(m).name); 
   elseif strcmpi(file_type, '.txt')
       textfile_name = fullfile(output_folder, temp(m).name); 
   end
end
if ~exist(matfile_name, 'file')
    error('no mat file in the selected zip file'); 
end 

%% load mat file 
load(matfile_name); 

% check whether we used the same computer again 
if exist(neuron.P.log_file, 'file')
    movefile(matfile_name, neuron.P.log_folder); 
    movefile(textfile_name, neuron.P.log_file); 
    rmdir(output_folder); 
else
    original_log_folder = neuron.P.log_folder;
    original_logfile = neuron.P.log_file;
    neuron.P.log_folder = output_folder;
    neuron.P.log_file = textfile_name;
end


fprintf('The results have been loaded. \n'); 