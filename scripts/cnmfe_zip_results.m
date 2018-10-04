%% 
if ~exist('original_log_folder', 'var') 
   original_log_folder = neuron.P.log_folder; 
   original_logfile = neuron.P.log_file; 
end

if ~exist('zip_file_path', 'var')
   error('zip_file_path has not been set. Please assign its value first.'); 
elseif exist(zip_file_path, 'file')
    fprintf('There a zip file with the same name. \nDo you want to overwrite it?\n'); 
    temp = input('Yes(Y)/No(N): ', 's');
    if ~strcmpi(temp, 'y')
        zip_file_path = fullfile(neuron.P.log_folder, [get_minute(), '.zip']); 
    end
end
neuron.save_workspace(zip_file_path, original_log_folder, original_logfile); 