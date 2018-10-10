%% get the current version of CNMF-E 
version_file = fullfile(cnmfe_folder, 'versions.log'); 

tmp_f = fopen(version_file, 'r'); 

version_line = fgetl(tmp_f); 
commit_line = fgetl(tmp_f); 

disp(version_line); 

%% in case you are using git to manage your CNMF-E pacage 
try 
    [s,git_hash_string] = system('git rev-parse HEAD'); 
    commit_line = sprintf('git commit: %s', git_hash_string); 
    disp(commit_line); 
catch 
    disp(commit_line); 
end 