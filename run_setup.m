%% add path 
cnmf_dir = fileparts(which('run_setup.m')); 
addpath(sprintf('%s%s%', cnmf_dir, filesep)); 
addpath(sprintf('%s%s%utilities', cnmf_dir, filesep)); 
addpath(sprintf('%s%s%endoscope', cnmf_dir, filesep)); 
