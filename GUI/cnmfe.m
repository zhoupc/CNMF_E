clear; clc; close all; 
default_options = CNMFSetParms(); 

%% PACKAGE FOLDER 
cur_cd = cd();      % current directory 
[temp, ~, ~] = fileparts(mfilename('fullpath')); 
cd(temp); cd ..;  
CNMFE_DIR = pwd(); 
cd(cur_cd); 

try
    load([CNMFE_DIR, filesep, '.dir.mat']); %load previous path
catch
    dir_nm = cd(); %use the current path
end

%% create main figure of the GUI window 
cnmfe_main; 

%% panel for including all information of the raw data 
cnmfe_data; 

%% panel for choosing neuron parameters
cnmfe_neuron; 

%% panel for loading the data 
cnmfe_load; 

%% initialization  
cnmfe_init; 

%% update background 
cnmfe_background; 