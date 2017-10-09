function nam_mat = tif2mat_sequence(nam)
%% convert tiff files into mat files
% inputs:
%   nam: file names
% output:
%   nam_mat: name of the *.mat file
% Author: Pengcheng Zhou, Carnegie Mellon University

%% tiff file information 
if isempty(nam)
    nam_mat = [];
    fprintf('empty tiff file');
    return;
else
    [tmp_dir, tmp_file, ~] = fileparts(nam);
    nam_mat = sprintf('%s%s%s.mat', tmp_dir, filesep, tmp_file);
    Y = smod_bigread2_trials(nam);
    Ysiz = size(Y);
    save(nam_mat, 'Y', 'Ysiz', '-v7.3');  
end