function nam_mat = tif2mat(nam)
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
end

info = imfinfo(nam);

      
if  isfield(info(1), 'ImageDescription') && ~isempty(strfind(info(1).ImageDescription,'ImageJ'))
    junk1=regexp(info(1).ImageDescription,'images=\d*','match');
    if isempty(junk1)
        junk1=regexp(info(1).ImageDescription,'frames=\d*','match');
    end
    junk2=strjoin(junk1);
    T=strread(junk2,'%*s %d','delimiter','=');
else
    T = length(info);   % number of frames
end
d1 = info.Height;   % height of the image
d2 = info.Width;    % width of the image
Ysiz = [d1, d2, T]'; 

fprintf('CNMF_E is converting TIFF file to *.mat file'); 
% create a mat file 
Tchunk = min(T, round(2^29/d1/d2)); %each chunk uses at most 4GB
Y = smod_bigread2(nam, 1, Tchunk);  %#ok<*NASGU>
save(nam_mat, 'Y', 'Ysiz', '-v7.3'); 
if Tchunk==T
    return; 
else
    data = matfile(nam_mat, 'Writable', true); 
    t0 = Tchunk+1; 
    while t0<=T
        num2read = min(t0+Tchunk-1, T) - t0 + 1; 
        tmpY = smod_bigread2(nam, t0, num2read); 
        data.Y(:, :, (1:num2read)+t0-1) = tmpY; 
        t0 = t0 + num2read; 
    end 
end 