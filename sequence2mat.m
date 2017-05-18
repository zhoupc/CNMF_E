%% Reads in sequence of .tiff files and writes them to .tif stack. Assumes images are in NeuralMapper output format:
% (i.e. unsigned 16-bit .tifs with naming format: 'image_XXX_XXX.tif', where XXX is
% an integer
% R. Conor Heins, 05032017

function [stack_dir, nam_mat, data, sortedTrialsandFrames] = sequence2mat(excludeFrame)

%choose image directory (folder with image sequence)
stack_dir = uigetdir();
imgs = dir([stack_dir, filesep,'*.tif']);

%choose which frames to throw away (e.g. first frame of every trial, which is often very low fluorescence)
frames2exclude = excludeFrame;

%read in trial and frame indexes from .tiff names
trialandframe = zeros(length(imgs),2);
for file = 1:length(imgs)
    fnam = imgs(file).name;
    IDs = regexp(fnam,'\d*\d*','match');
    trialandframe(file,1)=str2num(IDs{1});
    trialandframe(file,2)=str2num(IDs{2});
end

excludeIDs = find(trialandframe(:,2)==frames2exclude);
imgs(excludeIDs)=[]; trialandframe(excludeIDs,:) = [];

[sortedTrialsandFrames, sortInds] = sortrows(trialandframe,[1 2]);
sortedImgs = imgs(sortInds);

tic;
fprintf('converting the selected image sequence to *.mat version...\n');

nam_mat = [stack_dir, '.mat'];

info = imfinfo(fullfile(stack_dir,sortedImgs(1).name)); %read first image from sequence to get image info
d1 = info.Height;   % height of the image 
d2 = info.Width;    % width of the image 
T = length(sortedImgs);
Ysiz = [d1, d2, T]';

Tchunk = min(T, round(2^29/d1/d2)); %each chunk uses at most 4GB
Y = sequence_bigread2(stack_dir,sortedImgs,1,Tchunk);
save(nam_mat, 'Y', 'Ysiz','sortedTrialsandFrames','-v7.3');
if Tchunk==T
    return; 
else
    data = matfile(nam_mat, 'Writable', true); 
    t0 = Tchunk+1; 
    while t0<=T
        num2read = min(t0+Tchunk-1, T) - t0 + 1; 
        tmpY = sequence_bigread2(stack_dir,sortedImgs, t0, num2read); 
        data.Y(:, :, (1:num2read)+t0-1) = tmpY; 
        data.Y = [cat(3,data.Y,tmpY)];
        t0 = t0 + num2read; 
    end 
end

fprintf('Time cost in converting data to *.mat file:     %.2f seconds\n', toc);

end
    


