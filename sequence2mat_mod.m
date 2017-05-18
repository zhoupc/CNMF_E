%% Reads in sequence of .tif files and memory maps them to a .mat file. Assumes images are in NeuralMapper output format:
% (i.e. unsigned 16-bit .tifs with naming format: 'image_XXX_XXX.tif', where XXX is  an integer
% R. Conor Heins, 05022017
% edit 1.1, Conor and Tarun Madangopal, 05072017

function [data,Ysiz,sortedTrialsandFrames] = sequence2mat_mod(data,add_data_flag)

%choose image directory (folder with image sequence)
fprintf('Choose a folder to load new image data from...\n')
stack_dir = uigetdir();
imgs = dir([stack_dir, filesep,'*.tif']);

%choose which frames to throw away (e.g. first frame of every trial, which is often very low fluorescence)
frames2exclude = input('Are there any frames you would like to exclude?\n(enter single digit or vector, or press enter if you want to keep everything):\n');

%read in trial and frame indexes from .tiff names
trialandframe = zeros(length(imgs),2);
firstInd = input('Enter number of the first trial in this stack (default = 1):\n');
if isempty(firstInd)
    firstInd = 1;
end
for file = 1:length(imgs)
    fnam = imgs(file).name;
    IDs = regexp(fnam,'\d*\d*','match');
    trialandframe(file,1)=str2num(IDs{1})+firstInd-1;
    trialandframe(file,2)=str2num(IDs{2});
end

% throw out trials with less than 800 frames 
firstCol = trialandframe(:,1);
allTrials = unique(firstCol);
badTrials = zeros(size(firstCol));
for i = 1:length(allTrials)
    if length(trialandframe(firstCol==allTrials(i),:)) < 800
        badTrials = badTrials + firstCol==allTrials(i);
    end
end
trialandframe(logical(badTrials),:) = [];
imgs(logical(badTrials))=[];

excludeIDs = find(trialandframe(:,2)==frames2exclude);
imgs(excludeIDs)=[]; trialandframe(excludeIDs,:) = [];

[sortedTrialsandFrames, sortInds] = sortrows(trialandframe,[1 2]);
sortedImgs = imgs(sortInds);

tic;
fprintf('converting the selected image sequence to *.mat version...\n');

info = imfinfo(fullfile(stack_dir,sortedImgs(1).name)); %read first image from sequence to get image info
d1 = info.Height;   % height of the image 
d2 = info.Width;    % width of the image 
T = length(sortedImgs);
Ysiz = [d1, d2, T]';

if add_data_flag
    firstFrame = size(data.Y,3)+1;
else
    data.Y = uint16(zeros(d1,d2,2)); %clunky, but initialize length of 3rd dimension to 2 so that we can add to the third dimension later
    firstFrame = 1;
end

Tchunk = min(T, round(2^29/d1/d2)); %each chunk uses at most 4GB

data.Y(:,:,firstFrame:(firstFrame+Tchunk-1)) = sequence_bigread2(stack_dir,sortedImgs,1,Tchunk);
if Tchunk==T
    return; 
else
    t0 = Tchunk+1; 
    %t0 = Tchunk
    while t0<=T
        num2read = min(t0+Tchunk-1, T) - t0 + 1; 
        tmpY = sequence_bigread2(stack_dir,sortedImgs, t0, num2read); 
        data.Y(:, :, (firstFrame-1)+(1:num2read)+t0-1) = tmpY; 
        t0 = t0 + num2read; 
    end 
end

fprintf('Time cost in converting data to *.mat file:     %.2f seconds\n', toc);

end