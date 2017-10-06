function dims = get_data_dimension(path_to_file)
% get he dimension of a tiff file. now it only supports 2D

% Author: Pengcheng Zhou, Columbia University, 2017

%get image info


[~,~,ext] = fileparts(path_to_file);

if strcmpi(ext,'.tiff') || strcmpi(ext,'.tif');
    info = imfinfo(path_to_file);
    
    if isfield(info(1), 'ImageDescription') && ~isempty(strfind(info(1).ImageDescription,'ImageJ'))
        junk1=regexp(info(1).ImageDescription,'images=\d*','match');
        if isempty(junk1)
            junk1=regexp(info(1).ImageDescription,'frames=\d*','match');
        end
        junk2=strjoin(junk1);
        temp =textscan(junk2,'%*s %d','delimiter','=');
        numFrames = temp{1};
        
        d1 = info(1).Height;
        d2 = info(1).Width;
    else
        blah=size(info);
        numFrames= blah(1);
        d1 = info.Height;
        d2 = info.Width;
    end
    dims = [d1, d2, numFrames];
elseif strcmpi(ext,'.hdf5') || strcmpi(ext,'.h5');
    info = hdf5info(path_to_file);
    dims = info.GroupHierarchy.Datasets.Dims;
elseif strcmpi(ext, '.avi')
    obj = audiovideo.mmreader(path_to_file);
    
    frame_rate = obj.FrameRate;
    total = obj.Duration;
    numFrames = total*frame_rate;
    
    sizx = obj.Width;
    sizy = obj.Height;
    dims = [sizy, sizx, numFrames];
elseif strcmpi(ext, '.mat')
    data = matfile(path_to_file); 
    dims = data.Ysiz; 
elseif isempty(ext)
    % the input is a folder and data are stored as image sequences
    imgs = dir([path_to_file, filesep,'*.tif']);
    
    trialandframe = zeros(length(imgs),2);

    % fill out (trial and) frame list for sorting/excluding frames
    firstInd = 1;
    for file = 1:length(imgs)
        fnam = imgs(file).name;
        IDs = regexp(fnam,'\d*\d*','match'); % use a regular expression to find trial index and frame index within each file name
        trialandframe(file,1)=str2num(IDs{1})+firstInd-1;
        trialandframe(file,2)=str2num(IDs{2});
    end
    
    % throw out trials with less than trial_length frames 
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
    
%   find indices of frames to throw out (i.e. the first 5 frames of every
%   trial, etc.)
    excludeIDs = (trialandframe(:,2) == 0);
    
    % sort frames in ascending order before reading
    trialandframe(excludeIDs,:) = [];
    [~, srt] = sortrows(trialandframe,[1 2]);
    
    imgs(excludeIDs)=[]; 
    imgs = imgs(srt);
    
    %get image info from first image of sequence
    first_img= fullfile(path_to_file,imgs(1).name);
    
    info = imfinfo(first_img);
    sizx = info.Height;
    sizy = info.Width;
    
    numFrames=length(imgs);
    
    dims = [sizy, sizx, numFrames];

else
    error('Unknown file extension. Only .tiff and .hdf5 files are currently supported');
    dims = [];
end
dims = double(dims);