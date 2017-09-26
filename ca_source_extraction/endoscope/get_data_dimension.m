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
    data = matfile(file_to_file); 
    dims = data.Ysiz; 
else
    error('Unknown file extension. Only .tiff and .hdf5 files are currently supported');
    dims = [];
end
dims = double(dims);