% convert the data to mat file
nam_mat = [dir_nm, file_nm, '.mat'];
if strcmpi(file_type, '.mat')
    fprintf('The selected file is *.mat file\n');
elseif  exist(nam_mat', 'file')
    % the selected file has been converted to *.mat file already
    fprintf('The selected file has been replaced with its *.mat version\n');
elseif or(strcmpi(file_type, '.tif'), strcmpi(file_type, '.tiff'))
    % convert
    tic;
    fprintf('converting the selected file to *.mat version...\n');
    nam_mat = tif2mat(nam);
    fprintf('Time cost in converting data to *.mat file:     %.2f seconds\n', toc);
else
    fprintf('The selected file type was not supported yet! email me to get support (zhoupc1988@gmail.com)\n');
    return;
end

data = matfile(nam_mat);
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames

% fprintf('\nThe data has been mapped to RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1, d2, numFrame, prod(Ysiz)*8/(2^30));
