function neuron = fix_mat_files(root)

    %% Get path of all .mat files in root directory
    % get filenames to process
    files = dir([root, '/**/*-*_*_*_*.mat']);
    not = dir([root, '/**/*-*-*_*_*_*.mat']);
    [~,ia] = setdiff({files.name}, {not.name});
    files = files(ia);

    %% Open up each .mat file and save as a cleaned file (and format v6)
    for ff = 1:length(files)
        if ~exist([files(ff).folder, '/out.mat'], 'file')
            load([files(ff).folder, '/', files(ff).name]);
            neuron = struct(neuron);
            save([files(ff).folder, '/out.mat'],'-v6','-struct','neuron');
            clear neuron;
        end
    end

end