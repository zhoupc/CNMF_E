function selectedTracesForBatch(DataFolder, cnmfeFilePath, indSeqSort, ExtraMatrixToPlot, ExtraMatY, ExtraPlotYLabel)
% Input:
%   DataFolder should be where 'compiled.mat' and 'analysis.mat' is put.
%   cnmfeFilePath is where 

    color_palet = 1-[[1 0 0]; [1 .6 0]; [.7 .6 .4]; [.6 .8 .3]; [0 .6 .3]; [0 0 1]; [0 .6 1]; [0 .7 .7]; [.7 0 .7];  [.7 .4 1]]; 
    color_palet = color_palet([1:2:end 2:2:end],:); % scramble slightly
    
    if nargin<4; ExtraMatrixToPlot = []; end
    if nargin<5; ExtraMatY = 1:size(ExtraMatrixToPlot,1); end
    if nargin<6; ExtraPlotYLabel = ''; end
    
    %% Compute spectrogram and add time stamp
    % to make it faster, consider make song spectrogram ahead of time

    variableInfo = who('-file', fullfile(DataFolder, 'compiled.mat'));
    if and(~ismember('SongSpec', variableInfo),ismember('CompSoundSONG', variableInfo))
        display('need to compute song spectrogram... only need to do this once... could take about a minute')
        load(fullfile(DataFolder, 'compiled.mat'), 'CompSoundSONG', 'SOUNDfs')
        display('loaded song')
        tic; [SongSpec,SpecTime,SpecF] = spectrogramELM(CompSoundSONG,SOUNDfs,.005, 0); toc
        SongSpec = 10*log10(SongSpec); 
        display('computed spectrogram')
        save(fullfile(DataFolder, 'compiled.mat'), 'SongSpec', 'SpecTime', 'SpecF', ...
            '-append');
        display('saved spectrogram')
    elseif and(~ismember('SongSpec', variableInfo),~ismember('CompSoundSONG', variableInfo))
        fprintf('This file does not have sound information, could be sleep data.')
    end
    
    %%
    AddGcampDataTimestamps(DataFolder); % compute

    %% load data
    % (1) Recording data
    global istart w pressed patches h SpecIm ExtraIm Tit npat ...
        slines1 slines2 dbase FnumBnum Timestamps; 

    load(fullfile(DataFolder, 'compiled.mat'), 'Labels', 'segs', 'VIDEOfs',...
            'SOUNDfs', 'SongSpec','SpecTime','SpecF', 'FnumBnum', 'Timestamps');
    load(fullfile(DataFolder, 'analysis.mat')); 
    SongSpec = SongSpec + eps; 
    SongSpec(SongSpec(:)<prctile(SongSpec(:),75)) = prctile(SongSpec(:),75);   
    display('loaded data')
    
    % (2) Load extracted neurons activity from path, or take variable from function
    % input
    load(cnmfeFilePath, 'neuron_batch'); 
    PlotC=AllTraces(neuron_batch);
    if nargin<3; indSeqSort = 1:size(PlotC,1); end
    nColors = color_palet(mod(1:(size(PlotC,1)),size(color_palet,1))+1,:); 
    
    if ~exist('SOUNDfs')
        SOUNDfs = 40000; 
    end
    
    % for stripes between files
    borders = find((diff([0; FnumBnum(:,1)])~=0)|(diff([0; FnumBnum(:,1)])~=0))-1; 
    PlotC(:,borders(borders>0)) = nan; 
    
    clims = [0 prctile(PlotC(:),99)]; %[3 25]; 

    patches = {};
    istart = 1; 
    stepdur = 150; 
    npat = {}; 
    slines1 = {}; 
    slines2 = {}; 
%%    
    figure(1); clf; shg; 
    set(gcf, 'color', [1 1 1])
    InitializePlot();
    
    fpos = get(gcf, 'Position');
    set(gcf, 'WindowKeyPressFcn', @ButtonWasPressed, 'busyaction', 'cancel', 'interruptible', 'off')
    mTextBox = uicontrol('style','edit','callback', @TextboxCallback, 'Position', [fpos(3)-80 fpos(4)-40 60 20]);
    set(mTextBox,'String',num2str(istart));
end