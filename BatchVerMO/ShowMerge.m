picDir='/Users/gushijie/Documents/Fee/BatchedVerDebug/6922/';
picdir=fullfile(picDir,'*CaELM*');
piclist=dir(picdir);
C=imread(fullfile(picDir,piclist(1).name));
for i=2:numel(piclist) 
    pic=imread(fullfile(picDir,piclist(i).name));
    C= imfuse(C,pic,'blend');
end



