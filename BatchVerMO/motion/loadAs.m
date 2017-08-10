%% Load data
DatadirForA=fullfile(cnmfedir,'*PartI_Afinal*');
Alist=dir(DatadirForA);
AnumInFolder=numel(Alist);

AsfromDaysPic=[];
AsfromDaysCell=cell(1,AnumInFolder);
AsfromDaysCell_central=cell(1,AnumInFolder);
%sizes=[];
for ia=1:AnumInFolder
    Anext=load(fullfile(cnmfedir,Alist(ia).name));

    Atemp=Anext.Afinal;
    AsfromDaysCell{ia}=Atemp;   % individual column is each A, for actual motion correction

    Atemp=centralA(Atemp);
    AsfromDaysCell_central{ia}=Atemp;

    k = size(Atemp,2);
    %sizes=[sizes k];
    
    AsfromDaysPic=cat(3,AsfromDaysPic,A2image(Atemp,300,400)); % whole picture, for registering, getting shifts        
end