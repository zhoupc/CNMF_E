function Splitting(outputdir,comp_dir)
compiled=matfile(fullfile(comp_dir,'compiled.mat'));
load(fullfile(comp_dir, 'analysis.mat')); 

DIFF=[1,1; diff(compiled.FnumBnum)];
DIFF=sum(DIFF,2)>0;
Allframes=1:size(compiled.FnumBnum,1);
Boutframes=Allframes(DIFF)';
Boutframes=[Boutframes [Boutframes(2:end)-1;Allframes(end)]]; % 2 columns, start and end frame of each bout.
for i=1:size(Boutframes,1)
    Y=compiled.Y(:,:,Boutframes(i,1):Boutframes(i,2));
    Ysiz=Boutframes(i,2)-Boutframes(i,1)+1;
    
    Fnum=compiled.FnumBnum(Boutframes(i,1),1);
    sourcefilenam=dbase.SoundFiles(Fnum).name;
%    tmp = load(fullfile(comp_dir, dbase.SoundFiles(Fnum).name),'info');  
    filename=sprintf('Bout%.0fFile%s',i, sourcefilenam);
    save([outputdir filename '.mat'],'Y','Ysiz','-v7.3')
end
    % go to the original file where the bout is from to fill in info
%      Fnum=compiled.FnumBnum(Boutframes(i,1),1);
%      tmp = load(fullfile(dbasepath, dbase.SoundFiles(Fnum).name), 'Y', 'Ysiz', 'info');    
%     info= 
%     FnumBum=compiled.FnumBnum(Boutframes,:);
    
end
    