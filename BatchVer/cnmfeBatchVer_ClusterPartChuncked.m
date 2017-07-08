%% C Run All below on cluster!
%% 0. Get cluster ready
if running_on_cluster % some procedures making cluster use robust
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end

%% 3 Merge similar neurons
%%% Merge similar neurons based on spatial AND temporal correlation
% C_all=cat(2,ACS.Cin);
% Amask_temp=cat(2,A0s{1:2});
% Amask_temp=bsxfun(@gt,Amask_temp,quantile(Amask_temp,0.3)); %only use central part for merging.

% C_temp=C_all(1:size(Amask_temp,2),:);
% [Amask_temp,C_temp,ACS] = mergeAC(Amask_temp,C_temp,ACS,merge_thr_2);
% 
% merge2start=1+size(A0s{1},2)+size(A0s{2},2);
% for i=3:length(samplelist)
%     Amask_temp=cat(2,Amask_temp,A0s{i});
%     Amask_temp=bsxfun(@gt,Amask_temp,quantile(Amask_temp,0.3));
%     
%     C_temp=[C_temp;C_all(merge2start:merge2start+size(A0s{i},2)-1,:)];
% 
%     [Amask_temp,C_temp,ACS] = mergeAC(Amask_temp,C_temp,ACS,merge_thr_2);
%     merge2start=merge2start+size(A0s{i},2);
% end

Amask_temp=cat(2,A0s{:});
Amask_temp=bsxfun(@gt,Amask_temp,quantile(Amask_temp,0.3)); %only use central part for merging.
%C_all=cat(2,ACS.Cin);

[Afinal,MC,newIDs,merged_ROIs] = mergeAC(Amask_temp,ACS,merge_thr_2);
% [size1,~]=cellfun(@size,newIDs);
% ACS(size1~=1)
save([outputdir 'commonAcnmfeBatchVer.mat'],'-v7.3')
%% 4 Collapsing A's
% As=cell(1,length(samplelist));
% STDs=cell(1,length(samplelist));
% parfor i=1:length(samplelist)
%     As{i}=ACS(i).Ain;
%     STDs{i}=ACS(i).STD;    
% end
% Afinal=ReducingA(As,STDs); % the format for cell input is designed for potential other applications.
% save([outputdir 'uniqueAcnmfeBatchVer.mat'],'-v7.3')
%% 4.5 Determine Afinal that will be used to extract C's in each file.

%%% Some processes making Afinal nicer, modified from Pengcheng Zhou's
%%% idea.
for i=1:size(Afinal,2)
    ai=Afinal(:,i);
    temp = full(ai>quantile(ai, 0.5, 1));
    ai(~temp(:)) = 0;
    Afinal(:,i)=ai;
end

% Just in case some all zero A's got passed to this stage.
nz_ind=any(Afinal);
Afinal=Afinal(:,nz_ind);
newIDs=newIDs(nz_ind);

save([outputdir 'AfinalcnmfeBatchVer.mat'],'-v7.3')