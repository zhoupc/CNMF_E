%% C-2 Run All below on cluster!

% load one of the logistics and overwrite some folders
cnmfefolder='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/';
load(fullfile(cnmfefolder,'LogisticscnmfeBatchVer20170713.mat'));

%% 0. Get cluster ready
if running_on_cluster % some procedures making cluster use robust
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end
load(fullfile(outputdir,'RoughAfinalcnmfeBatchVerMOTION.mat'))
filelist_fulllist=filelist_fulllist(1:50);

%% 4.5 Determine Afinal that will be used to extract C's in each file.

%%% Some processes making Afinal nicer, modified from Pengcheng Zhou's
%%% idea.
M3=cell(1,numel(M2));
for c=1:numel(M2)
    Afinal=M2{c};
    for i=1:size(Afinal,2)
        ai=Afinal(:,i);
        temp = full(ai>quantile(ai, 0.5, 1));
        ai(~temp(:)) = 0;
        Afinal(:,i)=ai;
    end
    % Just in case some all zero A's got passed to this stage.
    nz_ind=any(Afinal);
    Afinal=Afinal(:,nz_ind);
    newIDs{c}=newIDs{c}(nz_ind);
    Apicname=sprintf('Day%.0fAFinal',num2str(totaldays(c)));
    ColorAllNeurons(Afinal,d1,d2,Apicname,[outputdir, num2str(totaldays(c)), '/']);
    M3{c}=Afinal;
end
Vars = {'newIDs';'close_ind';'M2';'M3'}; Vars=Vars';
eval(sprintf('save %sAfinalcnmfeBatchVerMotion %s -v7.3', outputdir, strjoin(Vars)));
save([outputdir 'NiceAfinalcnmfeBatchVerMOTION.mat'],'-v7.3')
%% 5 "massive" procedure: Extract A from each file
neuron_batchMO(length(filelist_fulllist)) = struct('ind_del',[],'rawsignal',[],'signal',[],'FileOrigin',[],'neuron',[]);

parfor i= 1:length(filelist_fulllist)
    mode='massive';
    nam=cell(1,2);
    nam{1}=fullfile(datadir,filelist_fulllist(i).name);
    nam{2}=File_fulllist(i).Ysignal;
    display([nam{1} num2str(size(nam{2},2))])
    k=find((eachfilenum_cumsum>=i),1);      
    %[~,neuron_batchMO(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,min_pixel,bd,FS,SSub,TSub,bg_neuron_ratio,nam,mode,[],M3{k},neuron_batchMO(i),convolveType,merge_thr);
    [~,neuron_batchMO(i)]=demo_endoscope2(bg_neuron_ratio,merge_thr,with_dendrites,K,sframe,num2read,...
                                   nam,neuron_full,mode,[],neuron_batchMO(i),M3{k},...
                                   thresh_detecting_frames);
    neuron_batchMO(i).FileOrigin=filelist_fulllist(i); % save origin(filelist)
end
fprintf('Massive extraction done.');
save([outputdir 'MassivecnmfeBatchVerMotion.mat'],'-v7.3')

%% 6 Save results.

for i= 1:length(filelist_fulllist)
    for j=1:size(neuron_batchMO(i).neuron.A,2)
        jA=neuron_batchMO(i).neuron.A(:,j);
        jC=neuron_batchMO(i).neuron.C_raw(j,:);
        neuron_batchMO(i).rawsignal(j,:)=median(jA(jA>0)*jC);
    end        
    fprintf('neuron_batch %.0f extraction done\n', i);
end

%% 5.5 deconvolve signal
[~, ~, ~, ~,neuron_batchMO]=PartTraces(neuron_batchMO);

fprintf('ALL moBatchVer extractions done.\n');
eval(sprintf('save %sCNMFE_moBatchVer.mat %s -v7.3', outputdir, 'neuron_batchMO'));
fprintf('ALL moBatchVer data saved, check them out!');