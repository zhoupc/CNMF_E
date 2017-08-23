%% C-2 Run All below on cluster!

% load one of the logistics and overwrite some folders
cnmfefolder='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/';
load(fullfile(cnmfefolder,'LogisticscnmfeBatchVer20170713.mat'));

%% 0. Get cluster ready
if running_on_cluster % some procedures making cluster use robust
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end
load(fullfile(outputdir,'NiceAfinalcnmfeBatchVerMOTION.mat'))
filelist_fulllist=filelist_fulllist(1:25);

%% 5 "massive" procedure: Extract A from each file
neuron_batchMO(length(filelist_fulllist)) = struct('ind_del',[],'rawsignal',[],'signal',[],'FileOrigin',[],'neuron',[],'C',[],'C_raw',[]);

parfor i= 1:length(filelist_fulllist)
    mode='massive';
    nam=cell(1,2);
    nam{1}=fullfile(datadir,filelist_fulllist(i).name);
    nam{2}=File_fulllist(i).Ysignal;
    display([nam{1} num2str(size(nam{2},2))])
    k=find((eachfilenum_cumsum>=i),1);  
    display('line25')
    %[~,neuron_batchMO(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,min_pixel,bd,FS,SSub,TSub,bg_neuron_ratio,nam,mode,[],M3{k},neuron_batchMO(i),convolveType,merge_thr);
    [~,neuron_batchMO(i)]=demo_endoscope2(bg_neuron_ratio,merge_thr,with_dendrites,K,sframe,num2read,...
                                   nam,neuron_full,mode,[],neuron_batchMO(i),M3{k},...
                                   thresh_detecting_frames);
    display('line30')
    neuron_batchMO(i).neuron.C=neuron_batchMO(i).C;
    neuron_batchMO(i).neuron.C_raw=neuron_batchMO(i).C_raw;
    neuron_batchMO(i).neuron.A=M3{k};
    display('line36')
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