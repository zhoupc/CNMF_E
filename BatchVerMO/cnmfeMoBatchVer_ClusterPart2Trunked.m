%% C-2 Run All below on cluster!

%% 0. Get cluster ready
if running_on_cluster % some procedures making cluster use robust
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end
load(fullfile(outputdir,'RoughAfinalcnmfeBatchVerMOTION.mat'))

%% 4 Determine Afinal that will be used to extract C's in each file.

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
    if c==1
        newIDs=newIDs(nz_ind);
    end
    Apicname=sprintf('Day%.0fAFinal',num2str(totaldays(c)));
    ColorAllNeurons(Afinal,d1,d2,Apicname,[outputdir, num2str(totaldays(c)), '/']);
    M3{c}=Afinal;
end
Vars = {'newIDs';'close_ind';'M2';'M3'}; Vars=Vars';
eval(sprintf('save %sAfinalcnmfeBatchVerMotion %s -v7.3', outputdir, strjoin(Vars)));
save([outputdir 'NiceAfinalcnmfeBatchVerMOTION.mat'],'-v7.3')

%% 5 "massive" procedure: Extract A from each file
neuron_batchMO(length(filelist_fulllist)) = struct('ind_del',[],'rawsignal',[],'signal',[],'FileOrigin',[],'neuron',[],'C',[],'C_raw',[]);

parfor i= 1:length(filelist_fulllist)
    mode='massive';
    nam=cell(1,2);
    nam{1}=fullfile(datadir,filelist_fulllist(i).name);
    nam{2}=File_fulllist(i).Ysignal;
    k=find((eachfilenum_cumsum>=i),1);  
    [~,neuron_batchMO(i)]=demo_endoscope2(bg_neuron_ratio,[],with_dendrites,K,sframe,num2read,...
                                   nam,neuron_full,mode,[],neuron_batchMO(i),M3{k},...
                                   thresh_detecting_frames);

    neuron_batchMO(i).FileOrigin=filelist_fulllist(i); % save origin(filelist)
end
neuron_batchMO = rmfield(neuron_batchMO,{'C','C_raw'});
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

%% 7 deconvolve signal
[~, ~, ~, ~,~,neuron_batchMO,boundary_raw]=PartTraces(neuron_batchMO);
Vars = {'neuron_batchMO';'boundary_raw'}; Vars=Vars';

fprintf('ALL moBatchVer extractions done.\n');
eval(sprintf('save %sCNMFE_moBatchVer.mat %s -v7.3', outputdir, strjoin(Vars)));
fprintf('ALL moBatchVer data saved, check them out!');