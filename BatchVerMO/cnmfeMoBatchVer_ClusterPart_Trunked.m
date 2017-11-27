%% C Run All below on cluster!

%% 0. Get cluster ready
if running_on_cluster % some procedures making cluster use robust
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end

%% 5 "massive" procedure: Extract A from each file
neuron_batch(length(samplelist)) = struct('ind_del',[],'rawsignal',[],'signal',[],'DeconvSpiketrain',[],'FileOrigin',[],'neuron',[],'C',[],'C_raw',[]);

parfor i= 1:length(samplelist)  
    mode='massive';
    nam=cell(1,2);
    nam{1}=fullfile(datadir,samplelist(i).name);
    nam{2}=File(i).Ysignal;
    [~,neuron_batch(i)]=demo_endoscope2(bg_neuron_ratio,[],with_dendrites,K,sframe,num2read,...
                                   nam,neuron_full,mode,[],neuron_batch(i),Afinal,...
                                   thresh_detecting_frames);
    neuron_batch(i).FileOrigin=filelist(i); % save origin(filelist)
end
neuron_batch = rmfield(neuron_batch,{'C','C_raw'});
fprintf('Massive extraction in each file done.');

%% 6 Save A*C

parfor i= 1:length(samplelist)
    for j=1:size(neuron_batch(i).neuron.A,2)
        jA=neuron_batch(i).neuron.A(:,j);
        jC=neuron_batch(i).neuron.C_raw(j,:);
        neuron_batch(i).rawsignal(j,:)=median(jA(jA>0)*jC);
    end        
    fprintf('neuron_batch %.0f extraction done\n', i);
end

%% 7 deconvolve signal
[~, ~, ~, ~,~,neuron_batch,boundary_raw]=PartTraces(neuron_batch);
Vars = {'neuron_batch';'boundary_raw'}; Vars=Vars';
fprintf('ALL extractions done.\n');
eval(sprintf('save %sCNMFE_BatchVer.mat %s -v7.3', outputdir, strjoin(Vars)));
fprintf('ALL data saved, check them out!');