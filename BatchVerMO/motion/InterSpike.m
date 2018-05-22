acf=ce;
for i=1:size(C,1)
    acf(i) = autocorr(C(i,:)) ;
end
    %%
    figure
    subplot(3,1,1)
    plot(lags{1}(3,:),r{1}(3,:));
    subplot(3,1,2)
    plot(lags{2}(6,:),r{2}(6,:));
    subplot(3,1,3)
    plot(lags{3}(5,:),r{3}(5,:));
    
    %%
    eachfilenum_cumsum=cumsum(eachdayfilenum);
    filenumsum = eachfilenum_cumsum(end);
    S_L=length(eachfilenum_cumsum);
    for i= 1:filenumsum % file
        k(i)=find((eachfilenum_cumsum>=i),1);
    end
    C_in_days=cell(k(end),1);
    for i =1:filenumsum
        C_in_days{k(i)}=[C_in_days{k(i)} ACS(i).Cin_raw];
    end
    %%
    r=cell(k(end),1);
    lags=r;
    for i =1:k(end)
        r_=zeros(size(C_in_days{i},1),(size(C_in_days{i},2)*2-1));
        lags_=r_;
        for j=1:size(C_in_days{i},1)
            [r_(j,:),lags_(j,:)] = xcorr(C_in_days{i}(j,:));
        end
        r{i}=r_;
        lags{i}=lags_;
    end
%% Auto-correlation
each_neuron_end=cumsum(cellfun(@(x) size(x,2),M1{1}));
each_neuron_start=[1,each_neuron_end(1:end-1)+1];
each_neuron_number=cellfun(@(x) size(x,2),M1{1});

C_raw=cell(numel(each_neuron_number),1);
for n=1:filenumsum % each day's temperal trace in each day's spatial trace
    d=find((eachfilenum_cumsum>=n),1);
    C_raw{d}=[C_raw{d} ACS(n).Cin_raw(each_neuron_start(d):each_neuron_end(d),:)];
end

r=zeros(sum(each_neuron_number),neuron_full.Fs*100*2+1);
lags=r;
for d=1:numel(each_neuron_number)% do autocorrelation
    for j=1:size(C_raw{d},1)
        k=each_neuron_start(d);
        [r(k+j-1,:),lags(k+j-1,:)] = xcorr(C_raw{d}(j,:),neuron_full.Fs*100);
    end
end
%% Cross correlation
R=zeros(sum(each_neuron_number),sum(each_neuron_number));
parfor j=1:sum(each_neuron_number)
    Rj=R(j,:);
    d=find((cumsum(each_neuron_number)>=j),1);
    if d==1
        for k=each_neuron_start(1):each_neuron_end(d+1)
            [R_,~] = xcorr(r(j,:),r(k,:),'coeff');Rj(k)=max(R_);
        end
    elseif d==numel(each_neuron_number)
        for k=each_neuron_start(d-1):each_neuron_end(d)
            [R_,~] = xcorr(r(j,:),r(k,:),'coeff');Rj(k)=max(R_);
        end
    else
        for k=each_neuron_start(d-1):each_neuron_end(d+1)
            [R_,~] = xcorr(r(j,:),r(k,:),'coeff');Rj(k)=max(R_);
        end
    end
    R(j,:)=Rj;
end

%% Print R
Rs_mean=zeros(1,numel(merged_ROIs));
Rs=cell(1,numel(merged_ROIs));
for n=1:numel(merged_ROIs)
    origins=merged_ROIs{n};
    Rs_=zeros(1,numel(origins)-1);
    for i=1:(numel(origins)-1)
        Rs_(i)=R(origins(i),origins(i+1));
    end
    Rs_mean(n)=mean(Rs_);
    Rs{n}=Rs_;
end