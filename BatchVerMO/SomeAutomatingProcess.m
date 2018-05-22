load('X:\elm\ProcessedCalciumData\AllRows\7030\ActuallyUsedInCNMFE\CaELM_row4444_b7030_20170612_145529.mat')
load('X:\shared\EmilyShijieShared\CaExtraction\7030\LogisticscnmfeBatchVer612.mat')
%%
neuron_full.options.d1=300;
neuron_full.options.d2=400;
[results, center, Cn, PNR, save_avi] = greedyROI_endoscope(double(Y), [], neuron_full.options,false,false);
%%
HY=results.HY;
pars=[];
localmax=results.ind_localmax;
for i=1:numel(results.ind_localmax)  
    y=HY(localmax(i),:);
    sn = GetSn(y);
    % estimate time constant
    switch neuron_full.options.deconv_options.type
        case 'ar1'
            par = estimate_time_constant(y, 1, sn);
            if length(par)~=1
                par=0;end
            pars(i)=par;
    end
end
%%
PNR0=results.PNR0; PNR_maxes=PNR0(localmax); PNR_ind=PNR_maxes>results.min_pnr;
test=[pars;PNR_ind'];
%%
corr=results.THRESH.Corr;
pnr=results.THRESH.PNR;
CorrOut=results.THRESH.CorrOut;
PNROut=results.THRESH.PNROut;
%%
PNR_all=reshape(PNR,300*400,1);
Cn=reshape(Cn,300*400,1);
%%
[coeff2,score2,latent2,tsquared2,explained2,mu2] = pca([PNR_all, Cn],'Centered',false);
t = score2(:,1)*coeff2(:,1)';
%%
figure
subplot(1,2,1)
scatter(score(:,1),score(:,2))
subplot(1,2,2)
scatter(Cn,PNR_all)
%%
t = score(:,1)*coeff1(:,1)' + repmat(mu,size(pixels,1),1);
%%
HY00=sort(HY0,2,'descend');
%MX = mean(HY00(:,floor(size(HY0,2).*0.02)),2);
MX = mean(HY00(:,8),2);
MEAN=mean(HY0,2);
STD=std(HY0,1,2);
deviation_ind=abs((MX-MEAN))./STD>4;
%deviation_sub=ind2sub([300,400],deviation_ind);
Corr_dev=Cn(deviation_ind);
PNR_dev=PNR(deviation_ind);
%%
nonzero=CorrOut>0;
plot(CorrOut(nonzero),PNROut(nonzero),'*r')
matr=[corr CorrOut(nonzero); pnr PNROut(nonzero)]';
[idx,C]=kmeans(matr,5);
plot(C(:,1),C(:,2),'*k','MarkerSize',15)
plot(matr(idx==1,1),matr(idx==1,2),'r.','MarkerSize',12)
hold on
plot(matr(idx==2,1),matr(idx==2,2),'b.','MarkerSize',12)
plot(matr(idx==3,1),matr(idx==3,2),'g.','MarkerSize',12)
plot(Corr_dev,PNR_dev,'*k','MarkerSize',8)
plot(matr(idx==4,1),matr(idx==4,2),'y.','MarkerSize',12)
plot(matr(idx==5,1),matr(idx==5,2),'y.','MarkerSize',12)