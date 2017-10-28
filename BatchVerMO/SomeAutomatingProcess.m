corr=results.THRESH.Corr;
pnr=results.THRESH.PNR;
CorrOut=results.THRESH.CorrOut;
PNROut=results.THRESH.PNROut;

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