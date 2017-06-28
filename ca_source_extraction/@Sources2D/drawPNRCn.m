function drawPNRCn(obj,min_pnr,min_corr)
global Picname outputdir

THRESH=obj.P.THRESH;
numberOfNeuron=length(THRESH.Corr);

x = linspace(0,1);
y = min_pnr*min_corr./x;
figure
scatter(THRESH.CorrOut,THRESH.PNROut,20,'blue')
hold on
scatter(THRESH.Corr,THRESH.PNR,20,'red')
plot(x,y,'k')
plot([min_corr min_corr],[0 max(max(THRESH.PNROut),max(THRESH.PNR))],'k')
plot([0 1],[min_pnr min_pnr],'k')
MAXPNR=max(max(THRESH.PNROut),max(THRESH.PNR));
axis([0 1 0 MAXPNR])
text(THRESH.Corr,THRESH.PNR,cellstr(num2str((1:numberOfNeuron)'))','Color','black')
title(strcat('Current PNR and Corr cutoff: ',num2str(numberOfNeuron),'neurons'),'interpreter','none');
xlabel('Corr');
ylabel(strcat('PNR=',num2str(min_pnr)));
legend('Those outside neurons','Neuron seeds','Corr*PnrThresh','CorrThresh','PnrThresh')

fignam=[outputdir,Picname,strcat('PNR=',num2str(min_pnr)),'.png'];
saveas(gcf,fignam);
close(gcf);