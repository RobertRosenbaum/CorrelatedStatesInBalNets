

clear
close all

lw=4;
fs=30;


load Fig5BDFsim.mat

npairplot=1e4;
Inds=randperm(numel(Cplot),npairplot);
msize=10;
lwscatter=3;

figure
scatter(Cplot(Inds),Cthplot(Inds),msize,1000*GeoMeanRatesplot(Inds),'LineWidth',lwscatter)
hold on
colormap jet
colorbar
hold on
axis tight
axis(axis+[-1 1 -1 1])
plot([min(axis) max(axis)],[min(axis) max(axis)],'k--','LineWidth',lw)
box off
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
xlabel('empirical spike count cov.')
ylabel('theoretical')
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)


figure
scatter(LHSplot(Inds),RHSplot(Inds),msize,1000*GeoMeanRatesplot(Inds),'LineWidth',lwscatter)
colormap jet
colorbar
%plot(LHSplot(Inds),RHSplot(Inds),'.','MarkerSize',msize,'Color',[.5 .5 .5])%,'LineWidth',lwscatter)
hold on
axis tight
axis(axis+[-1.5 3.5 -1.5 1.5])
plot([min(axis) max(axis)],[min(axis) max(axis)],'k--','LineWidth',lw)
box off
set(gca,'XTick',[0 20 40])
set(gca,'YTick',[0 20 40])
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
xlabel('LHS')
ylabel('RHS')
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)



figure
scatter(LHSmgplot(Inds),RHSmgplot(Inds),msize,1000*GeoMeanRatesplot(Inds),'LineWidth',lwscatter)
colormap jet
colorbar
hold on
axis tight
axis(axis+[-1.5 1.5 -1.5 1.5])
plot([min(axis) max(axis)],[min(axis) max(axis)],'k--','LineWidth',lw)
box off
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
xlabel('LHS')
ylabel('RHS')
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)

figure
normfactor=1;
simtime=(1:numtrials)*T/1000;
plot(simtime(isfinite(LRerrr)),LRerrr(isfinite(LRerrr))./normfactor,'Color',[.4 .4 .4],'LineWidth',lw)
axis([0 Inf -Inf Inf])
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
box off
xlabel('total sim. time (s)')
ylabel('mean abs. error')
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
