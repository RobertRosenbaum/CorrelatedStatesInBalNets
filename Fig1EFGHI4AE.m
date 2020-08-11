
clear
%close all

lw=4;
fs=30;

RRclr=[.6 .2 .8];
XXclr=[.2 .6 .2];
XRclr=[.8 .6 .1];
IIclr=[0 0 0];

eeclr=[.8 .2 .2];
eiclr=[.6 .2 .8];
iiclr=[.2 .2 .8];

load Fig1EFGHI4AEsimSmallerT.mat


% Mean inputs
figure
errorbar(Nrange,mean(mX),std(mX)/sqrt(numtrials),'Color',XXclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(mR),std(mR)/sqrt(numtrials),'Color',RRclr,'LineWidth',lw)
errorbar(Nrange,mean(mI),std(mI)/sqrt(numtrials),'Color',IIclr,'LineWidth',lw)
xlabel('network size (N)')
ylabel('mean current')
set(gca,'LineWidth',lw)
box off
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'FontSize',fs)

% Rates
figure
errorbar(Nrange,1000*mean(mre),1000*std(mre)/sqrt(numtrials),'Color',eeclr,'LineWidth',lw)
hold on
errorbar(Nrange,1000*mean(mri),1000*std(mri)/sqrt(numtrials),'Color',iiclr,'LineWidth',lw)
plot(Nrange,1000*rT(1)+0*Nrange,'--','Color',eeclr,'LineWidth',lw)
plot(Nrange,1000*rT(2)+0*Nrange,'--','Color',iiclr,'LineWidth',lw)
axis([0 Inf 0 16.1])
set(gca,'XTick',[1e4 2e4])
xlabel('network size (N)')
ylabel('mean rate (Hz)')
set(gca,'LineWidth',lw)
box off
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'FontSize',fs)

% Mean Input Covs
figure
errorbar(Nrange,mean(mCRR),std(mCRR)/sqrt(numtrials),'Color',RRclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(mCXX),std(mCXX)/sqrt(numtrials),'Color',XXclr,'LineWidth',lw)
errorbar(Nrange,mean(mCXR),std(mCXR)/sqrt(numtrials),'Color',XRclr,'LineWidth',lw)
errorbar(Nrange,mean(mCII),std(mCII)/sqrt(numtrials),'Color',IIclr,'LineWidth',lw)
xlabel('network size (N)')
ylabel('mean input cov.')
set(gca,'LineWidth',lw)
box off
set(gca,'FontSize',fs)
axis tight
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)


figure
errorbar(Nrange,mean(mCII),std(mCII)/sqrt(numtrials),'Color',IIclr,'LineWidth',lw)
hold on
plot(Nrange,mean(mCII(:,end))*Nrange(end)./Nrange,'--','Color',[.4 .4 .4],'LineWidth',lw)
set(gca,'LineWidth',lw)
box off
axis([-Inf Inf -Inf mean(mCII(:,1))*1.1])
set(gca,'FontSize',fs)
set(gca,'XScale','log')
set(gca,'XTick',[1e3 1e4])
set(gca,'YTick',[1e0 1e1])
set(gca,'YScale','log')
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
xlabel('network size (N)')
ylabel('mean input cov.')

% Mean spike count covs
figure
errorbar(Nrange,mean(mSCovEE),std(mSCovEE)/sqrt(numtrials),'Color',eeclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(mSCovEI),std(mSCovEI)/sqrt(numtrials),'Color',eiclr,'LineWidth',lw)
errorbar(Nrange,mean(mSCovII),std(mSCovII)/sqrt(numtrials),'Color',iiclr,'LineWidth',lw)
plot(Nrange,TSCovEE,'--','Color',eeclr,'LineWidth',lw)
plot(Nrange,TSCovEI,'--','Color',eiclr,'LineWidth',lw)
plot(Nrange,TSCovII,'--','Color',iiclr,'LineWidth',lw)
%legend('ee','','ei','','ii','')
xlabel('network size (N)')
ylabel('mean spike count cov')
set(gca,'LineWidth',lw)
set(gca,'XTick',[1e3 1e4])
box off
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'FontSize',fs)
set(gca,'XScale','log')
set(gca,'YScale','log')


% Mean Input Covs (exc and inh)
figure
errorbar(Nrange,mean(mCee),std(mCee)/sqrt(numtrials),'Color',eeclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(mCii),std(mCii)/sqrt(numtrials),'Color',iiclr,'LineWidth',lw)
errorbar(Nrange,mean(mCei),std(mCei)/sqrt(numtrials),'Color',eiclr,'LineWidth',lw)
errorbar(Nrange,mean(mCII),std(mCII)/sqrt(numtrials),'Color',IIclr,'LineWidth',lw)
xlabel('network size (N)')
ylabel('mean input cov.')
set(gca,'LineWidth',lw)
box off
axis tight
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)

% Var Input Covs (exc and inh)
figure
errorbar(Nrange,mean(vCee),std(vCee)/sqrt(numtrials),'Color',eeclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(vCii),std(vCii)/sqrt(numtrials),'Color',iiclr,'LineWidth',lw)
errorbar(Nrange,mean(vCei),std(vCei)/sqrt(numtrials),'Color',eiclr,'LineWidth',lw)
errorbar(Nrange,mean(vCII),std(vCII)/sqrt(numtrials),'Color',IIclr,'LineWidth',lw)
xlabel('network size (N)')
ylabel('var. input cov.')
set(gca,'LineWidth',lw)
box off
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'FontSize',fs)

% % Var input covs (X and R)
% figure
% errorbar(Nrange,mean(vCRR),std(vCRR)/sqrt(numtrials),'Color',RRclr,'LineWidth',lw)
% hold on
% errorbar(Nrange,mean(vCXX),std(vCXX)/sqrt(numtrials),'Color',XXclr,'LineWidth',lw)
% errorbar(Nrange,mean(vCXR),std(vCXR)/sqrt(numtrials),'Color',XRclr,'LineWidth',lw)
% errorbar(Nrange,mean(vCII),std(vCII)/sqrt(numtrials),'Color',IIclr,'LineWidth',lw)
% xlabel('network size (N)')
% ylabel('var. input cov.')
% %legend('RR','XX','XR','II')
% set(gca,'LineWidth',lw)
% box off
% temp=get(gcf,'Position');
% temp(3:4)=[500 400];
% set(gcf,'Position',temp)
% set(gca,'FontSize',fs)
