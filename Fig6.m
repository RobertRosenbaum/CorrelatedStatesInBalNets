
clear
close all

lw=4;
fs=30;

RRclr=[.6 .2 .8];%[.8 .6 .1];
XXclr=[.2 .6 .2];
XRclr=[.8 .6 .1];%[.6 .6 .6];
TTclr=[0 0 0];

eeclr=[.8 .2 .2];
eiclr=[.6 .2 .8];
iiclr=[.2 .2 .8];

%% Vs N stuff
load Fig6sim.mat

% 2x2 theoretical values 
qe=Ne/N;qi=Ni/N;qx=Nx/N;
Wbar =P.*(Jm*sqrt(N)).*[qe qi; qe qi];
Wxbar=Px.*(Jxm*sqrt(N)).*[qx;qx];
CXXth=2*(rx/qx)*(Wxbar*Wxbar')*winsize;
rT=-inv(Wbar)*Wxbar*rx;

% 4x4 theoretical total input cov
CTTth=.5*[CXXth -CXXth; -CXXth CXXth];


% Firing rates
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


% Mean input covs for same-population pairs
figure
errorbar(Nrange,mean(mCRRw),std(mCRRw)/sqrt(numtrials),'Color',RRclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(mCXXw),std(mCXXw)/sqrt(numtrials),'Color',XXclr,'LineWidth',lw)
errorbar(Nrange,mean(mCXRw),std(mCXRw)/sqrt(numtrials),'Color',XRclr,'LineWidth',lw)
errorbar(Nrange,mean(mCIIw),std(mCIIw)/sqrt(numtrials),'Color',TTclr,'LineWidth',lw)
plot(Nrange,0*Nrange+CTTth(1,1),'--','Color',[.65 .65 .65],'LineWidth',lw)
xlabel('network size (N)')
ylabel('mean input cov.')
set(gca,'LineWidth',lw)
box off
set(gca,'XTick',[1e4 2e4])
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)


% Mean input covs for opposite-population pairs
figure
errorbar(Nrange,mean(mCRRb),std(mCRRb)/sqrt(numtrials),'Color',RRclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(mCXXb),std(mCXXb)/sqrt(numtrials),'Color',XXclr,'LineWidth',lw)
errorbar(Nrange,mean(mCXRb),std(mCXRb)/sqrt(numtrials),'Color',XRclr,'LineWidth',lw)
errorbar(Nrange,mean(mCIIb),std(mCIIb)/sqrt(numtrials),'Color',TTclr,'LineWidth',lw)
plot(Nrange,0*Nrange+CTTth(3,1),'--','Color',[.65 .65 .65],'LineWidth',lw)
xlabel('network size (N)')
ylabel('mean input cov.')
set(gca,'LineWidth',lw)
box off
set(gca,'XTick',[1e4 2e4])
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)


% Spike count covs for same-population pairs
figure
errorbar(Nrange,mean(mSCovEEw),std(mSCovEEw)/sqrt(numtrials),'Color',eeclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(mSCovEIw),std(mSCovEI)/sqrt(numtrials),'Color',eiclr,'LineWidth',lw)
errorbar(Nrange,mean(mSCovIIw),std(mSCovIIw)/sqrt(numtrials),'Color',iiclr,'LineWidth',lw)
xlabel('network size (N)')
ylabel('mean spike count cov')
set(gca,'LineWidth',lw)
box off
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'FontSize',fs)


% Spike count covs for opposite-population pairs
figure
errorbar(Nrange,mean(mSCovEEb),std(mSCovEEb)/sqrt(numtrials),'Color',eeclr,'LineWidth',lw)
hold on
errorbar(Nrange,mean(mSCovEIb),std(mSCovEIb)/sqrt(numtrials),'Color',eiclr,'LineWidth',lw)
errorbar(Nrange,mean(mSCovIIb),std(mSCovIIb)/sqrt(numtrials),'Color',iiclr,'LineWidth',lw)
xlabel('network size (N)')
ylabel('mean spike count cov')
set(gca,'LineWidth',lw)
box off
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'FontSize',fs)


