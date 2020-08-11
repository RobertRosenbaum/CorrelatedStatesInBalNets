

clear
close all

lw=4;
fs=30;

RRclr=[.6 .2 .8];
XXclr=[.2 .6 .2];
XRclr=[.8 .6 .1];
IIclr=[0 0 0];

eeclr=[.8 .2 .2];
eiclr=[.6 .2 .8];
iiclr=[.2 .2 .8];

load Fig2CD4Dsim.mat

mSCorrAll

% Plot spike train CSDs
figure
plot(1000*F,1000*real(mCSDEE), 'LineWidth',lw, 'Color',eeclr+[.2 .2 .2])
hold on
plot(1000*F,1000*real(mCSDEI), 'LineWidth',lw, 'Color',eiclr+[.2 .2 .2])
plot(1000*F,1000*real(mCSDII), 'LineWidth',lw, 'Color',iiclr+[.2 .2 .2])
plot(1000*F,1000*real(SSee),'--', 'LineWidth',lw, 'Color',eeclr-[.1 .1 .1])
plot(1000*F,1000*real(SSei),'--', 'LineWidth',lw, 'Color',eiclr-[.1 .1 .1])
plot(1000*F,1000*real(SSii),'--', 'LineWidth',lw, 'Color',iiclr-[.1 .1 .1])
axis([0 100 -Inf Inf])
ylabel('mean CSD (Hz)')
xlabel('frequency (Hz)')
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)


figure
plot(X_edges/winsize,X_hist, 'LineWidth',lw, 'Color',XXclr)
hold on
plot(R_edges/winsize,R_hist, 'LineWidth',lw, 'Color',RRclr)
plot(I_edges/winsize,I_hist, 'LineWidth',lw, 'Color',IIclr)
axis tight
xlabel('time-averaged input current')
ylabel('count')
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'YTick',[])
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)


figure
plot( ee_edges,ee_hist, 'LineWidth',lw, 'Color',eeclr)
hold on
plot( ei_edges,ei_hist, 'LineWidth',lw, 'Color',eiclr)
plot( ii_edges,ii_hist, 'LineWidth',lw, 'Color',iiclr)
plot( II_edges,II_hist, 'LineWidth',lw, 'Color',IIclr)
%legend('XX','XR','RR','II (total)')
%set(gca,'xlim',[-802 802])
axis([-9e4 9e4 0 max(ee_hist)*1.2])
xlabel('input current covariance')
ylabel('count')
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'YTick',[])
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)



