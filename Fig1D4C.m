

clear
close all

lw=4;
fs=30;

RRclr=[.6 .2 .8];%[.8 .6 .1];
XXclr=[.2 .6 .2];
XRclr=[.8 .6 .1];%[.6 .6 .6];
IIclr=[0 0 0];

eeclr=[.8 .2 .2];
eiclr=[.6 .2 .8];
iiclr=[.2 .2 .8];

load Fig1D4Csim.mat

mSCorrAll

figure
plot(X_edges/winsize,X_hist, 'LineWidth',lw, 'Color',XXclr)
hold on
plot(R_edges/winsize,R_hist, 'LineWidth',lw, 'Color',RRclr)
plot(I_edges/winsize,I_hist, 'LineWidth',lw, 'Color',IIclr)
%axis([-9e4 9e4 0 max(ee_hist)*1.2])
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
axis([-802 802 0 max(ee_hist)*1.2])
xlabel('input current covariance')
ylabel('count')
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)
set(gca,'YTick',[])
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
