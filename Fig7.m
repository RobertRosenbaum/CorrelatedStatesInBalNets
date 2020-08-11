
clear
close all

lw=4;
fs=30;


load Fig7sim.mat;

EXclr=[.5 .1 .1];
Iclr=[.1 .1 .5];


figure
plot(plottime/1000,E1plot+X1plot ,'LineWidth',lw,'Color',EXclr)
hold on
plot(plottime/1000, I1plot,'LineWidth',lw,'Color',Iclr) 
xlabel('time (s)')
axis([t1plot/1000 t2plot/1000 -Inf Inf])
set(gca,'XTick',[1 1.5 2])
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500 400*.6];
set(gcf,'Position',temp)


figure
plot(plottime/1000,E2plot+X2plot ,'LineWidth',lw,'Color',EXclr)
hold on
plot(plottime/1000, I2plot,'LineWidth',lw,'Color',Iclr) 
xlabel('time (s)')
axis([t1plot/1000 t2plot/1000 -Inf Inf])
set(gca,'XTick',[1 1.5 2])
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500 400*.6];
set(gcf,'Position',temp)

figure
plot(plottime/1000,E3plot+X3plot ,'LineWidth',lw,'Color',EXclr)
hold on
plot(plottime/1000, I3plot,'LineWidth',lw,'Color',Iclr) 
xlabel('time (s)')
axis([t1plot/1000 t2plot/1000 -Inf Inf])
set(gca,'XTick',[1 1.5 2])
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500 400*.6];
set(gcf,'Position',temp)


