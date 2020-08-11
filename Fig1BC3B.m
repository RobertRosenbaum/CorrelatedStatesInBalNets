
clear
close all

load Fig1BC3Bsim.mat

lw=4;
fs=30;

% Line colors
eeclr=[.8 .2 .2];
iiclr=[.2 .2 .8];

% Marker size
msize=10;

% Raster color
rasterclr=[0 0 0];

figure
plot(Splot(1,:)/1000,Splot(2,:),'.','Color',rasterclr,'MarkerSize',msize)
xlabel('time (s)')
ylabel('E neuron index')
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
set(gca,'XTick',[1 1.5 2])
set(gca,'YTick',[1 100 200])
axis([t1plot/1000 t2plot/1000 1 nsplot])
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)

figure
plot(Sxplot(1,:)/1000,Sxplot(2,:),'.','Color',rasterclr,'MarkerSize',msize)
xlabel('time (s)')
ylabel('X neuron index')
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
set(gca,'XTick',[1 1.5 2])
set(gca,'YTick',[1 100 200])
axis([t1plot/1000 t2plot/1000 1 nsplot])
temp=get(gcf,'Position');
temp(3:4)=[500 400];
set(gcf,'Position',temp)


figure
plot(plottime/1000,E1plot+X1plot ,'LineWidth',lw,'Color',eeclr)
hold on
plot(plottime/1000, I1plot,'LineWidth',lw,'Color',iiclr)
xlabel('time (s)')
axis([t1plot/1000 t2plot/1000 -Inf Inf])
set(gca,'XTick',[1 1.5 2])
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500*.9 400*.6];
set(gcf,'Position',temp)

figure
plot(plottime/1000,E2plot+X2plot ,'LineWidth',lw,'Color',eeclr)
hold on
plot(plottime/1000, I2plot,'LineWidth',lw,'Color',iiclr)
xlabel('time (s)')
axis([t1plot/1000 t2plot/1000 -Inf Inf])
set(gca,'XTick',[1 1.5 2])
box off
set(gca,'LineWidth',lw)
set(gca,'FontSize',fs)
temp=get(gcf,'Position');
temp(3:4)=[500*.9 400*.6];
set(gcf,'Position',temp)

