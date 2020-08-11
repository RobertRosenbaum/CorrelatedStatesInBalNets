
clear
close all

lw=4;
fs=30;
eeclr=[.8 .2 .2];
iiclr=[.2 .2 .8];

nskip=10;

% Hyperpolarized-Depolarized mode
load hd.mat

v1=v1(1:nskip:end);
v2=v2(1:nskip:end);

% Total time is 1 sec
time=(1:numel(v1))*1000/numel(v1);

figure

% Plot V for first cell
plot(time,v1+37,'LineWidth',lw,'Color',eeclr)
hold on

% Plot V for second cell
plot(time,v2,'LineWidth',lw,'Color',iiclr)

% Plot scale bars
plot([100 300],[-60 -60],'k','LineWidth',lw)
plot([800 800],[-63 -43],'k','LineWidth',lw)

axis tight
axis off
temp=get(gcf,'Position');
temp(3:4)=[500*.9 400*.6];
set(gcf,'Position',temp)
set(gcf,'Color',[1 1 1])

%%% HH mode 
load hh.mat

v1=v1(1:nskip:end);
v2=v2(1:nskip:end);

figure
plot(time,v1,'LineWidth',lw,'Color',eeclr)
hold on
plot(time,v2+23,'LineWidth',lw,'Color',eeclr)
plot([100 300],[-60 -60],'k','LineWidth',lw)
plot([800 800],[-65 -45],'k','LineWidth',lw)
axis tight
axis off
temp=get(gcf,'Position');
temp(3:4)=[500*.9 400*.6];
set(gcf,'Position',temp)
set(gcf,'Color',[1 1 1])


%%% DD mode
load dd.mat

v1=v1(1:nskip:end);
v2=v2(1:nskip:end);

figure
plot(time,v1,'LineWidth',lw,'Color',iiclr)
hold on
plot(time,v2+25,'LineWidth',lw,'Color',iiclr)
plot([100 300],[-50 -50],'k','LineWidth',lw)
plot([800 800],[-55 -35],'k','LineWidth',lw)
axis tight
axis off
temp=get(gcf,'Position');
temp(3:4)=[500*.9 400*.6];
set(gcf,'Position',temp)
set(gcf,'Color',[1 1 1])

