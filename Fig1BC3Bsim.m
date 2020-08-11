
clear
close all

% Number of neurons in each pop
Ne=8000;
Ni=2000;
Nx=2000;
N=Ne+Ni;

% Mean synaptic weights
Jm=[25 -150; 112.5 -250]/sqrt(N);
Jxm=[180;135]/sqrt(N);

% Connection probabilities
P=[.1 .1; .1 .1];
Px=[.1; .1];

% Synaptic timescales
taue=8;
taui=4;
taux=10;

% Rate of external population
rx=10/1000;

% Correlation between external spike trains
% and timescale of correlation
taujitter=5;
c=0.1;

% Neuron parameters
Cm=1; 
gL=1/15;
EL=-72;
Vth=-50;
Vre=-75;
Vlb=-100;
DeltaT=1;
VT=-55;  

% Sim time and time step
T=2000;
dt=.1;
time=dt:dt:T;
Nt=numel(time);

% Indices of neurons from which to record currents, voltages
Irecord=1:4;
numrecord=numel(Irecord);

% Time discretization of recordings. This can be coarser 
% than the dt for the Euler solver. If so, it will record
% the average current across each coarser time bin
dtRecord=1;
nBinsRecord=round(dtRecord/dt);
timeRecord=dtRecord:dtRecord:T;
Ntrec=numel(timeRecord);

% Generate connectivity matrices
tic
J=[sparse(Jm(1,1)*binornd(1,P(1,1),Ne,Ne)) sparse(Jm(1,2)*binornd(1,P(1,2),Ne,Ni)); ...
   sparse(Jm(2,1)*binornd(1,P(2,1),Ni,Ne)) sparse(Jm(2,2)*binornd(1,P(2,2),Ni,Ni))];
Jx=[sparse(Jxm(1)*binornd(1,Px(1),Ne,Nx)); sparse(Jxm(2)*binornd(1,Px(2),Ni,Nx))];
tGen=toc
disp(sprintf('\nTime to generate connections: %.2f sec',tGen))

% Make Poisson spike times for ffwd layer
nspikeX=poissrnd(Nx*rx*T);
st=rand(nspikeX,1)*T;
sx=zeros(2,numel(st));
sx(1,:)=sort(st);
sx(2,:)=randi(Nx,1,numel(st)); % neuron indices
clear st;

% Random initial voltages
V0=rand(N,1)*(VT-Vre)+Vre;

% Maximum number of spikes for all neurons
% in simulation. Make it 50Hz across all neurons
% If there are more spikes, the simulation will
% terminate
maxns=ceil(.05*N*T);

% Integer division function
IntDivide=@(n,k)(floor((n-1)/k)+1);

% Initialize stuff
V=V0;
Ie=zeros(N,1);
Ii=zeros(N,1);
Ix=zeros(N,1);
IeRec=zeros(numrecord,Ntrec);
IiRec=zeros(numrecord,Ntrec);
IxRec=zeros(numrecord,Ntrec);
VRec=zeros(numrecord,Ntrec);
iXspike=1;
s=zeros(2,maxns);
nspike=0;
TooManySpikes=0;

% Euler solver
tic
for i=1:numel(time)

    % Store recorded variables
    ii=IntDivide(i,nBinsRecord); 
    IeRec(:,ii)=IeRec(:,ii)+Ie(Irecord);
    IiRec(:,ii)=IiRec(:,ii)+Ii(Irecord);
    IxRec(:,ii)=IxRec(:,ii)+Ix(Irecord);
    VRec(:,ii)=VRec(:,ii)+V(Irecord);    
    

    % Euler update to V
    V=V+(dt/Cm)*(Ie+Ii+Ix+gL*(EL-V)+gL*DeltaT*exp((V-VT)/DeltaT));
    V=max(V,Vlb);
    
    % Find which neurons spiked
    Ispike=find(V>=Vth);    


    % Euler update to synaptic currents
    Ie=Ie-dt*Ie/taue;
    Ii=Ii-dt*Ii/taui;
    Ix=Ix-dt*Ix/taux;    
    
    % If there are spikes
    if(~isempty(Ispike))

        % Store spike times and neuron indices
        if(nspike+numel(Ispike)<=maxns)
            s(1,nspike+1:nspike+numel(Ispike))=time(i);
            s(2,nspike+1:nspike+numel(Ispike))=Ispike;
        else
            TooManySpikes=1;
            break;
        end
                
        % Reset mem pot.
        V(Ispike)=Vre;
        
        % Update synaptic currents
        Ie=Ie+sum(J(:,Ispike(Ispike<=Ne)),2)/taue;    
        Ii=Ii+sum(J(:,Ispike(Ispike>Ne)),2)/taui;            
              
        % Update cumulative number of spikes
        nspike=nspike+numel(Ispike);
    end
    
    % Propogate ffwd spikes
    while(sx(1,iXspike)<=time(i) && iXspike<nspikeX)
        jpre=sx(2,iXspike);
        Ix=Ix+Jx(:,jpre)/taux;
        iXspike=iXspike+1;
    end
        
end
if(TooManySpikes)
    warning('\nMax number spikes exceeded, simulation terminated at time t=%1.1f.\n',dt*i)
end
IeRec=IeRec/nBinsRecord; % Normalize recorded variables by # bins
IiRec=IiRec/nBinsRecord;
IxRec=IxRec/nBinsRecord;
VRec=VRec/nBinsRecord;
s=s(:,1:nspike); % Get rid of padding in s
tSim=toc;
disp(sprintf('\nTime for simulation: %.2f min',tSim/60))


nsplot=200; % Num of neurons to plot for raster
t1plot=1000; % Plot start time
t2plot=2000; % Plot end time
plottime=t1plot:dtRecord:t2plot;

% Store data to plot

E1plot=IeRec(1,timeRecord>=t1plot & timeRecord<=t2plot);
E2plot=IeRec(2,timeRecord>=t1plot & timeRecord<=t2plot);
I1plot=IiRec(1,timeRecord>=t1plot & timeRecord<=t2plot);
I2plot=IiRec(2,timeRecord>=t1plot & timeRecord<=t2plot);
X1plot=IxRec(1,timeRecord>=t1plot & timeRecord<=t2plot);
X2plot=IxRec(2,timeRecord>=t1plot & timeRecord<=t2plot);

Splot=s(:,s(1,:)>=t1plot & s(1,:)<=t2plot & s(2,:)<nsplot);
Sxplot=sx(:,sx(1,:)>=t1plot & sx(1,:)<=t2plot & sx(2,:)<nsplot);

% Save data to plot
save Fig1BC3Bsim.mat nsplot plottime E1plot E2plot I1plot I2plot X1plot X2plot Splot Sxplot t1plot t2plot;



