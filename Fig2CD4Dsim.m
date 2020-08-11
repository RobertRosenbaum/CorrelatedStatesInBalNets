
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
T=5e4;
dt=.1;
time=dt:dt:T;
Nt=numel(time);

% Indices of neurons from which to record currents, voltages
Irecord=1:Ne;
numrecord=numel(Irecord);

% Time discretization of recordings. This can be coarser 
% than the dt for the Euler solver. If so, it will record
% the average current across each coarser time bin
dtRecord=250;
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


% Make correlated Poisson spike times for ffwd layer
tic
rm=rx/c; % Firing rate of mother process
nstm=poissrnd(rm*T); % Number of mother spikes
stm=rand(nstm,1)*T; % spike times of mother process    
maxnsx=T*rx*Nx*1.2; % Max num spikes
sx=zeros(2,maxnsx);
ns=0;
for j=1:Nx  % For each ffwd spike train
    ns0=binornd(nstm,c); % Number of spikes for this spike train
    st=randsample(stm,ns0); % Sample spike times randomly
    st=st+taujitter*randn(size(st)); % jitter spike times
    st=st(st>0 & st<T); % Get rid of out-of-bounds times
    ns0=numel(st); % Re-compute spike count
    sx(1,ns+1:ns+ns0)=st; % Set the spike times and indices        
    sx(2,ns+1:ns+ns0)=j;
    ns=ns+ns0;
end
% Get rid of padded zeros
sx = sx(:,sx(1,:)>0);
% Sort by spike time
[~,I] = sort(sx(1,:));
sx = sx(:,I);
nspikeX=size(sx,2);
tGenSx=toc;
disp(sprintf('\nTime to generate X spike times: %.2f sec',tGenSx))


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

% Window size over which to count spikes
% and integrate currents for low freq
% coherence
winsize=dtRecord; 

% Vector external and recurrent input
% currents integrated over blocks of 
% size winsize.
XRCurr = [IxRec; IeRec+IiRec]*winsize;
XRCurr = XRCurr(:,2:end-1); % Get rid of first and last 250 ms

C=cov(XRCurr'); % Compute covariance matrix between all pairs

% Isolate each block to get XX, XR, RR covs
[II,JJ]=meshgrid(1:2*numrecord,1:2*numrecord);II=II(:);JJ=JJ(:);
CXX=C(II<numrecord & JJ<II);
CXR=C(II<numrecord & JJ>numrecord);
CRR=C(II>numrecord & JJ>II);
clear C;

% Comptue means
mCXX=mean(CXX(:));
mCXR=mean(CXR(:));
mCRR=mean(CRR(:));

% Now get total input covariances
TotalCurr=(IeRec+IiRec+IxRec)*winsize;
TotalCurr = TotalCurr(:,2:end-1); % Get rid of first and last 250 ms
C=cov(TotalCurr');
[II,JJ]=meshgrid(1:numrecord,1:numrecord);II=II(:);JJ=JJ(:);
CTT=C(II>JJ);

% Comptue mean
mCTT=mean(CTT(:));

% Vector external and recurrent input
% currents integrated over blocks of 
% size winsize.
EICurr = [IeRec+IxRec; IiRec]*winsize;
EICurr = EICurr(:,2:end-1); % Get rid of first and last 250 ms
C=cov(EICurr'); % Compute covariance matrix between all pairs

% Isolate each block to get XX, XR, RR covs
[II,JJ]=meshgrid(1:2*numrecord,1:2*numrecord);II=II(:);JJ=JJ(:);
Cee=C(II<numrecord & JJ<II);
Cei=C(II<numrecord & JJ>numrecord);
Cii=C(II>numrecord & JJ>II);
clear C;

% Comptue means
mCee=mean(Cee(:));
mCei=mean(Cei(:));
mCii=mean(Cii(:));



% Compute histograms of current covariances
nbins = 1e2;

% Histogram current covariances
[XX_hist,XX_edges] = histcounts(CXX,nbins, 'Normalization','pdf');
[XR_hist,XR_edges] = histcounts(CXR,nbins, 'Normalization','pdf');
[RR_hist,RR_edges] = histcounts(CRR,nbins, 'Normalization','pdf');
[II_hist,II_edges] = histcounts(CTT,nbins, 'Normalization','pdf');
[ee_hist,ee_edges] = histcounts(Cee,nbins, 'Normalization','pdf');
[ei_hist,ei_edges] = histcounts(Cei,nbins, 'Normalization','pdf');
[ii_hist,ii_edges] = histcounts(Cii,nbins, 'Normalization','pdf');
[X_hist,X_edges] = histcounts(mean(XRCurr(1:numrecord,:),2),nbins, 'Normalization','pdf');
[R_hist,R_edges] = histcounts(mean(XRCurr(numrecord+1:end,:),2),nbins, 'Normalization','pdf');
[I_hist,I_edges] = histcounts(mean(XRCurr(1:numrecord,:)+XRCurr(numrecord+1:end,:),2),nbins, 'Normalization','pdf');
[e_hist,e_edges] = histcounts(mean(EICurr(1:numrecord,:),2),nbins, 'Normalization','pdf');
[i_hist,i_edges] = histcounts(mean(EICurr(numrecord+1:end,:),2),nbins, 'Normalization','pdf');

% Center the bins
XX_edges = movmean(XX_edges,2);
XX_edges(1) = [];
XR_edges = movmean(XR_edges,2);
XR_edges(1) = [];
RR_edges = movmean(RR_edges,2);
RR_edges(1) = [];
II_edges = movmean(II_edges,2);
II_edges(1) = [];
ee_edges = movmean(ee_edges,2);
ee_edges(1) = [];
ei_edges = movmean(ei_edges,2);
ei_edges(1) = [];
ii_edges = movmean(ii_edges,2);
ii_edges(1) = [];
X_edges = movmean(X_edges,2);
X_edges(1) = [];
R_edges = movmean(R_edges,2);
R_edges(1) = [];
I_edges = movmean(I_edges,2);
I_edges(1) = [];
e_edges = movmean(e_edges,2);
e_edges(1) = [];
i_edges = movmean(i_edges,2);
i_edges(1) = [];


% Compute spike count correlations
% for neurons with sufficiently high rates
AllRates=hist(s(2,s(1,:)>250),1:N)/(T-250);
Igood=find(AllRates>1/1000); % Find which neurons are >1Hz
C=SpikeCountCorr(s,N,250,T-250,winsize,Igood);
NeGood=find(Igood<=Ne,1,'last'); % How many are exc neurons

% Get spike count corrs over each sub-pop
[II,JJ]=meshgrid(1:numel(Igood),1:numel(Igood));
SCall=C(JJ<II & isfinite(C));
SCee=C(II<=NeGood & JJ<II & isfinite(C));
SCei=C(II<=NeGood & JJ>NeGood & isfinite(C));
SCii=C(II>NeGood & JJ>II & isfinite(C));

% Compute histograms of spike count corrs
[SCEE_hist,SCEE_edges] = histcounts(SCee,nbins, 'Normalization','pdf');
[SCEI_hist,SCEI_edges] = histcounts(SCei,nbins, 'Normalization','pdf');
[SCTT_hist,SCTT_edges] = histcounts(SCii,nbins, 'Normalization','pdf');

% Center the bins
SCEE_edges = movmean(SCEE_edges,2);
SCEE_edges(1) = [];
SCEI_edges = movmean(SCEI_edges,2);
SCEI_edges(1) = [];
SCTT_edges = movmean(SCTT_edges,2);
SCTT_edges(1) = [];

% Comptue means
mSCorrAll=mean(SCall(:))
mSCorrEE=mean(SCee(:));
mSCorrEI=mean(SCei(:));
mSCorrII=mean(SCii(:));

% Compute spike count covariances
% for all neurons
C=SpikeCountCov(s,N,250,T-250,winsize);

% Get mean spike count covs over each sub-pop
[II,JJ]=meshgrid(1:N,1:N);
mSCovEE=mean(C(II<=Ne & JJ<II));
mSCovEI=mean(C(II<=Ne & JJ>Ne));
mSCovII=mean(C(II>Ne & JJ>II));

% Compute Spike Train CSDs
nCSD=150;
nreps=12;
winrad=2000;
dtsc=4;
tic
mCSDEE=0;mCSDEI=0;mCSDII=0;
for jjj=1:nreps

    [jjj nreps]
    
    % Edges for histogram
    edgest=0:dtsc:(T);
    edgesi=(1:N+1)-.01;
    edges={edgest,edgesi};
    % Get 2D histogram of spike indices and times
    counts=hist3(s','Edges',edges);
    % Get rid of edges, which are 0
    counts=counts(1:end-1,1:end-1);

    % Burn in one winrad at beginning and end
    counts=counts(2:end-1,:);

    % EE
    CSDInds1=(1:nCSD)+(jjj-1)*nCSD;
    CSDInds2=nCSD+CSDInds1;
    [CSD,F]=GetCrossSpec(counts(:,CSDInds1)/dtsc,counts(:,CSDInds2)/dtsc,winrad,dtsc);
    mCSDEE=mCSDEE+mean(mean(CSD,2),3)/nreps;
    disp('EECSD')

    % EI
    CSDInds1=(1:nCSD)+(jjj-1)*nCSD;
    CSDInds2=Ne+CSDInds1;
    [CSD,F]=GetCrossSpec(counts(:,CSDInds1)/dtsc,counts(:,CSDInds2)/dtsc,winrad,dtsc);
    mCSDEI=mCSDEI+mean(mean(CSD,2),3)/nreps;
    disp('EICSD')

    % II
    CSDInds1=Ne+(1:nCSD)+(jjj-1)*nCSD;
    CSDInds2=nCSD+CSDInds1;
    [CSD,F]=GetCrossSpec(counts(:,CSDInds1)/dtsc,counts(:,CSDInds2)/dtsc,winrad,dtsc);
    mCSDII=mCSDII+mean(mean(CSD,2),3)/nreps;
    disp('IICSD')
end

tCSD=toc;
disp(sprintf('\nTime for CSD calc: %.2f min',tCSD/60))

% Compute theoretical CSD matrix
qe=Ne/N;qi=Ni/N;qx=Nx/N;
Wbar =P.*(Jm*sqrt(N)).*[qe qi; qe qi];
Wxbar=Px.*(Jxm*sqrt(N)).*[qx;qx];

etaetilde=1./(1+2*pi*sqrt(-1)*taue*F);
etaitilde=1./(1+2*pi*sqrt(-1)*taui*F);
etaxtilde=1./(1+2*pi*sqrt(-1)*taux*F);
SxSx=c*rx*exp(-4*pi^2*F.^2*taujitter^2);

wee=Wbar(1,1)*etaetilde;
wei=Wbar(1,2)*etaitilde;
wie=Wbar(2,1)*etaetilde;
wii=Wbar(2,2)*etaitilde;
wex=Wxbar(1)*etaxtilde;
wix=Wxbar(2)*etaxtilde;

SSee=zeros(size(F));
SSei=zeros(size(F));
SSii=zeros(size(F));
for i=1:numel(F)
    Wtemp=[wee(i) wei(i);wie(i) wii(i)];
    Wxtemp=[wex(i); wix(i)];
    SStemp=inv(Wtemp)*(Wxtemp*SxSx(i)*Wxtemp')*inv(Wtemp');
    SSee(i)=SStemp(1,1);
    SSei(i)=abs(SStemp(1,2));
    SSii(i)=SStemp(2,2);
end

clear EICurr SCall Cee Cei Cii counts edges edgest edgesi I CSD time C CTT CRR CXR CXX IeRec II JJ IiRec IxRec VRec J Jx s sx SCee SCei SCii TotalCurr XRCurr;
save Fig2CD4Dsim.mat;

