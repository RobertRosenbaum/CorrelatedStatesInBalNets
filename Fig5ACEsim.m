
clear
close all

numtrials=50;

% Number of neurons in each pop
Ne=4000;
Ni=1000;
Nx=1000;
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

% Neuron parameters
Cm=1; 
gL=1/15;
EL=-72;
Vth=-50;
Vre=-75;
Vlb=-100;
DeltaT=1;
VT=-55;     



% Generate connectivity matrices
tic
J=[sparse(Jm(1,1)*binornd(1,P(1,1),Ne,Ne)) sparse(Jm(1,2)*binornd(1,P(1,2),Ne,Ni)); ...
   sparse(Jm(2,1)*binornd(1,P(2,1),Ni,Ne)) sparse(Jm(2,2)*binornd(1,P(2,2),Ni,Ni))];
Jx=[sparse(Jxm(1)*binornd(1,Px(1),Ne,Nx)); sparse(Jxm(2)*binornd(1,Px(2),Ni,Nx))];
tGen=toc
disp(sprintf('\nTime to generate connections: %.2f sec',tGen))

% Integer division function
IntDivide=@(n,k)(floor((n-1)/k)+1);

% Sim time and time step
T=1e5;
dt=.1;
time=dt:dt:T;
Nt=numel(time);

% Indices of neurons from which to record currents, voltages
Irecord=1:N;
numrecord=numel(Irecord);

% Time discretization of recordings. This can be coarser 
% than the dt for the Euler solver. If so, it will record
% the average current across each coarser time bin
dtRecord=250;
nBinsRecord=round(dtRecord/dt);
timeRecord=dtRecord:dtRecord:T;
Ntrec=numel(timeRecord);

% Random initial voltages
% only for first iteration
V0=rand(N,1)*(VT-Vre)+Vre;
V=V0;

for iii=1:numtrials
    
    [iii numtrials]

    % Make Poisson spike times for ffwd layer
    nspikeX=poissrnd(Nx*rx*T);
    st=rand(nspikeX,1)*T;
    sx=zeros(2,numel(st));
    sx(1,:)=sort(st);
    sx(2,:)=randi(Nx,1,numel(st)); % neuron indices


    % Maximum number of spikes for all neurons
    % in simulation. Make it 50Hz across all neurons
    % If there are more spikes, the simulation will
    % terminate
    maxns=ceil(.05*N*T);


    % Initialize stuff
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

    %CCC{iii}=SpikeCountCov(s,N,250,T-250,winsize);
    sss{iii}=s;
    
    AllRates0{iii}=hist(s(2,s(1,:)>250),1:N)/(T-250);
    
    mTotalInput0{iii}=mean(IxRec+IeRec+IiRec,2);

    
    %TTT{iii}=IeRec+IiRec+IxRec;
    XXX{iii}=IxRec;

end

SxSxth=eye(Nx)*rx*winsize;
CXXth=Jx*SxSxth*Jx';
[II,JJ]=meshgrid(1:N,1:N);II=II(:);JJ=JJ(:);
ccoef=zeros(numtrials,1);
errr=zeros(numtrials,1);
errrmg=zeros(numtrials,1);
LRerrr=zeros(numtrials,1);
LRerrrmg=zeros(numtrials,1);
mC=zeros(numtrials,1);
mCth=zeros(numtrials,1);
vC=zeros(numtrials,1);
vCth=zeros(numtrials,1);
mCthmg=zeros(numtrials,1);
vCthmg=zeros(numtrials,1);
mLHS=zeros(numtrials,1);
mRHS=zeros(numtrials,1);
vLHS=zeros(numtrials,1);
vRHS=zeros(numtrials,1);
mLHSmg=zeros(numtrials,1);
mRHSmg=zeros(numtrials,1);
vLHSmg=zeros(numtrials,1);
vRHSmg=zeros(numtrials,1);
for jjj=1:numtrials
    tic
    [jjj numtrials]
    CSS=0;
    CXX=0;
    AllRates=0;
    mTotalInput=0;
    for iii=1:jjj
        CSS=CSS+SpikeCountCov(sss{iii},N,250,T-250,winsize)/jjj;
        CXX=CXX+cov(XXX{iii}(:,2:end)'*winsize)/jjj;        
        AllRates=AllRates+AllRates0{iii}/jjj;
        mTotalInput=mTotalInput+mTotalInput0{iii}/jjj;
    end
       
    f = fittype('(a*(I-th)+b*(I-th).^2).*(I>th)','independent','I');
    efit=fit(mTotalInput(1:Ne),AllRates(1:Ne)',f,'StartPoint',[0 .001 -.1]);
    fe=@(I)((efit.a*(I-efit.th)+efit.b*(I-efit.th).^2).*(I>efit.th));

    ifit=fit(mTotalInput(Ne+1:N),AllRates(Ne+1:N)',f,'StartPoint',[0 .001 -.1]);
    fi=@(I)((ifit.a*(I-ifit.th)+ifit.b*(I-ifit.th).^2).*(I>ifit.th));

    feprime=@(I)((efit.a*I+2*efit.b*(I-efit.th)).*(I>efit.th));
    ge=feprime(mTotalInput(1:Ne));

    fiprime=@(I)((ifit.a*I+2*ifit.b*(I-ifit.th)).*(I>ifit.th));
    gi=fiprime(mTotalInput(Ne+1:N));

    g=[ge;gi];   
   
    G=diag(g);
    
    LHS=(eye(N)-G*J)*CSS*(eye(N)-G*J)';
    RHS=G*CXXth*G';    
    LRerrr(jjj)=mean(abs(LHS(II<JJ)-RHS(II<JJ)));    
    mLHS(jjj)=mean(LHS(II<JJ));
    mRHS(jjj)=mean(RHS(II<JJ));
    vLHS(jjj)=var(LHS(II<JJ));
    vRHS(jjj)=var(RHS(II<JJ));
    
    mg=mean(g);
    LHSmg=(eye(N)-mg*J)*CSS*(eye(N)-mg*J)';
    RHSmg=mg*CXXth*mg;
    LRerrrmg(jjj)=mean(abs(LHSmg(II<JJ)-RHSmg(II<JJ)));    
    mLHSmg(jjj)=mean(LHSmg(II<JJ));
    mRHSmg(jjj)=mean(RHSmg(II<JJ));
    vLHSmg(jjj)=var(LHSmg(II<JJ));
    vRHSmg(jjj)=var(RHSmg(II<JJ));
    
    invDJmg=inv(eye(N)/mg-full(J));
    Cthmg=invDJmg*CXXth*invDJmg';
    errrmg(jjj)=mean(abs(CSS(II<JJ)-Cthmg(II<JJ)));    
    
    mCthmg(jjj)=mean(Cthmg(II<JJ));
    vCthmg(jjj)=var(Cthmg(II<JJ));

    
    D=diag(1./g);
    invDJ=inv(D-full(J));
    Cth=invDJ*CXXth*invDJ';
    
    tmp=corrcoef(CSS(II>JJ),Cth(II>JJ));
    ccoef(jjj)=tmp(1,2);
    
    errr(jjj)=mean(abs(CSS(II<JJ)-Cth(II<JJ)));
    
    mC(jjj)=mean(CSS(II>JJ));
    mCth(jjj)=mean(Cth(II>JJ));
    
    vC(jjj)=var(CSS(II>JJ));
    vCth(jjj)=var(Cth(II>JJ));
    
    toc
    
end


%%%

GeoMeanRates=sqrt(AllRates'*AllRates);

nbins=5000;
[Chist,Cedges] = histcounts(CSS(II>JJ),nbins, 'Normalization','pdf');
Cedges = movmean(Cedges,2);
Cedges(1) = [];

nbins=5000;
[Cthhist,Cthedges] = histcounts(Cth(II>JJ),nbins, 'Normalization','pdf');
Cthedges = movmean(Cthedges,2);
Cthedges(1) = [];

nplot=300;
Cplot=CSS(II<JJ & JJ<nplot);
Cthplot=Cth(II<JJ & JJ<nplot);
LHSplot=LHS(II<JJ & JJ<nplot);
RHSplot=RHS(II<JJ & JJ<nplot);
LHSmgplot=LHSmg(II<JJ & JJ<nplot);
RHSmgplot=RHSmg(II<JJ & JJ<nplot);
Cthmgplot=Cthmg(II<JJ & JJ<nplot);
GeoMeanRatesplot=GeoMeanRates(II<JJ & JJ<nplot);

save Fig5ACEsim.mat N T LRerrrmg LRerrr mCthmg vCthmg mLHS mRHS mLHSmg mRHSmg vLHS vRHS vLHSmg vRHSmg vC vCth Cthmgplot RHSmgplot LHSmgplot RHSplot LHSplot mg mC mCth Cthhist Cthedges Chist Cedges T numtrials errr ccoef Cplot Cthplot GeoMeanRatesplot nplot;


