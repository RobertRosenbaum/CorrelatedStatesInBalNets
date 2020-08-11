
clear
close all

rand('seed',10);

T=5e4; % Length of sim in ms
Nmax = 2e4; % largest N
Nmin = 500; % smallest N
NN = 12; % number of points in N mesh
Nrange = round(logspace(log10(Nmin),log10(Nmax),NN)/NN)*NN;
numN = numel(Nrange);

numtrials=10;

tstart=tic;
mre=zeros(numtrials,numN);
mri=zeros(size(mre));
mX=zeros(size(mre));
mR=zeros(size(mre));
mE=zeros(size(mre));
mI=zeros(size(mre));
mCXX=zeros(size(mre));
mCXR=zeros(size(mre));
mCRR=zeros(size(mre));
mCII=zeros(size(mre));
vCXX=zeros(size(mre));
vCXR=zeros(size(mre));
vCRR=zeros(size(mre));
vCII=zeros(size(mre));
mCee=zeros(size(mre));
mCei=zeros(size(mre));
mCii=zeros(size(mre));
vCee=zeros(size(mre));
vCei=zeros(size(mre));
vCii=zeros(size(mre));
mSCorrEE=zeros(size(mre));
mSCorrEI=zeros(size(mre));
mSCorrII=zeros(size(mre));
vSCorrEE=zeros(size(mre));
vSCorrEI=zeros(size(mre));
vSCorrII=zeros(size(mre));
mSCorrAll=zeros(size(mre));
mSCovEE=zeros(size(mre));
mSCovEI=zeros(size(mre));
mSCovII=zeros(size(mre));
vSCovEE=zeros(size(mre));
vSCovEI=zeros(size(mre));
vSCovII=zeros(size(mre));
mSVarE=zeros(size(mre));
mSVarI=zeros(size(mre));
mSVar=zeros(size(mre));
mSCovAll=zeros(size(mre));
vSCovAll=zeros(size(mre));
mFFE=zeros(size(mre));
mFFI=zeros(size(mre));
mFFAll=zeros(size(mre));
for iiii=1:numtrials
    
    [iiii numtrials round(toc(tstart)/60)]
    %save VsNAsynchTemp.mat iii mSCorrEE;
    drawnow;
for iii=1:numN

    %[iiii numtrials iii numN round(toc(tstart)/60)]
    
    N=Nrange(iii);
    
    % Number of neurons in each pop
    Ne=round(.8*N);
    Ni=round(.2*N);
    Nx=round(.2*N);
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

    % Sim time and time step    
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

    % Window size over which to count spikes
    % and integrate currents for low freq
    % coherence
    winsize=dtRecord; 

    % Get mean input currents
    mX(iiii,iii)=mean(IxRec(:));
    mR(iiii,iii)=mean(IeRec(:)+IiRec(:));
    mE(iiii,iii)=mean(IeRec(:));
    mI(iiii,iii)=mX(iiii,iii)+mR(iiii,iii);
        
    % Vector external and recurrent input
    % currents integrated over blocks of 
    % size winsize.
    XRCurr = [IxRec; IeRec+IiRec]*winsize;
    XRCurr = XRCurr(:,3:end-1); % Get rid of first and last 250 ms
    C=cov(XRCurr'); % Compute covariance matrix between all pairs

    % Isolate each block to get XX, XR, RR covs
    [II,JJ]=meshgrid(1:2*numrecord,1:2*numrecord);II=II(:);JJ=JJ(:);
    CXX=C(II<numrecord & JJ<II);
    CXR=C(II<numrecord & JJ>numrecord);
    CRR=C(II>numrecord & JJ>II);
    clear C;

    % Comptue means
    mCXX(iiii,iii)=mean(CXX(:));
    mCXR(iiii,iii)=mean(CXR(:));
    mCRR(iiii,iii)=mean(CRR(:));
    vCXX(iiii,iii)=var(CXX(:));
    vCXR(iiii,iii)=var(CXR(:));
    vCRR(iiii,iii)=var(CRR(:));

    % Vector exc and inh input
    % currents integrated over blocks of 
    % size winsize.
    eiCurr = [IxRec+IeRec; IiRec]*winsize;
    eiCurr = eiCurr(:,3:end-1); % Get rid of first and last 250 ms
    C=cov(eiCurr'); % Compute covariance matrix between all pairs

    % Isolate each block to get XX, XR, RR covs
    [II,JJ]=meshgrid(1:2*numrecord,1:2*numrecord);II=II(:);JJ=JJ(:);
    Cee=C(II<numrecord & JJ<II);
    Cei=C(II<numrecord & JJ>numrecord);
    Cii=C(II>numrecord & JJ>II);
    clear C;

    % Comptue means
    mCee(iiii,iii)=mean(Cee(:));
    mCei(iiii,iii)=mean(Cei(:));
    mCii(iiii,iii)=mean(Cii(:));
    vCee(iiii,iii)=var(Cee(:));
    vCei(iiii,iii)=var(Cei(:));
    vCii(iiii,iii)=var(Cii(:));
        
    
    % Now get total input covariances
    TotalCurr=(IeRec+IiRec+IxRec)*winsize;
    TotalCurr = TotalCurr(:,3:end-1); % Get rid of first and last 250 ms
    C=cov(TotalCurr');
    [II,JJ]=meshgrid(1:numrecord,1:numrecord);II=II(:);JJ=JJ(:);
    CII=C(II>JJ);

    % Comptue mean
    mCII(iiii,iii)=mean(CII(:));
    vCII(iiii,iii)=var(CII(:));
    
    % Compute spike count correlations
    % for neurons with sufficiently high rates
    AllRates=hist(s(2,s(1,:)>500),1:N)/(T-500);
    Igood=find(AllRates>1/1000); % Find which neurons are >1Hz
    C=SpikeCountCorr(s,N,250,T-250,winsize,Igood);
    NeGood=find(Igood<=Ne,1,'last'); % How many are exc neurons

    mre(iiii,iii)=mean(AllRates(1:Ne));
    mri(iiii,iii)=mean(AllRates(Ne+1:N));
    
%     % Get spike count corrs over each sub-pop
    [II,JJ]=meshgrid(1:numel(Igood),1:numel(Igood));
    SCee=C(II<=NeGood & JJ<II & isfinite(C));
    SCei=C(II<=NeGood & JJ>NeGood & isfinite(C));
    SCii=C(II>NeGood & JJ>II & isfinite(C));
    SCall=C(JJ>II & isfinite(C));

    % Comptue means
    mSCorrAll(iiii,iii)=mean(SCall(:));
    mSCorrEE(iiii,iii)=mean(SCee(:));
    mSCorrEI(iiii,iii)=mean(SCei(:));
    mSCorrII(iiii,iii)=mean(SCii(:));
    vSCorrEE(iiii,iii)=var(SCee(:));
    vSCorrEI(iiii,iii)=var(SCei(:));
    vSCorrII(iiii,iii)=var(SCii(:));

    % Compute spike count covariances
    % for all neurons
    C=SpikeCountCov(s,N,500,T-250,winsize);

    % Spike count variances
    AllSCVars=diag(C);
    mSVarE(iiii,iii)=mean(AllSCVars(1:Ne));
    mSVarI(iiii,iii)=mean(AllSCVars(Ne+1:N));
    mSVar(iiii,iii)=mean(AllSCVars);
    
    % Get mean Fano factors
    minrate=1/1000;
    Igood=find(AllRates(1:Ne)>=minrate);
    mFFE(iiii,iii)=mean(AllSCVars(Igood)./(winsize*AllRates(Igood)'));
    Igood=find(AllRates(Ne+1:N)>=minrate);
    mFFI(iiii,iii)=mean(AllSCVars(Igood+Ne)./(winsize*AllRates(Igood+Ne)'));
    Igood=find(AllRates>=minrate);
    mFFAll(iiii,iii)=mean(AllSCVars(Igood)./(winsize*AllRates(Igood)'));
    
    % Get mean spike count covs over each sub-pop
    [II,JJ]=meshgrid(1:N,1:N);
    mSCovEE(iiii,iii)=mean(C(II<=Ne & JJ<II));
    mSCovEI(iiii,iii)=mean(C(II<=Ne & JJ>Ne));
    mSCovII(iiii,iii)=mean(C(II>Ne & JJ>II));
    mSCovAll(iiii,iii)=mean(C(II<JJ));
    vSCovEE(iiii,iii)=var(C(II<=Ne & JJ<II));
    vSCovEI(iiii,iii)=var(C(II<=Ne & JJ>Ne));
    vSCovII(iiii,iii)=var(C(II>Ne & JJ>II));
    vSCovAll(iiii,iii)=var(C(II<JJ));
    
end
end
t0=toc(tstart);

qe=Ne/N;qi=Ni/N;qx=Nx/N;
Wbar =P.*(Jm*sqrt(N)).*[qe qi; qe qi];
Wxbar=Px.*(Jxm*sqrt(N)).*[qx;qx];

SStemp=(rx/qx)*inv(Wbar)*(Wxbar*Wxbar')*inv(Wbar');
TSCovEE=(winsize./Nrange)*SStemp(1,1);
TSCovEI=(winsize./Nrange)*SStemp(2,1);
TSCovII=(winsize./Nrange)*SStemp(2,2);

rT=-inv(Wbar)*Wxbar*rx;


clear AllSCVars eiCurr counts edges edgest edgesi I CSD time C CII Cii Cei Cee CRR CXR CXX IeRec II JJ IiRec IxRec VRec J Jx s sx SCee SCei SCii TotalCurr XRCurr;
save Fig1EFGHI4AEsim.mat

