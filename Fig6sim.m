
clear
close all

rand('seed',10);

T=5e4; % Length of sim in ms
Nmax = 2e4; % largest N
Nmin = 500; % smallest N
NN = 12; % number of points in N mesh
Nrange = round(logspace(log10(Nmin),log10(Nmax),NN)/(NN*10))*NN*10;
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

mCXXw=zeros(size(mre));
mCXRw=zeros(size(mre));
mCRRw=zeros(size(mre));
mCIIw=zeros(size(mre));
vCXXw=zeros(size(mre));
vCXRw=zeros(size(mre));
vCRRw=zeros(size(mre));
vCIIw=zeros(size(mre));
mCXXb=zeros(size(mre));
mCXRb=zeros(size(mre));
mCRRb=zeros(size(mre));
mCIIb=zeros(size(mre));
vCXXb=zeros(size(mre));
vCXRb=zeros(size(mre));
vCRRb=zeros(size(mre));
vCIIb=zeros(size(mre));

mSCorrEE=zeros(size(mre));
mSCorrEI=zeros(size(mre));
mSCorrII=zeros(size(mre));
vSCorrEE=zeros(size(mre));
vSCorrEI=zeros(size(mre));
vSCorrII=zeros(size(mre));
mSCovEE=zeros(size(mre));
mSCovEI=zeros(size(mre));
mSCovII=zeros(size(mre));
vSCovEE=zeros(size(mre));
vSCovEI=zeros(size(mre));
vSCovII=zeros(size(mre));

mSCorrEEw=zeros(size(mre));
mSCorrEIw=zeros(size(mre));
mSCorrIIw=zeros(size(mre));
vSCorrEEw=zeros(size(mre));
vSCorrEIw=zeros(size(mre));
vSCorrIIw=zeros(size(mre));
mSCovEEw=zeros(size(mre));
mSCovEIw=zeros(size(mre));
mSCovIIw=zeros(size(mre));
vSCovEEw=zeros(size(mre));
vSCovEIw=zeros(size(mre));
vSCovIIw=zeros(size(mre));

mSCorrEEb=zeros(size(mre));
mSCorrEIb=zeros(size(mre));
mSCorrIIb=zeros(size(mre));
vSCorrEEb=zeros(size(mre));
vSCorrEIb=zeros(size(mre));
vSCorrIIb=zeros(size(mre));
mSCovEEb=zeros(size(mre));
mSCovEIb=zeros(size(mre));
mSCovIIb=zeros(size(mre));
vSCovEEb=zeros(size(mre));
vSCovEIb=zeros(size(mre));
vSCovIIb=zeros(size(mre));

for iiii=1:numtrials
    
    [iiii numtrials round(toc(tstart)/60)]
    %save VsNCorrTemp.mat iii mSCorrEE;
    drawnow;
for iii=1:numN

    %[iiii numtrials iii numN round(toc(tstart)/60)]
    
    N=Nrange(iii);
    
    % Number of neurons in each pop
    Ne=round(.8*N);
    Ni=round(.2*N);
    Nx=round(.2*N);
    
    % Number of neurons in each sub-pop
    Ne1=round(Ne/2);
    Ni1=round(Ni/2);
    Nx1=round(Nx/2);
    
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

    % Scaling of connection probability
    % across populations
    qconn=0;

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
    Jx=[sparse(Jxm(1)*binornd(1,Px(1)*2*(1-qconn),Ne1,Nx1)) ...
        sparse(Jxm(1)*binornd(1,Px(1)*2*qconn,Ne1,Nx1)); ...
        sparse(Jxm(1)*binornd(1,Px(1)*2*qconn,Ne1,Nx1)) ...
        sparse(Jxm(1)*binornd(1,Px(1)*2*(1-qconn),Ne1,Nx1)); ...    
        sparse(Jxm(2)*binornd(1,Px(2)*2*(1-qconn),Ni1,Nx1)) ...
        sparse(Jxm(2)*binornd(1,Px(2)*2*qconn,Ni1,Nx1)); ...
        sparse(Jxm(2)*binornd(1,Px(2)*2*qconn,Ni1,Nx1)) ...
        sparse(Jxm(2)*binornd(1,Px(2)*2*(1-qconn),Ni1,Nx1))];
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
    XRCurr = XRCurr(:,2:end-1); % Get rid of first and last 250 ms
    C=cov(XRCurr'); % Compute covariance matrix between all pairs

    % Isolate each block to get XX, XR, RR covs
    [II,JJ]=meshgrid(1:2*Ne,1:2*Ne);II=II(:);JJ=JJ(:);
    CXX=C(II<=Ne & JJ<=Ne & II~=JJ);
    CXR=C(II<=Ne & JJ>Ne & II~=JJ+Ne);
    CRR=C(II>Ne & JJ>Ne & II~=JJ);
    CXXw=C(II<=Ne1 & JJ<=Ne1 & II~=JJ);
    CXXb=C(II<=Ne1 & JJ>Ne1 & JJ<=Ne & II~=JJ);
    CXRw=C(II<=Ne1 & JJ>Ne & JJ<=Ne1+Ne & II~=JJ+Ne);
    CXRb=C(II<=Ne1 & JJ>Ne+Ne1 & II~=JJ);
    CRRw=C(II>Ne & II<=Ne+Ne1 & JJ>Ne & JJ<=Ne1+Ne & II~=JJ);
    CRRb=C(II>Ne & II<=Ne+Ne1 & JJ>Ne+Ne1 & II~=JJ);
    clear C;

    % Comptue means
    mCXX(iiii,iii)=mean(CXX);
    mCXR(iiii,iii)=mean(CXR);
    mCRR(iiii,iii)=mean(CRR);
    mCXXw(iiii,iii)=mean(CXXw);
    mCXXb(iiii,iii)=mean(CXXb);
    mCXRw(iiii,iii)=mean(CXRw);
    mCXRb(iiii,iii)=mean(CXRb);
    mCRRw(iiii,iii)=mean(CRRw);
    mCRRb(iiii,iii)=mean(CRRb);
    
    % and vars
    vCXX(iiii,iii)=var(CXX);
    vCXR(iiii,iii)=var(CXR);
    vCRR(iiii,iii)=var(CRR);
    vCXXw(iiii,iii)=var(CXXw);
    vCXXb(iiii,iii)=var(CXXb);
    vCXRw(iiii,iii)=var(CXRw);
    vCXRb(iiii,iii)=var(CXRb);
    vCRRw(iiii,iii)=var(CRRw);
    vCRRb(iiii,iii)=var(CRRb);
    

    % Now get total input covariances
    TotalCurr=(IeRec+IiRec+IxRec)*winsize;
    TotalCurr = TotalCurr(:,2:end-1); % Get rid of first and last 250 ms
    C=cov(TotalCurr');
    [II,JJ]=meshgrid(1:Ne,1:Ne);II=II(:);JJ=JJ(:);
    CIIw=C(II<=Ne1 & JJ<=Ne1 & II~=JJ);
    CIIb=C(II<=Ne1 & JJ>Ne1);
    CII=C(JJ~=II);

    % Compute means
    mCII(iiii,iii)=mean(CII);
    mCIIw(iiii,iii)=mean(CIIw);
    mCIIb(iiii,iii)=mean(CIIb);

    vCII(iiii,iii)=var(CII);
    vCIIw(iiii,iii)=var(CIIw);
    vCIIb(iiii,iii)=var(CIIb);

    AllRates=hist(s(2,s(1,:)>500),1:N)/(T-500);
    
    mre(iiii,iii)=mean(AllRates(1:Ne));
    mri(iiii,iii)=mean(AllRates(Ne+1:N));
    
    % Compute spike count covariances
    % for all neurons
    C=SpikeCountCov(s,N,250,T-250,winsize);

    % Get spike count corrs over each sub-pop
    [II,JJ]=meshgrid(1:N,1:N);
    SCee=C(II<=Ne & JJ<=Ne & II~=JJ);
    SCeew=C(II<=Ne1 & JJ<=Ne1 & JJ~=II & isfinite(C));
    SCeeb=C(II<=Ne1 & JJ>Ne1 & JJ<=Ne & JJ~=II & isfinite(C));

    SCei=C(II<=Ne & JJ>Ne & JJ~=II & isfinite(C));
    SCeiw=C(II<=Ne1 & JJ>Ne & JJ<=Ne+Ni1 & JJ~=II & isfinite(C));
    SCeib=C(II<=Ne1 & JJ>Ne+Ni1 & JJ~=II & isfinite(C));

    SCii=C(II>Ne & JJ>Ne & JJ~=II & isfinite(C));
    SCiiw=C(II>Ne & II<=Ne+Ni1 & JJ>Ne & JJ<=Ne+Ni1 & JJ~=II & isfinite(C));
    SCiib=C(II>Ne & II<=Ne+Ni1 & JJ>Ne+Ni1 & JJ~=II & isfinite(C));

    % Comptue means
    mSCovEE(iiii,iii)=mean(SCee(:));
    mSCovEI(iiii,iii)=mean(SCei(:));
    mSCovII(iiii,iii)=mean(SCii(:));
    mSCovEEw(iiii,iii)=mean(SCeew(:));
    mSCovEIw(iiii,iii)=mean(SCeiw(:));
    mSCovIIw(iiii,iii)=mean(SCiiw(:));
    mSCovEEb(iiii,iii)=mean(SCeeb(:));
    mSCovEIb(iiii,iii)=mean(SCeib(:));
    mSCovIIb(iiii,iii)=mean(SCiib(:));

    vSCovEE(iiii,iii)=var(SCee(:));
    vSCovEI(iiii,iii)=var(SCei(:));
    vSCovII(iiii,iii)=var(SCii(:));
    vSCovEEw(iiii,iii)=var(SCeew(:));
    vSCovEIw(iiii,iii)=var(SCeiw(:));
    vSCovIIw(iiii,iii)=var(SCiiw(:));
    vSCovEEb(iiii,iii)=var(SCeeb(:));
    vSCovEIb(iiii,iii)=var(SCeib(:));
    vSCovIIb(iiii,iii)=var(SCiib(:));

    
end
end
t0=toc(tstart);

clear stm counts edges edgest edgesi I CSD time C CII CRR CXR CXX IeRec II JJ IiRec IxRec VRec J Jx s sx SCee SCei SCii TotalCurr XRCurr;
clear counts edges edgest edgesi I CSD time C CXX CRR CXR CII CIIw CRRw CXRw CXXw CIIb CRRb CXRb CXXb IeRec II JJ IiRec IxRec VRec J Jx s sx SCee SCei SCii SCeew SCeiw SCiiw SCeeb SCeib SCiib TotalCurr XRCurr;
save Fig6sim.mat
