% Compute cross-spectral density between x and y
% which are discretized by dt
% The processes are broken into blocks of size winrad
% which sets the frequency discretization 
% (dF is proportional to 1/winrad)
% 
% To compute pairwise spectra between a bunch of processes
% let x and y be NtxM and NtxN where M and N are the number 
% of processes in each. The resulting S is NFxMxN
% This only works in Matlab 2017 or newer

function [CSD,F]=GetCrossSpec(x,y,winrad,dt)

    % x and y should be NtxN where Nt is the number
    % of time bins and N is the number of processes   
    if(size(x,1)==1)
        x=x';
    end
    if(size(y,1)==1)
        y=y';
    end


   if(size(x,2)==1 && size(y,2)==1) % If each is just one process
        [CSD,F]=cpsd(x-mean(x),y-mean(y),round(winrad/dt),[],[],1/dt);        
        CSD(2:end)=CSD(2:end)/2;
   else 
       [CSD,F]=cpsd(x-repmat(mean(x),size(x,1),1),y-repmat(mean(y),size(y,1),1),round(winrad/dt),[],[],1/dt,'mimo');
       CSD(2:end,:)=CSD(2:end,:)/2;
   end
   
end

