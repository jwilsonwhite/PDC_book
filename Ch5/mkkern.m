function kxy = mkkern(x,y,F,fixparm,T)

% This function creates the IPM kernel.
% Adapted from Easterling et al. (2001) for use with the red urchin model
% Note that this one uses a gamma distribution for Linf

%Define which sizes can be fished:
isjuv = 1 - normcdf(x,fixparm(5),fixparm(7)); 

%define pr(reproductive) using a maturity ogive
%not needed now because open population
%ismat = 1 - normcdf(x,fixparm(5),diff(x(1,1:2))/2);

%SURVIVAL PART OF KERNEL
M = fixparm(4); % natural mortality rate
if isnan(M) % declining mortality rate for urchin model
m = 5.81*exp(-0.043*x);    
    
else
m = ones(size(x)).*M + (1-isjuv).*F; %this is a matrix size x, mortality for each size 
end

p1 = exp(-m*T); % convert mortality rate to survivorship, iterate over time steps


%GROWTH PART OF KERNEL
Linf = fixparm(1);
k = fixparm(2);

%add variability around von Bertalanffy growth with a distribution of Linfs
nn = 1e3;
Linfs = normrnd(Linf,fixparm(7),[1,1,nn]);

pmean1 = Linfs - (Linfs - repmat(x,[1,1,nn])).*exp(-k*T);

%evaluate growth part of kernel
p2 = mean( normpdf(repmat(y,[1 1 nn]),pmean1,2*diff(x(1,1:2))),3);

%ensure integrates to 1:
p2 = p2./(repmat(sum(p2),[size(y,1),1])*diff(y(1:2))); % midpoint rule

p1 = max(0,p1);    %to make sure no negatives
p2 = max(0,p2);    %to make sure no negatives

kxy = p1.*p2;

%add fecundity to kernel if a closed population
if length(fixparm)>7
   Rvec = normpdf(y,fixparm(8),fixparm(9));
   Fec = fixparm(10)*x.^fixparm(11); 
   Fec = Fec/max(Fec(:)); % rescale for visualization
   Q = Fec.*Rvec;
   kxy = kxy + Q;
end

