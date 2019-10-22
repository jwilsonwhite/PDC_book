function invasion

% Simulate basic invasion models. 

X = 400; % size of spatial domain
D = 4; % diffusion coefficient
r = 0.3; % growth rate
T = 40; % time scale

figure(1)
clf
set(gcf,'units','cent','position',[2 30 18 21])

%---------------------------------
% 1) Reaction-diffusion with exponential growth
K = Inf; % no carrying capacity
N0 = 0.5; % initial seed population
N = zeros(1,X); % 
N(round(X/2)) = N0; % initial state
dt = 0.1;
v = 0; % velociy (left to right)
Intw = 400; % interval width for plotting

% Run the model
Nt = reactiondiffusion(N,r,K,D,X,dt,T,v);

% Plotting
subplot(3,2,1) % spatial spread
Intervals = 1:dt*Intw:size(Nt,1);
plot(Nt(Intervals,:)','k-')
prettyplot_space(gca,false,false)

subplot(3,2,2) % rate of spread (to the right):
hold on
ASR = 2*sqrt(r*D); % asymptotic rate of spread
thresh = 1e-3; % threshold density for enumerating invasion
plot_ASR(Nt,X,ASR,Intervals,thresh,dt)
prettyplot_ASR(gca,T,false,false)


% End reaction-diffusion with exponential growth
%---------------------------------

%---------------------------------
% 1) Reaction-diffusion with logistic growth
K = 0.5; % no carrying capacity
N0 = 0.5; % initial seed population
N = zeros(1,X); % 
N(round(X/2)) = N0; % initial state
dt = 0.1;
v = 0; % velociy (left to right)
Intw = 400; % interval width for plotting

% Run the model
Nt = reactiondiffusion(N,r,K,D,X,dt,T,v);

% Plotting
subplot(3,2,3) % spatial spread
Intervals = 1:dt*Intw:size(Nt,1);
plot(Nt(Intervals,:)','k-')
prettyplot_space(gca,true,false)

subplot(3,2,4) % rate of spread (to the right):
hold on
ASR = 2*sqrt(r*D); % asymptotic rate of spread
thresh = 1e-3; % threshold density for enumerating invasion
plot_ASR(Nt,X,ASR,Intervals,thresh,dt)
prettyplot_ASR(gca,T,true,false)
% End reaction-diffusion with logistic growth
%---------------------------------


%---------------------------------
% Leptokurtic dispersal using integrodifference equation
type = 'lepto';
A = 3; % diffusion coefficient
%r = 0.5; % growth rate
K = 0.5; % carrying capacity
N0 = 0.5; % initial seed population
N = zeros(1,X); % 
N(round(X/2)) = N0; % initial state
v = 0; % velociy (left to right)

% Run the model
switch type
    case 'gauss'
Nt = integrodiff(N,r,K,sqrt(2*D),X,T,v,type); % use the same D as in the PDE version above
    case 'lepto'
Nt = integrodiff(N,r,K,A,X,T,v,type);
end

% Plotting
subplot(3,2,5) % spatial spread
Intervals = 1:2:T;
plot(Nt(:,Intervals),'k-')
prettyplot_space(gca,false,true)

subplot(3,2,6) % rate of spread (to the right):
hold on
switch type
    case 'gaussian'
ASR = 2*sqrt(r*D); % asymptotic rate of spread
    case 'lepto'
        % Find the SD of the distribution to get the equivalent D
        xx = -100:0.1:100;
        yy = A^2*exp(-A*sqrt(abs(xx)))/4;
        sd = sqrt(sum(yy.*abs(xx).^2.)/sum(yy));
        ASR = sqrt(2)*sd*sqrt(r);
end % end switch type
thresh = 1e-3; % threshold density for enumerating invasion
plot_ASR(Nt',X,ASR,Intervals,thresh,1)
prettyplot_ASR(gca,T,false,true)
% End leptokurtic
%---------------------------------

%---------------------------------
% plot the kernels:
figure(2)
clf
set(gcf,'units','cent','position',[12 13 4 2.5])
hold on
X = -10:0.01:10;
%D = 4;
%A = 2;
Gauss = normpdf(X,0,sqrt(2*D));
%Lepto = A^2*exp(-A*sqrt(abs(X)))/4;
Lepto = (1/A)*exp(-abs(X)/A)/2;
plot(X,Gauss,'k-')
plot(X,Lepto,'k:')
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',-10:10:10,'ytick',[])
%close gcf

%---------------------------------
% reaction-diffusion model (with logistic growth)
function Nt = reactiondiffusion(N,r,K,D,X,dt,T,v)
Nt = [];
for t = 0:dt:T
    Nt = [Nt; N];    
    % calculate spatial gradient in n
    nL = [N(X), N(1:(X-1))];
    nR = [N(2:X), N(1)];
    dNdt = r*N.*(1-N/K) + D.*(nL + nR - 2*N)+ v.*(nL-nR); % (the second-to-last term is the 2nd derivative of n wrt x)
    N = N + dNdt*dt;
end % end T
%---------------------------------

%---------------------------------
function plot_ASR(Nt,X,ASR,Intervals,thresh,dt)
Nd = Nt(:,round(X/2):end);
Nd = Nd>thresh;
Nd = Nd.*repmat(1:(round(X/2)+1),[size(Nd,1),1]);
Nd = max(Nd,[],2);
Yint = Nd(Intervals(2)); % approximately when the asymptotic rate of spread (ASR) is reached

plot((Intervals-1)*dt,Nd(Intervals),'ko')
plot((Intervals-1)*dt,ASR.*(Intervals-1)*dt+Yint-ASR,'k-')
%---------------------------------

%---------------------------------
function N = integrodiff(N0,r,K,D,X,T,v,type)

N = zeros(X,T);
N(:,1) = N0; % initial conditions

% calculate dispersal matrix:
x = repmat(1:X,[X,1]);
y = x';
d = x-y;
switch type
    case 'gauss'
        KK = normcdf(d+0.5,v,D) - normcdf(d-0.5,v,D); % normal dispersal kernel
    case 'lepto'
        Sy = makeSimpVec(1,X); % Simpson's integration (midpoint rule is too coarse)
      %  KK = D^2*exp(-D*sqrt(abs(d)))/4;
        KK = (1/D)*exp(-abs(d)/D)./2;
        KK = KK.*repmat(Sy(:)',[X,1]); % integrate each column
        
end %end switch

for t = 2:T
    fN = N(:,t-1).*exp(r*(1-N(:,t-1)/K)); % discrete-time logistic
    N(:,t) = (KK)*fN;    % division by 2 integrates using midpoint rule 
end % end T
%---------------------------------

%---------------------------------
function prettyplot_space(h,YL,XL)

set(h,'tickdir','out','ticklength',[0.02 0.02])
set(h,'xlim',[100 300],'xtick',100:50:300,'xticklabel',-100:50:100)
set(h,'ylim',[0 1],'ytick',0:0.2:1)
set(h,'fontsize',10)
if YL
    ylabel(h,'Population density (N)','fontsize',14)
end
if XL
    xlabel(h,'Distance (km)','fontsize',14)
end

function prettyplot_ASR(h,T,YL,XL)

set(h,'tickdir','out','ticklength',[0.02 0.02])
set(h,'xlim',[0 T],'xtick',0:(T/4):T)
set(h,'fontsize',10)
%set(h,'ylim',[0 1],'ytick',0:0.2:1)
if YL
    ylabel(h,'Distance spread (km)','fontsize',14)
end
if XL
    xlabel(h,'Time (y)','fontsize',14)
end
    




