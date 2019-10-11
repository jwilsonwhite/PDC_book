function example_dndm_figs

% Some plots showing the components of dndm

figure(1)
clf
set(gcf,'units','cent','position',[10 10 18 10])


% Parameters for red sea urchin:
k = 0.3;
Linf = 100;
sig_k = 0.01;
sig_Linf = 100;

A = 0:50;
L = 0:0.1:130;
Rmu = 60;

% plot the different growth curves:
Nc = 1e4;
gamk1 = k^2/sig_k; % shape param
gamk2 = k/gamk1; % scale param 
ks = gamrnd(gamk1,gamk2,[Nc,1]);
gamL1 = Linf^2/sig_Linf; % shape param 
gamL2 = Linf/gamL1; % scale param
Linfs = gamrnd(gamL1,gamL2,[Nc,1]);

% Panel 1 ...........................................
% Size distribution: n vs m
subplot(2,2,1)
hold on

mu =0;
sd = NaN;
for n = 1:Nc
p = [k, Linfs(n), mu, sd, Rmu];
N(:,n) = SizeStruct_VarR_EqVB(1,p,A,L,'case2');
end

N = mean(N,2);
plot(L,N,'k-')

% Plot the tangent line and growth rate (g) at m = 95 cm
Lx = 95;
Lind = find(L == Lx);
Lind = [Lind Lind+1];
dndm = diff(N(Lind))./diff(Lind);



plot([Lx Lx+diff(L(Lind))*30],[N(Lind(1)),N(Lind(1))+diff(N(Lind))*30],'r-')

% Growth rate:
Gs = mean(k*(Linfs - Lx));

plot([Lx Lx+20],[N(Lind(1)) N(Lind(1))*(1 + Gs)],'g-')

set(gca,'xcolor','k','ycolor','k','xtick',[],'ytick',[],'ylim',[0 600],'xlim',[0 130])
xlabel('Size (m)','fontsize',12)
ylabel('Abundance [n(m,t)]','fontsize',12)
%............................................................

%............................................................
% Subplot 2: dg/dm example
subplot(2,2,2)

G = L.^-0.05;
plot(L,G,'k-')

set(gca,'xcolor','k','ycolor','k','xtick',[],'ytick',[],'ylim',[0.7 1.2],'xlim',[0 130])
xlabel('Size (m)','fontsize',12)
ylabel('Growth rate (g)','fontsize',12)
%............................................................

%............................................................
% Subplot 3: alternative dg/dm example
subplot(2,2,3)

G = L.^0.3;
plot(L,G,'k-')

set(gca,'xcolor','k','ycolor','k','xtick',[],'ytick',[],'ylim',[0 5],'xlim',[0 130])
xlabel('Size (m)','fontsize',12)
ylabel('Growth rate (g)','fontsize',12)
%............................................................

%............................................................
% Subplot 4: characteristic lines.
% Assume von Bert for simpliciy
subplot(2,2,4)
hold on

Lags = [0 30 60];
k = 0.05;
for l = 1:3

Ls = Linf*(1 - exp(-k*A));

plot(A+Lags(l),Ls,'k-')
end

plot([30 30],[-10 110],'k--')


set(gca,'xcolor','k','ycolor','k','xtick',[],'ytick',[],'ylim',[-10 110],'xlim',[-10 120])
xlabel('Time (t)','fontsize',12)
ylabel('Size (m)','fontsize',12)
%............................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = SizeStruct_VarR_EqVB(t,p,A,L,mortcase)
% size-based von Foerster equation
% Pulsed recruitment, no variability in von Bertalanffy equation

k = p(1);
Linf = p(2);
mu = p(3);
sd = p(4);
Rmu = p(5);

% von Bertalanffy growth
L0 = 0;
A=log(1-L/Linf)./-k; % inverse of vB equation
A(L>Linf) = 0;
A(A==0) = max(A); % substitute the maximum value


gl0 = k*(Linf-L0);
gl =  k*(Linf-L);

if isnan(sd) % is recruitment pulsed or constant?
% Constant recruitment
R = Rmu;
else
% Pulsed recruitment. Mean = Jun 30. Time is fractions of a year
tt = mod(t-A,1); % take fraction of the year for all ages going into the past
R = normpdf(tt,mu,sd);
dA = diff(A);
dA = [dA(1) dA]; % pad with a dummy first entry to extend length
R = R./mean(R);
R = R.*Rmu;
end

% Mortality rate: (several different options)
switch mortcase
    case 'case1'
% unimodal mortality rate:
D = 0.60 - 0.59/50^2*(L-50).^2; % min at 0.01, peak of 0.6 at 50 mm 
D = max(D,0.01);
    case 'case2'
% constant mortality rate:
D = 0.1;
    case 'case3'
% low mortality for small individuals:
D = ones(size(L)).*0.1;
D(L<40) = 0;
    case 'case4'
% high mortality for small individuals:
D = ones(size(L)).*0.1;
D(L<80) = 1;
    case 'case5'
% high mortality for medium individuals:
D = ones(size(L)).*1.0;
D(L<40) = 0;
D(L>80) = 0.1;
end

% Integrate D/g(l) wrt l
Int = cumsum(D./gl).*diff(L(1:2));

% Solution:
N = R.*gl0./gl.*exp(-Int);
N(L>0.99*Linf)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%