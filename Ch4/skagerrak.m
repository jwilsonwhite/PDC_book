function skagerrak

% Make plots of Skagerrak cod dynamics from Bjornstad et al. 2001

D = importdata('skagerrak_clean.csv');
Y = D.data(:,1);
X = D.data(:,2);
X = log10(X);

% only consider post-WWII for FFT
X2 = X(Y>=1946);
N = length(X2);

% FFT
nfft= 2^nextpow2(N);
y = fft(detrend(X2),nfft)/N;
f = 1/2*linspace(0,1,nfft/2+1);
y = y(1:nfft/2+1);
y2 = (abs(y));
y2 = 2*y2.^2; % variance

y3 = fit(f(:),y2(:),'poly4'); % polynomial fit to the FFT spectrum


% Plotting
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 21])

% Subplot 1: log10 abundance data
subplot(3,1,1)
plot(Y,X,'marker','o','color','k','markersize',4)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'fontsize',10,'xcolor','k','ycolor','k')
set(gca,'ytick',1:4,'yticklabel',[10 100 1000 10000])
xlabel(gca,'Year','fontsize',14)
ylabel(gca,'Age-0 abundance','fontsize',14)


% Subplot 2: Spectrum of population variability:
subplot(3,1,2)
hold on
plot(f,y2,'k-')
h=plot(y3);
set(h,'color',[0.5 0.5 0.5],'linewidth',2);

legend hide
%set(gca,'xlim',[1 6],'xtick',1:6,'xticklabel',1./(2.^(-1:-1:-6)))
set(gca,'xlim',[0 0.5],'xtick',0:0.1:1,'xcolor','k','ycolor','k')
set(gca,'ytick',0:0.01:1)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
xlabel('Frequency','fontsize',14)
ylabel('Variance','fontsize',14)


% Transfer functions, from Bjornstad et al. 2004:

Gam = 0.3;
Bet = 0.4;
Lam = 0.4;

% Result for the limit as d -> Inf
%Tr = (1 + Gam*exp(-2*pi*1i*f))./...
%     (1 + Gam*1i*exp(-2*pi*f)-exp(-2*pi*1i*f).*(1-Bet).*(1-Lam)./...
%     (exp(2*pi*1i*f)-Lam) );

% Result for the case with d = 6 (realistic lifespan for Skagerrak)
f = linspace(0,0.5,1e3);
d = 6; m = 2;
k = 1:d;
Pk = Lam.^(k-1)./cumsum(Lam.^(k-1));
Pk = Pk(m:end);
Pk = Pk./sum(Pk);

k = k(m:end);

% To calculate summation in denominator:
fmat = repmat(f(:)',[length(k),1]);
kmat = repmat(k(:),[1,length(f)]);
Pmat = repmat(Pk(:),[1,length(f)]);
Denom = sum(exp(-2*pi*1i*kmat.*fmat).*Pmat);

Tr = ((1-Bet) .* Denom)./...
     (1 + Gam*exp(-2*pi*1i*f)-...
     (1-Bet).*Denom);


% Plot transfer function
 subplot(3,1,3)
 plot(f,abs(Tr),'k-','linewidth',2)
set(gca,'xlim',[0 0.5],'xtick',0:0.1:1,'xcolor','k','ycolor','k')
set(gca,'tickdir','out','ticklength',[0.015 0.015])
xlabel('Frequency','fontsize',14)
ylabel('Variance','fontsize',14)
 

