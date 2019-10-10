function filter_example

% Make figure showing example of population filtering colored noise:

% Inspired in part by code shared by Hristo Zhivomirov
% (https://www.mathworks.com/matlabcentral/fileexchange/42919-pink--red--blue-and-violet-noise-generation-with-matlab-implementation/content/bluenoise.m)
rng(8) %4
N = 2^7;% number of samples (must be power of 2)
f = linspace(0,0.5,N/2); % frequencies
NumUniquePts = N/2+1;

%-----------------------------------------
% White noise:
Wh = normrnd(0,1,[N,1]);

% FFT
y = fft(detrend(Wh),N)/N;
Wh_f1 = y; % the raw FFT
y = y(2:(N/2+1));
Wh_f = abs(y); % magnitude
%------------------------------------------

%------------------------------------------
% Red noise:
x = normrnd(0,1,[N,1]);
X = fft(x); 
n = 1:NumUniquePts;

% multiply to make power spectrum proportional to 1/(f^2)
X(1:NumUniquePts) = X(1:NumUniquePts)./(n(:));
% Do the same for the right half:
X((NumUniquePts+1):N) = real(X(N/2:-1:2)) -1i*imag(X(N/2:-1:2));

Red = real(ifft(X));
Red = Red(1:N);
Red = (Red-mean(Red))./std(Red);

% FFT
y = fft(Red,N)/N;
Red_f1 = y; % raw FFT
y = y(2:N/2+1);
Red_f = abs(y); % magnitude
%------------------------------------------

%------------------------------------------
% A filter
Filt = normpdf(f,0.1,0.05);
Filt = Filt./sum(Filt);
%------------------------------------------

%------------------------------------------
% Apply filter
Wh_fn = Wh_f1(:).*[Filt(:);flipud(Filt(:))];
Red_fn = Red_f1(:).*[Filt(:);flipud(Filt(:))];

% Recover the new signal:
Wh_n = real(ifft(Wh_fn,N))*N;
Wh_n = Wh_n./std(Wh_n);
Red_n = real(ifft(Red_fn,N))*N;
Red_n = Red_n./std(Red_n);


% Check that the FFT matches
Wh_nf = abs(fft(Wh_n,N))/N;
Red_nf = abs(fft(Red_n,N))/N;

Wh_nf = Wh_nf(1:N/2);
Red_nf = Red_nf(1:N/2);

figure(1)
clf
set(gcf,'units','cent','position',[10 10 21 10])
Gray = [0.5 0.5 0.5];

subplot(2,3,1)
set(gcf,'units','cent','position',[10 10 21 10])
hold on
plot(2+Wh,'color',Gray)
plot(-2+Red,'-','color','k')
set(gca,'ytick',[],'xtick',0:20:160)
set(gca,'xlim',[0 124],'xtick',0:40:160)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
xlabel('Years','fontsize',12)
ylabel('Population size','fontsize',12)

subplot(2,3,4)
hold on
plot(-log2(f),2*Wh_f.^2,'color',Gray)
plot(-log2(f),2*Red_f.^2,'color','k')
set(gca,'xlim',[1 6],'xtick',1:6,'xticklabel',1./(2.^(-1:-1:-6)))
set(gca,'ylim',[0 0.1],'ytick',0:0.05:0.4)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
xlabel('Period','fontsize',14)
ylabel('Variance','fontsize',12)

subplot(2,3,5)
plot(-log2(f),Filt,'k')
set(gca,'ytick',[])
set(gca,'xlim',[1 6],'xtick',1:6,'xticklabel',1./(2.^(-1:-1:-6)))
set(gca,'tickdir','out','ticklength',[0.015 0.015])
xlabel('Period','fontsize',14)
ylabel('Variance','fontsize',12)

subplot(2,3,3)
hold on
plot(2+Wh_n,'color',Gray)
plot(-2+Red_n,'color','k')
set(gca,'ytick',[],'xtick',0:20:160)
set(gca,'xlim',[0 124],'xtick',0:40:160)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
xlabel('Years','fontsize',12)
ylabel('Population size','fontsize',12)

subplot(2,3,6)
hold on
plot(-log2(f),2*(Wh_nf).^2,'color',Gray)
plot(-log2(f),2*(Red_nf).^2,'color','k')
set(gca,'xlim',[1 6],'xtick',1:6,'xticklabel',1./(2.^(-1:-1:-6)))
set(gca,'ylim',[0 0.25],'ytick',0:0.1:0.4)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
xlabel('Period','fontsize',14)
ylabel('Variance','fontsize',12)


% Define a transfer function:


