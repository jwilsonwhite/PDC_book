function central_limit_theorem

% Demonstrate the properties of the central limit theorem

% parameters of a gamma distribution
A = 2;
B = 3;

% sample size
N = 25;

% number of simulated samples
T = 1e7;

% abscissa coordinates
X = linspace(0,1e2,1e4);


% Setup figure
figure(3)
set(gcf,'units','cent','position',[30,30,18,15])
clf


% plot the gamma distribution
sp(1) = subplot(3,2,1:2);
plot(X,gampdf(X,A,B),'k','linewidth',2)


% simulate sampling the distribution
XX = gamrnd(A,B,[N,T]);
Edges = linspace(0,15,1e3);

% Plot the mean for n =5 (only take first 5 rows of XX)
sp(2) = subplot(3,2,3);
histogram(mean(XX(1:5,:)),Edges,'normalization','pdf','facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5])

% Plot the mean for n = 25 (take all rows of XX)
sp(3) = subplot(3,2,5);
histogram(mean(XX(1:25,:)),Edges,'normalization','pdf','facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5])

% Sum for n = 5
sp(4) = subplot(3,2,4);
histogram(sum(XX(1:5,:)),Edges*5,'normalization','pdf','facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5])

% Sum for n = 25
sp(5) = subplot(3,2,6);
histogram(sum(XX(1:25,:)),Edges*25,'normalization','pdf','facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5])


% Format plots:
set(sp(1),'ylim',[0 0.15],'ytick',[],'tickdir','out','ticklength',[0.015 0.015],'xlim',[0 40])
set(sp(2),'ylim',[0 0.4],'ytick',[],'tickdir','out','ticklength',[0.02 0.02])
set(sp(3),'ylim',[0 0.6],'ytick',[],'tickdir','out','ticklength',[0.02 0.02])
set(sp(4),'ylim',[0 0.08],'ytick',[],'tickdir','out','ticklength',[0.02 0.02])
set(sp(5),'ylim',[0 0.025],'ytick',[],'tickdir','out','ticklength',[0.02 0.02])
set(sp(2),'xtick',[0:3:15],'xlim',[0 15])
set(sp(3),'xtick',[0:3:15],'xlim',[0 15])
set(sp(4),'xtick',[0:3:15]*5,'xlim',[0 75])
set(sp(5),'xtick',[0:3:15]*25,'xlim',[0 375])
ylabel(sp(2),'Probability density','fontsize',12)
xlabel(sp(1),'X','fontsize',12)
xlabel(sp(3),'mean','fontsize',12)
xlabel(sp(5),'sum','fontsize',12)

% get actual values to print on the screen:
axes(sp(1))
text(30,0.1,strcat('m = ',num2str(mean(XX(:)))))
text(30,0.08,strcat('s2 = ',num2str(var(XX(:)))))

axes(sp(2))
text(10,0.3,strcat('m = ',num2str(mean(mean(XX(1:5,:))))))
text(10,0.26,strcat('s2 = ',num2str(var(mean(XX(1:5,:))))))

axes(sp(3))
text(10,0.4,strcat('m = ',num2str(mean(mean(XX(1:25,:))))))
text(10,0.32,strcat('s2 = ',num2str(var(mean(XX(1:25,:))))))

axes(sp(4))
text(60,0.06,strcat('m = ',num2str(mean(sum(XX(1:5,:))))))
text(60,0.05,strcat('s2 = ',num2str(var(sum(XX(1:5,:))))))

axes(sp(5))
text(300,0.02,strcat('m = ',num2str(mean(sum(XX(1:25,:))))))
text(300,0.015,strcat('s2 = ',num2str(var(sum(XX(1:25,:))))))

