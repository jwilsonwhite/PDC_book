function Leslie_transients

% Project Leslie matrices to illustrate transient dynamics. Same matrix as
% used in Fig 3.5

% Assemble the matrix:
F = [0 0 0 0 0.8 1 1 1 0.1]; % fecundity
S = [0.7 0.7 0.9 0.9 0.9 0.9 0.3 0.3]; % survival
Smat = [diag(S), zeros(length(S),1)]; % lower half of matrix
L = [F;Smat];  % the full matrix

% display the eigenvalues & magnitudes:
eig(L) % eigenvalues
abs(eig(L)) % magnitudes
[W,Lambda]=eig(L); % V is right eigenvectors, Lambda is eigenvalues
Lambda = diag(Lambda); % extract diagonal vector from the matrix

% Calculate generation time:
A = 1:length(F);
G = F.*[1,cumprod(S)];
Gt = G*A'/sum(G) % generation time

% Period of oscillation:
Pd = acos(real(Lambda(2))/abs(Lambda(2)))
2.*pi./Pd

% Initial conditions:
N0 = [1000, 600, 400, 300, 200, 0 0 0 0]';

% Project:
T = 24;
N = nan(length(N0),T);
N(:,1) = N0;
for t = 2:T
    N(:,t) = L*N(:,t-1);
end

% Express transient dynamics in terms of eigenvalues:
C = W\N0; % coefficients c in the eigenvector equation


% Now for a population with much constrained reproduction:
S2 = S;
S2(4:end) = S2(4:end)*0.25;
F2 = F*9.5;
Smat2 = [diag(S2), zeros(length(S2),1)]; % lower half of matrix
L2 = [F2;Smat2];  % the full matrix
% New generation time:
G = F.*[1,cumprod(S2)];
Gt = G*A'/sum(G)

[W2,Lambda2]=eig(L2); % V is right eigenvectors, Lambda is eigenvalues
Lambda2 = diag(Lambda2); % extract diagonal vector from the matrix

max(Lambda2)
[atan(imag(Lambda)./real(Lambda)),atan(imag(Lambda2)./real(Lambda2))]

% Project:
T = 24;
N2 = nan(length(N0),T);
N2(:,1) = N0;
for t = 2:T
    N2(:,t) = L2*N2(:,t-1);
end


    % Calculate distance metrics
    T2 = 100;
    NN = nan(length(N0),T2);
    NN(:,1) = N0;
    for t = 2:T2
    NN(:,t) = L*NN(:,t-1);
    end
    WW = W(:,1);
    WW = abs(WW)./sum(abs(WW));
    CC = abs(C(1));
    
    Delta = abs(sum(NN./(max(Lambda).^(1:T2)) - repmat(WW,[1,T2]))); % Duncan & Duncan Delta
    
  %  s = sum(NN./(max(Lambda).^(1:T2)) - repmat(CC.*NN(:,1),[1,T2]));
    s = sum(NN./(max(Lambda).^(1:T2)) - repmat(abs(CC.*W(:,1)),[1,T2]));

    
    D1 = flipud(abs(cumsum(flipud(s(:)))));
    D2 = flipud(cumsum(flipud(abs(s(:)))));
    

% Plots:



% Figure 1: timeseries of total population abundance
figure(1)
clf
set(gcf,'units','cent','position',[10 10 21 15])

% Plot eigenvalues in complex plane:
subplot(2,3,1) 
hold on
do_unit_circle
for l = 1:length(Lambda)
    plot(real(Lambda(l)),imag(Lambda(l)),'ko','markerfacecolor','k')
    text(real(Lambda(l))+0.075,imag(Lambda(l)),strcat('\lambda','_',num2str(l)),'fontsize',10)
end
format_complex_plane


% Plot transient dynamics:
subplot(2,3,2:3)
hold on
plot(1:T,sum(N),'k-')
set(gca,'xtick',1:5:100,'xticklabel',0:5:100,'ytick',0:400:10000)
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xminorgrid','on','xminortick','on','xgrid','on')
set(gca,'ylim',[2000 3300])
set(gca,'xlim',[0 24])
    

% Plot eigenvalues in complex plane:
subplot(2,3,4) 
hold on
do_unit_circle
for l = 1:length(Lambda2)
    plot(real(Lambda2(l)),imag(Lambda2(l)),'ko','markerfacecolor','k')
    text(real(Lambda2(l))+0.075,imag(Lambda2(l)),strcat('\lambda','_',num2str(l)),'fontsize',10)
end
format_complex_plane
abs(Lambda2)


% Plot transient dynamics
subplot(2,3,5:6)
hold on
plot(1:T,sum(N2),'k-')

set(gca,'xtick',1:5:100,'xticklabel',0:5:100,'ytick',0:200:10000)
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xminorgrid','on','xminortick','on','xgrid','on')
set(gca,'ylim',[2000 3300])
set(gca,'xlim',[0 24])

%----------------------------------------------
% Figure 2: 
% Plot the iteration of 2nd & 3rd eigenvalues
figure(2)
clf
set(gcf,'units','cent','position',[10 10 21 10])

% Plot eigenvalues in complex plane:
subplot(2,4,1) 
hold on
do_unit_circle
Lam2s = Lambda(2).^(1:T);
plot(Lam2s,'ko')

for t = 1:(T-1)
    plot([real(Lam2s(t)),real(Lam2s(t+1))],[imag(Lam2s(t)),imag(Lam2s(t+1))],'k-')
end
format_complex_plane

subplot(2,4,2) 
hold on
do_unit_circle
Lam3s = Lambda(3).^(1:T);
plot(Lam3s,'kd')

for t = 1:(T-1)
    plot([real(Lam3s(t)),real(Lam3s(t+1))],[imag(Lam3s(t)),imag(Lam3s(t+1))],'k-')
end
format_complex_plane


% Plot sequence of eigenvector terms
subplot(2,4,3:4)
hold on

E = Lam2s.*C(2).*sum(W(:,2)) + Lam3s.*C(3).*sum(W(:,3)); % 2nd & 3rd terms in eigen expansion


plot(1:(T+1),[0, real(E)],'k-')

set(gca,'xtick',1:5:100,'xticklabel',0:5:100,'ytick',-1000:200:1000)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xminorgrid','on','xminortick','on','xgrid','on')
set(gca,'ylim',[-400 400])
set(gca,'xlim',[0 24],'xcolor','k','ycolor','k')

% plot 2&3 eigenvalues back on Fig 1b
figure(1)
subplot(2,3,2:3)
plot(1:(T+1),sum(N(:,1))+[0,real(E)],'k--')

% Second row: 4th & 5th:
figure(2)
subplot(2,4,5) 
hold on
do_unit_circle
Lam4s = Lambda(4).^(1:T);
plot(Lam4s,'ko')

for t = 1:(T-1)
    plot([real(Lam4s(t)),real(Lam4s(t+1))],[imag(Lam4s(t)),imag(Lam4s(t+1))],'k-')
end
format_complex_plane

subplot(2,4,6) 
hold on
do_unit_circle
Lam5s = Lambda(5).^(1:T);
plot(Lam5s,'kd')

for t = 1:(T-1)
    plot([real(Lam5s(t)),real(Lam5s(t+1))],[imag(Lam5s(t)),imag(Lam5s(t+1))],'k-')
end
format_complex_plane


% Plot sequence of eigenvector terms
subplot(2,4,7:8)
hold on

E = Lam4s.*C(4).*sum(W(:,4)) + Lam5s.*C(5).*sum(W(:,5)); % 2nd & 3rd terms in eigen expansion


plot(1:(T+1),[0, real(E)],'k-')

set(gca,'xtick',1:5:100,'xticklabel',0:5:100,'ytick',-1000:50:1000)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xminorgrid','on','xminortick','on','xgrid','on')
%set(gca,'ylim',[-400 400])
set(gca,'xlim',[0 24],'xcolor','k','ycolor','k')

% plot 4&5 eigenvalues back on Fig 1b
figure(1)
subplot(2,3,2:3)
plot(1:(T+1),sum(N(:,1))+[0,real(E)],'k-.')

% Figure 3: transient convergence metrics:

figure(3)
clf
set(gcf,'units','cent','position',[20 10 9 18])
subplot(2,1,1)
plot(1:T,sum(N),'k')
set(gca,'xlim',[0 24],'xcolor','k','ycolor','k')
set(gca,'tickdir','out','ticklength',[0.02 0.02])

subplot(2,1,2)
hold on
plot(1:T,Delta(1:T),'k-')
Pos = get(gca,'position');
set(gca,'color','none')
set(gca,'xlim',[0 24],'xcolor','k','ycolor','k')
set(gca,'ytick',[])
axes
set(gca,'position',Pos)
plot(1:T,D1(1:T),'k--')
plot(1:T,D2(1:T),'k--')

set(gca,'xlim',[0 24],'xcolor','k','ycolor','k')
set(gca,'yaxislocation','right','color','none')
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'ytick',[])
set(gca,'position',Pos)


function do_unit_circle

% plot unit circle
Theta = linspace(0,2*pi,100);
[x,y] = pol2cart(Theta,1);
plot(x,y,'k-')
plot([-2 2],[0 0],'k--')
plot([0 0],[-2 2],'k--')

function format_complex_plane
axis equal
xlim([-1.2 1.2])
ylim([-1.2 1.2])
set(gca,'tickdir','out','ticklength',[0.02 0.02])
ylabel('Im','fontsize',12)
xlabel('Re','fontsize',12)



