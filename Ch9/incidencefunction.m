function incidencefunction

% Simulate incidence function model used by Pvaskainen & Hanski (2003, TPB)

% 1) Define patch locations and sizes:
Xc = [-0.5, 0.5,0, 0, 0, 0];
Yc = [0, 0, -2, 2, 0.25, -0.25];
S = [1, 0.5, 2, 0.5, 2, 2];
n = length(S);

% 2) Distance matrix
D = dist([Xc(:),Yc(:)]');

% 3) Metapopulation model parameters:
c = 0.5;
e = 0.5;
alfa = 0.5;

% 4) Connectivity matrix:
Ai = repmat(S(:),[1,n]); % destination
Aj = repmat(S(:)',[n,1]); % origin
C = exp(-alfa*D).*Ai.*Aj .* (1-eye(n));

%  Patch values
[V,L]=eig(C);
L = diag(L);
V1 = V(:,L==max(L));
[W,~]=eig(C');
W1 = W(:,L==max(L));

Val = V1.*W1.*max(L);
Val = Val./sum(Val);
max(L)

% Perturbation patch measures:
for i = 1:n
    OK = true(n,1);
    OK(i) = false;
    Ct = C(OK,OK);
    L_p(i) = max(eig(Ct));
end
L_p = max(L) - L_p;
L_p = L_p./sum(L_p);

% Run the model, deleting patches:
% Baseline:
Nsim = 3;
T = 100;
N = zeros(n,T,Nsim);
N(:,1,:) = 1; % round(rand([n,1,Nsim]));
for i = 1:Nsim
for t = 2:T
    Col = c*exp(-alfa*D).*Aj .* (1-eye(n));
    Ext = e/S(:);
    Col = min(max(Col,0),1); 
    Ext = min(max(Ext,0),1);
    P = Col*N(:,t-1,i) - Ext(:).*N(:,t-1,i);
    P = min(max(P,0),1);
    N(:,t,i) = binornd(1,P); 
end
end
Nbase = N;


% Remove patch 6:
for i = 1:Nsim
for t = 2:T
    Col = c*exp(-alfa*D).*Aj .* (1-eye(n));
    Ext = e/S(:);
    Col = min(max(Col,0),1); 
    Ext = min(max(Ext,0),1);
    P = Col*N(:,t-1,i) - Ext(:).*N(:,t-1,i);
    P = min(max(P,0),1);
    P(6) = 0;
    N(:,t,i) = binornd(1,P); 
end
end
N6 = N(1:5,:,:);
C6 = C(1:5,1:5);
max(eig(C6))

% Remove patch 1:
for i = 1:Nsim
for t = 2:T
    Col = c*exp(-alfa*D).*Aj .* (1-eye(n));
    Ext = e/S(:);
    Col = min(max(Col,0),1); 
    Ext = min(max(Ext,0),1);
    P = Col*N(:,t-1,i) - Ext(:).*N(:,t-1,i);
    P = min(max(P,0),1);
    P(2) = 0;
    N(:,t,i) = binornd(1,P); 
end
end
N2 = N([1, 3:6],:,:);
C2 = C([1, 3:6],[1, 3:6]);
max(eig(C2))
e/c

% Figure:
figure(1)
clf
set(gcf,'unit','cent','position',[20 2 9, 25])

subplot(3,1,1)
hold on
for i = 1:length(Xc)
plot(Xc(i),Yc(i),'ko','markersize',20*S(i)^0.5/pi,'markerfacecolor','k')
end
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xlim',[-2.25 2.25],'ylim',[-2.25 2.25])
axis square

subplot(3,1,2)
bar([L_p(:),Val(:)],'grouped','facecolor',[0.5 0.5 0.5])
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xlim',[0 7],'ylim',[0 0.5])
xlabel('Patch','fontsize',10)
ylabel('Relative patch value','fontsize',10)
axis square

subplot(3,1,3)
hold on
for i = 1:Nsim
    plot(mean(Nbase(:,:,i),1),'k-','color',[0.8 0.8 0.8])
    plot(mean(N2(:,:,i),1),'k-','color',[0.5 0.5 0.5])
    plot(mean(N6(:,:,i),1),'k-')
end
set(gca,'ylim',[0 1])
set(gca,'tickdir','out','ticklength',[0.02 0.02])
xlabel('Time (y)','fontsize',10)
ylabel('Mean occupancy','fontsize',10)
set(gca,'xtick',0:20:100,'ytick',0:0.2:1)
axis square


