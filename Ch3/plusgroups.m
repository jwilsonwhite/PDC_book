function plusgroups

% Matrix models with plus groups:

% Example 1: Hebblewith et al, black bears

L1 = [0 0 0 0 0.39;...
      0.64 0 0 0 0;...
      0 0.67 0 0 0;...
      0 0 0.765 0 0; ...
      0 0 0 0.765 0.835];
 
  
 % Full model. What if fecundity increased with age
 A = 21;
 F = zeros(1,A);
 F(5:end) = 0.39;
 
 S = [0.64, 0.67, 0.765, 0.765, 0.835*ones(1,A-5)];
 L2 = [diag(S),zeros(A-1,1)];

 L2 = [F;L2];

 
% Elasticities:
% first get the dominant eigenvectors:
% model 1:
[W,Ls,V]=eig(L1);
Ls = diag(Ls);
maxL = Ls==max(Ls);
W = W(:,maxL);
W1 = W;
V = V(:,maxL);
S = L1(2:end,1:end-1);
S = diag(S);
l = Ls(maxL);

% fecundity elasticities:
eF1 = V(1).*W'.*L1(1,:)./l./dot(V,W);

% survival elasticities:
eS1 = V(2:end).*W(1:end-1).*S(:)./l./dot(V,W);

% plus-group:
eP1 = V(end)*W(end).*L1(end,end)./l./dot(V,W);

E1 = [eF1(:);eS1(:);eP1];

% model 2:
[W,Ls,V]=eig(L2);
Ls = diag(Ls);
maxL = Ls==max(Ls);
W = W(:,maxL);
W2 = W;
V = V(:,maxL);
S = L2(2:end,1:end-1);
S = diag(S);
l = Ls(maxL);

% fecundity elasticities:
eF2 = V(1).*W'.*L2(1,:)./l./dot(V,W);

% survival elasticities:
eS2 = V(2:end).*W(1:end-1).*S(:)./l./dot(V,W);

E2 = [eF2(:);eS2(:)];

% Make the figure:
figure(1)
clf
set(gcf,'units','cent','position',[2 30 9 16])

% elasticity for plus-group model
subplot(2,1,1)
set(gca,'position',[0.15 0.56 0.64 0.35])

hold on
bar(E1,'facecolor',[0.6 0.6 0.6])
set(gca,'color','none',...
        'ytick',[0:0.2:1],'ylim',[0 1],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xlim',[0 11],'xtick',[])
ylabel('Elasticity','fontsize',12)
axes('position',get(gca,'position'))
plot(cumsum(E1),'k-')
set(gca,'yaxislocation','right',...
        'color','none',...
        'ytick',0:0.2:1,'ylim',[0 1],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xlim',[0 11],'xtick',1:10,'xticklabel',[1:5,1:5])
ylabel('Cumulative elasticity','fontsize',12)
xlabel('Age (a)','fontsize',12)

subplot(2,1,2)
set(gca,'position',[0.15 0.11 0.64 0.35])
hold on
bar(E2,'facecolor',[0.6 0.6 0.6])
set(gca,'color','none',...
        'ytick',[0:0.02:1],'ylim',[0 0.1],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xlim',[0 43],'xtick',[])
ylabel('Elasticity','fontsize',12)
axes('position',get(gca,'position'))
plot(cumsum(E2),'k-')
set(gca,'yaxislocation','right',...
        'color','none',...
        'ytick',0:0.2:1,'ylim',[0 1],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xlim',[0 43],'xtick',[5:5:45],'xticklabel',[5:5:20])
ylabel('Cumulative elasticity','fontsize',12)
xlabel('Age (a)','fontsize',12)
 

 
 
 % Compare transient dynamics:
 N01 = [1; zeros(4,1)];
 N02 = [1 zeros(1,A-1)]';
 N01 = W1./sum(W1);
 N02 = W2./sum(W2);
 
 T = 100;
  N1 = nan(5,T);
  N2 = nan(A,T);
  N1(:,1) = N01;
  N2(:,1) = N02;
  
 for t = 2:T
    N1(:,t)=L1*N1(:,t-1);
    N2(:,t)=L2*N2(:,t-1);
    if t == 10;
        N1(1:4,t) = 0;
        N2(1:4,t) = 0;
    end
        

 end
 
%% Uncomment below if a figure is desired:
% figure(2);
% clf
% hold on
% plot(1:T,(sum(N1)),'k--')
% plot(1:T,(sum(N2)),'k-')

% keyboard
     