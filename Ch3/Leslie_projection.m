function Leslie_projection

% Project Leslie matrices to illustrate their behavior (Fig. 3.5)

% Assemble the matrix (Table 3.2):
F = [0 0 0 0 0.8 1 1 1 0.1]; % fecundity
S = [0.7 0.7 0.9 0.9 0.9 0.9 0.3 0.3]; % survival
Smat = [diag(S), zeros(length(S),1)]; % lower half of matrix
L = [F;Smat];  % the full matrix

% Initial conditions:
N0 = [1000, 600, 400, 300, 200, 0 0 0 0]';

% Project:
T = 61;
N = nan(length(N0),T);
N(:,1) = N0;
for t = 2:T
    N(:,t) = L*N(:,t-1);
end

% Now for a declining population:
L2 = L;
L2(1,5:8) = 0.7;
% Project:
T = 61;
N2 = nan(length(N0),T);
N2(:,1) = N0;
for t = 2:T
    N2(:,t) = L2*N2(:,t-1);
end

% Plot:
figure(1)
clf
set(gcf,'units','cent','position',[10 10 18 15])
Col = 1-repmat(linspace(0.15,1,9)',[1,3]);

subplot(2,2,1)
hold on
b = bar(N(:,1:5)','grouped','edgecolor','none');
for i = 1:length(b)
    b(i).FaceColor=Col(i,:);
end
set(gca,'xtick',1:5,'xticklabel',0:4,'ytick',0:200:1000)
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'ylim',[0 1010])
    

subplot(2,2,2)
hold on
b=bar(N(:,:)','stacked','edgecolor','none');
xlim([0 T])
for i = 1:length(b)
    b(i).FaceColor=Col(i,:);
end
set(gca,'xtick',(0:20:60)+1,'xticklabel',0:20:60,'ytick',0:2000:10000)
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'ylim',[0 6000],'xlim',[1 61])

subplot(2,2,3)
hold on
b=bar(N2(:,1:5)','grouped','edgecolor','none');
for i = 1:length(b)
    b(i).FaceColor=Col(i,:);
end
set(gca,'xtick',1:5,'xticklabel',0:4,'ytick',0:200:1000)
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'ylim',[0 1010])

subplot(2,2,4)
hold on
b=bar(N2(:,:)','stacked','edgecolor','none');
for i = 1:length(b)
    b(i).FaceColor=Col(i,:);
end
xlim([0 T])
set(gca,'xtick',(0:20:60)+1,'xticklabel',0:20:60,'ytick',0:1000:5000)
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'ylim',[0 2600],'xlim',[1 61])

max(eig(L))
max(eig(L2))



