function instarfig

% Make plots showing insect instars of varying length

% Based on data from Anopheles gambiae (Bayoh & Lindsay 2003)

% Two temps (18 & 32)

% Instars:
Instar = 1:6;

% Times (cumulative):
T18 = [0 4 7 11 14 19 22];
T32 = [0 1 4 5 7 8 9];
T = [T18; T32];

d = 0.5;
Offset = [-0.1 0.1];

% Figure:
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 9])
hold on

for t = 1:2
   for i = 1:length(Instar)
      patch([T(t,i), T(t,i),T(t,i+1),T(t,i+1)],...
            [Instar(i), Instar(i)+d, Instar(i)+d, Instar(i)]+Offset(t),....
            [0.5 0.5 0.5])
        if i < length(Instar)
      plot([T(t,i+1),T(t,i+1)],...
          [Instar(i)+d/2,Instar(i+1)+d/2]+Offset(t),'k-')
        end
   end 
end

set(gca,'tickdir','out','ticklength',[0.02 0.02],'fontsize',12)
set(gca,'xtick',0:5:30,'ytick',[])
xlabel('Age (d)','fontsize',14)
ylabel('Stage','fontsize',14)
set(gca,'xlim',[0 25],'ylim',[0.5 7])


