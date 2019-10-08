function ricker_cobweb

% Plot the Ricker function with cobwebbing

%----------------------------------------------------------------
% Plot #1 : Example (Figure 2.9):
% Ricker parameters
a = 5;
b = 1;

Af = 0.1; %arrowhead factor, to scale position of arrowhead

% The plot will have four panels:
% Upper left will show all lines to axes
% Lower left will just show shortcut
% Right hand will show time series

maxN = 10;
N0 = 0:0.01:maxN;
N1 = ricker(a,b,N0);

figure(1)
set(gcf,'units','cent','position',[30,40,18,15])
clf

% First plot the function & replacement line
subplot(2,2,1)
hold on
plot(N0,N1,'k-') % function
plot([0 maxN],[0,maxN],'k--') % replacement line

% Now run the model:
T = 10;
Nt(1) = 2.5;
Ft(1) = ricker(a,b,Nt(1));
for t = 2:T
    Nt(t) = Ft(t-1);
    Ft(t) = ricker(a,b,Nt(t));
end

% cobweb:
do_cobweb(Nt,Ft,2,Af,6)

% Plot square:
plot([0, Nt(2)],[Nt(2),Nt(2)],'k:')
plot([Nt(2), Nt(2)],[0,Nt(2)],'k:')

% Format plot:
ylim([0 max(N1)*1.1])
axis equal
xlim([0 3])
ylim([0 3])
set(gca,'xtick',[],'ytick',[])
xlabel('N_t','fontsize',12); 
ylabel('N_t_+_1','fontsize',12);

% Top-right: time series
subplot(2,2,2)
plot(0:2,Nt(1:3),'ko','linestyle','-','color',[0.5 0.5 0.5],'markeredgecolor','k')

ylim([0 3])
xlim([-0.5 10])
set(gca,'xtick',[],'ytick',[])
xlabel('Time','fontsize',12);
ylabel('N','fontsize',12);


% Lower left - full cobweb:
subplot(2,2,3)
hold on
plot(N0,N1,'k-') % function
plot([0 maxN],[0,maxN],'k--') % replacement line

do_cobweb(Nt,Ft,10,Af,eps) % the full thing

% Add points:
plot(Nt,Ft,'ko','markerfacecolor','none')

xlim([0 3])
ylim([0 3])
set(gca,'xtick',[],'ytick',[])
xlabel('N_t','fontsize',12); 
ylabel('N_t_+_1','fontsize',12);

% Lower right: full timeseries
subplot(2,2,4)
hold on
plot(0:T-1,Nt,'ko','linestyle','-','color',[0.5 0.5 0.5],'markeredgecolor','k')

ylim([0 3])
xlim([-0.5 10])
set(gca,'xtick',[],'ytick',[])
xlabel('Time','fontsize',12);
ylabel('N','fontsize',12);

% End set of plots #1
%----------------------------------------------------------------

%----------------------------------------------------------------
% Plots #2: Examples of Ricker stability analysis, with different values of
% F'(X) (Figure 2.12)

figure(2)
set(gcf,'units','cent','position',[20,40,18,28])
clf

% Ricker parameters.
% Behavior depends on F'(X), which is equal to 1 - log(a)
A = [exp(0.85), exp(1.8), exp(2.0), exp(2.2)]; % stable, dampened, neutral limit cycles, undamp osc
b = 1;

Ninit = [2 1 1.5e-1 2e-1];

% subplot coordinates:
SP = [1 2; 3 4; 5 6; 7 8];

for aa = 1:length(A)

a = A(aa);
maxN = 10;
N0 = 0:0.01:maxN;
N1 = ricker(a,b,N0);

% Now run the model:
T = 10;
Nt(1) = log(a)/b + Ninit(aa);
Ft(1) = ricker(a,b,Nt(1));
for t = 2:T
    Nt(t) = Ft(t-1);
    Ft(t) = ricker(a,b,Nt(t));
end

% Left - cobweb:
subplot(4,2,SP(aa,1))
hold on
plot(N0,N1,'k-') % function
plot([0 maxN],[0,maxN],'k--') % replacement line

do_cobweb(Nt,Ft,T,Af,eps) % the full thing

xlim([0 5])
ylim([0 max(N1)*1.2])
%set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',[],'ytick',[])
if aa == 4; xlabel('N_t','fontsize',12); end
ylabel('N_t_+_1','fontsize',12);

% Right: full timeseries
subplot(4,2,SP(aa,2))
hold on
plot(0:T-1,Nt,'ko','linestyle','-','color',[0.5 0.5 0.5],'markeredgecolor','k')

ylim([min(Nt)*0.8 max(Nt)*1.1])
xlim([-0.5 10])
%set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',[],'ytick',[])
if aa == 4; xlabel('Time','fontsize',12); end
ylabel('N','fontsize',12);

end % end loop

% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = ricker(a,b,N)
N = a*N.*exp(-b.*N);

function do_cobweb(Nt,Ft,T,Af,Ms)

% Ms is markersize

% Initial line from origin to Ft(1):
plot([Nt(1) Nt(1)],[0 Ft(1)],'k-')
plot(Nt(1),Ft(1)-Af,'marker','^','markerfacecolor','k','markeredgecolor','k','markersize',Ms)
for t = 1:T-1
    
   % Segment from Ft to replacement line
   plot([Nt(t),Ft(t)],[Ft(t),Ft(t)],'k-')
   if Nt(t) > Ft(t) % if the cobweb line is pointing to the left
   plot(Ft(t)+Af,Ft(t),'marker','<','markerfacecolor','k','markeredgecolor','k','markersize',Ms)
   else % pointing to the right
   plot(Ft(t)-Af,Ft(t),'marker','>','markerfacecolor','k','markeredgecolor','k','markersize',Ms)    
   end
   
   % Segment from replacement line to curve:
   plot([Nt(t+1),Nt(t+1)],[Ft(t),Ft(t+1)],'k-')
   if Ft(t+1) > Nt(t+1) % line is point up:
   plot(Nt(t+1),Ft(t+1)-Af,'marker','^','markerfacecolor','k','markeredgecolor','k','markersize',Ms)
   else % line is pointing down:
   plot(Nt(t+1),Ft(t+1)+Af,'marker','v','markerfacecolor','k','markeredgecolor','k','markersize',Ms) 
   end
 
end


