function logistic_versions

% Illustrate differences in the continuous vs. discrete versions of the
% logistic map

Rs = exp([0.5 1.5 2]); % growth rates
K = 10; % carrying capacity
N0 = 0.1; % initial conditions

tspan = [0 50]; % time interval for solution

% Setup figure
figure(3)
set(gcf,'units','cent','position',[40,40,8,15])
clf

% Plot solution for each value of R
for i = 1:length(Rs)
    
subplot(3,1,i)
hold on
R = Rs(i);

% Continuous version:
[TT,N] = ode45(@(t,N) Logistic(t,N,log(R),K),tspan,N0);
ph=plot(TT,N,'k-');


N1 = zeros(tspan(2),1);
N1(1) = N0;
N2 = N1;
% Discrete version:
for t = 2:(tspan(2)+1)
    
    % Map 1:
    N1(t) = N1(t-1).*exp(log(R).*(1-N1(t-1)/K));
    
    % Map 2:
    N2(t) = N2(t-1)*(1 + log(R).*(1 - N2(t-1)/K));
   
    
end % end loop over time

% Create plots
plot(tspan(1):tspan(2),N1,'ko','linestyle','-')
plot(tspan(1):tspan(2),N2,'k^','linestyle','-','color',[0.5 0.5 0.5])

% Format plots
set(gca,'ylim',[0 15],'ytick',0:5:15,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:5:50,'xlim',[0 30])
if i == 2
ylabel(gca,'Population abundance (N)','fontsize',12)
end
if i == 3
xlabel(gca,'Time','fontsize',12)
end

end

% Define the continuous time model
function dN = Logistic(t,N,r,K)
% Exponential growth
dN = r.*N.*(1-N./K);