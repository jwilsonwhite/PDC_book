function cannibal

% Implement Tribolium cannibalism model, as in Hastings & Costatino (1987,
% Am Nat)

% (discrete-time version is given in book chapter)

% Parameters from H&C Fig 1:
mu_e = 0;
c = 0.005;

Aes = 4:2:8; % egg-stage durations
Bs = 0:1:1200; % birth rates
Als = 2:0.1:40; % Larval stage durations

% create grids of values on which to evaluate model
Bsmat = repmat(Bs(:),[1,length(Als)]);
Alsmat = repmat(Als(:)',[length(Bs),1]);

for a = 1:length(Aes)

% stability conditions, from the paper
omega = 2*pi./(Aes(a) + Alsmat);
Gam = omega.^2./(2 - cos(Aes(a).*omega) - cos(Alsmat.*omega));
Cond1 = Bsmat.*c.*exp(-mu_e.*Aes(a));
Cond2 = Gam.*exp(Gam.*Alsmat.*Aes(a));

Stab(:,:,a) = double(Cond1 < Cond2); % determine if it is stable or not at this point

end

% Plot stability contours
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 11])
hold on
for a = 1:length(Aes)
contour(Als,Bs,Stab(:,:,a),1,'k')
end

% Plot formatting
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',0:10:40,'ytick',0:400:1200,'fontsize',10)
ylabel(gca,'Birth rate (B, eggs per day)','fontsize',12)
xlabel(gca,'Duration of larval stage (AL, days)','fontsize',12)
axis square

