function [x,N2] = urchin_IPM

% Make at IPM that matches the red urchin von Foerster model

% parameter definitions:
%[Linf      k   x0  M  Lfish Lmat Lvar];
fixparm = [112 0.3 0 NaN 0 0 10.1]; 
Rmu = 60; % mean recruitment

% IPM integration mesh parameters
meshsize = 200;
meshmin = 0;
meshmax = fixparm(1)*2;
x = linspace(meshmin,meshmax,meshsize);
meshdiff = diff(x(1:2));
edges = x - meshdiff/2;
dy=x(2)-x(1);
Sy=makeSimpVec(dy,meshsize);

% this is useful for integrating the kernel later
Symat= repmat(Sy(:)',[length(Sy),1]);

% Recruit size distribution
Rmean = 5;
Rsd = dy;
Rvec = normpdf(x,Rmean,Rsd);
Rvec = Rvec./(Sy*Rvec');


% Create the kernel:
kmat = kernmatSimp(x,0,fixparm,1/36);

% Simpson's integration:
kmat = Symat.*kmat;

% Initial size distribution
N = zeros(size(kmat,1),100);
N(:,1) = Rvec.*Rmu;

% Run with only one cohort
for t = 2:100
N(:,t) = kmat*N(:,t-1);
end

% Run with only annual cohorts
N2 = N;
for t = 2:1e4
N2(:,t) = kmat*N2(:,t-1) + Rmu*Rvec(:);
end

% Plot the distribution, for insertion onto the earlier figure with the
% M-vF model result
N2(x<=17.5,:) = NaN;
N2 = N2./(repmat(nansum(N2),[length(x),1])*dy); % get the distribution

doplot =  false;
if doplot
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 6])
hold on

plot(x,N2(:,end));
set(gca,'xlim',[0 150],'xtick',0:20:160,'ylim',[0 0.06])
end
   
% Now create a kernel that also has reproduction:
Fec1 = 1; Fec2 = 1.5;
Rmean = 10; Rsd = 5;
fixparm(8:11) = [Rmean, Rsd, Fec1, Fec2];
kmat2 = kernmatSimp(x,0,fixparm,1);

if doplot
figure(2)
clf
set(gcf,'units','cent','position',[30 10 9 6])
hold on
contourf(flipud(kmat2))
%mesh(flipud(kmat2))
colormap(1-gray)
caxis([-0.01 0.07])
axis square
%set(gca,'view',[-25, 72])
end


