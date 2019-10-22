function MPA_spatial_map_habitat

% Make a plot similar to the DPR paper showing effect of both habitat and
% MPA placement on pop density for short- & long-distance dispersers

% Load in data
D = csvread('Habitats&Reserves.csv');

OKhab = 1:155;
Hab = D(OKhab,1);
Res = D(OKhab,2);

% Run a simple population model to equilibrium
Dispdist = [2, 15];
alpha = 1/0.25;
LEP_out = 0.0;
T = 100;

for d = 1:2
    D = dispmatrix(length(Hab),'gaussian',Dispdist(d));
R = zeros(length(Hab),T);
R(:,1) = 1;
LEPvec = ones(length(Hab),1)*LEP_out;
LEPvec(Res==1) = 1;
for t = 2:T
   L(:,t) = R(:,t-1).*LEPvec; % reproduction
   S(:,t) = D*L(:,t); % dispersal
   R(:,t) = alpha*S(:,t)./(1 + S(:,t)); % density-dependent survival
   R(Hab==0,t) = 0; % no habitat = settlers die
end
R_out(:,:,d) = R;

end % end loop over d

% Plotting
figure
clf
hold on

Col = [repmat(0.2,[1,3]); repmat(0.7,[1,3])];
for j = 1:(length(Hab))
    
    % plot hab
    if Hab(j)
    patch([j-0.5 j-0.5 j+0.5 j+0.5],[0 Hab(j) Hab(j) 0]*-0.5,Col(1,:),'edgecolor','none')
    end
    
    % plot Res
    if Res(j)
    patch([j-0.5 j-0.5 j+0.5 j+0.5],-0.5+-0.5*([0 Res(j) Res(j) 0]),Col(2,:),'edgecolor','none')
    end
    
end

Lt = {'-','--'};
for d = 1:2
   plot(R_out(:,end,d),'k-','linestyle',Lt{d}) 
end

set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',0:20:160,'ytick',[])
set(gcf,'units','cent','position',[10 10 19 6])

keyboard