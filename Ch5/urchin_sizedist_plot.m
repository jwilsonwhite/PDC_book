function urchin_sizedist_plot

% Plot figure from Botsford et al. (1994) showing Tegner's bimodal size
% distribution

% Read in data:
D= importdata('Urchin_data_tegner.csv');
D = D.data;
U = D(1,2)-D(:,2);
U = U(2:end);

% Plotting
figure(2)
clf
set(gcf,'units','cent','position',[10 10 9 6])
hold on

% length bins
L = 6.25:2.5:148.75;

% standardize to integrate to 1
U = U./sum(U)./diff(L(1:2));

bar(L,U,'facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5])

% Fit a model.
% Parameters:
k = 0.3; % vB growth rate
Linf = 112; % asymptotic size
sig_Linf = 10.1; % SD of Linf
A = 0:50; % age classes
L = 0:0.1:150; % length classes
Rmu = 60; % mean recruitment 

Nc = 1e4; % Number of simulated curves
Linfs = normrnd(Linf,sig_Linf,[Nc,1]);

% Iterate size-structured model with variable Linf and constant recruitment
for i = 1:Nc
p = [k, Linfs(i), 0, NaN, Rmu];
N(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case6');
end

N(:,L<=17.5) = NaN; % truncate small, unobserved sizes

% standardize to integrate to 1
N = N./repmat(nansum(N,2),[1,length(L)]);
% Plot mean of size structured model results
plot(L,nanmean(N)./nansum(nanmean(N))./diff(L(1:2)),'k-','linewidth',1)

set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k','fontsize',10)
ylim([0 0.06])
set(gca,'xtick',0:10:150,'xlim',[0 155])
ylabel('Frequency','fontsize',12)
xlabel('Test diameter (mm)','fontsize',12)

% Add IPM result:
[x,N2]=urchin_IPM;
plot(x,N2(:,end),'r-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = SizeStruct_VarR_EqVB(t,p,A,L,mortcase)
% size-based von Foerster equation
% Pulsed recruitment, no variability in von Bertalanffy equation

k = p(1);
Linf = p(2);
mu = p(3);
sd = p(4);
Rmu = p(5);

% von Bertalanffy growth
L0 = 0;
A=log(1-L/Linf)./-k; % inverse of vB equation
A(L>Linf) = 0;
A(A==0) = max(A); % substitute the maximum value


gl0 = k*(Linf-L0);
gl =  k*(Linf-L);

if isnan(sd) % is recruitment pulsed or constant?
% Constant recruitment
R = Rmu;
else
% Pulsed recruitment. Mean = Jun 30. Time is fractions of a year
tt = mod(t-A,1); % take fraction of the year for all ages going into the past
R = normpdf(tt,mu,sd);
dA = diff(A);
dA = [dA(1) dA]; % pad with a dummy first entry to extend length
R = R./mean(R);
R = R.*Rmu;
end

% Mortality rate: (several different options)
switch mortcase
    case 'case1'
% unimodal mortality rate:
D = 0.60 - 0.59/50^2*(L-50).^2; % min at 0.01, peak of 0.6 at 50 mm 
D = max(D,0.01);
    case 'case2'
% constant mortality rate:
D = 0.1;
    case 'case3'
% low mortality for small individuals:
D = ones(size(L)).*0.1;
D(L<40) = 0;
    case 'case4'
% high mortality for small individuals:
D = ones(size(L)).*0.1;
D(L<80) = 1;
    case 'case5'
% high mortality for medium individuals:
D = ones(size(L)).*1.0;
D(L<40) = 0;
D(L>80) = 0.1;
    case 'case6'
 % Declining mortality:
D = 5.81.*exp(-0.043.*L);
end

% Integrate D/g(l) wrt l
Int = cumsum(D./gl).*diff(L(1:2));

% Solution:
N = R.*gl0./gl.*exp(-Int);
%N(L>0.99*Linf)=0;
N(N<0) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
