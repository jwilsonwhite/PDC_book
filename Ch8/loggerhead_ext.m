function loggerhead_ext(doruns)

% Simulate extinction probabilities for the loggerhead age & stage models
% Input argument 'doruns': logical, if true executes simulations, if false
% it reads in stored values.

figure(2)
clf

if doruns
rng(22) % set random seed to obtain consistent results

%--------------------------
% Set up the models

% Assemble the matrix (for loggerhead turtles, from Crowder et al. 1994):
% Stage model
F = [0 0 0 4.665 61.896]; % fecundity
D = [0 0.703 0.657 0.682 0.8091]; % diagonal (stage retention)
S = [0.675 0.047 0.019 0.061]; % survival
Smat = [diag(S), zeros(length(S),1)]; % lower half of matrix
L = [F;Smat]+diag(D);  % the full matrix


L = L/max(eig(L));

% Manipulate 3 different stages: egg, small juv, large juv
Ind = {2, [7 8], [13, 14]};

%Ind = 13; % entry for small adult survival (to be varied)
%Ind = [13 14];
%Ind = 2; % juveniles only
%Ind = [7 8];
%Ind = find(L>0 & L<1); 
%Ind = [4 5]; % entry for fecundity

[W,Ls,V]=eig(L);
Ls = diag(Ls);
maxL = Ls==max(Ls);
W = abs(W(:,maxL));
V = V(:,maxL);
l = max(Ls);


% Age model
% Stage durations:
Ts = [1 7 8 6 32];
% Survival
sig = [0.6747, 0.75, 0.6758, 0.7425, 0.8091];
Fec = [zeros(1,sum(Ts(1:3))), zeros(1,Ts(4))*F(4), ones(1,Ts(5))*76.5];
Surv = [ones(1,Ts(1))*sig(1),...
    ones(1,Ts(2))*sig(2),...
    ones(1,Ts(3))*sig(3),...
    ones(1,Ts(4))*sig(4),...
    ones(1,Ts(5))*sig(5)];
Surv = Surv(1:end-1);
LA = diag(Surv);
LA = [LA, zeros(length(Surv),1)];
LA = [Fec; LA];

% Manipulate 3 different stages: small adult, juvenile, fecundity
IndA = {2,... % egg
        [size(LA,1)*(1:7) + (3:9)],...    % sm juveniles
        [size(LA,1)*(9:15) + (11:17)]}; % large juvs

%IndA = size(LA,1)*(9:15) + (11:17); % small adults
%IndA = size(LA,1)*(1:7) + (3:9); % large juveniles
%IndA = find(LA>0 & LA<1);
%IndA = 1:54;
%IndA = 2;

LA = LA/max(eig(LA));

[WA,LAs,VA]=eig(LA);
LAs = diag(LAs);
maxLA = LAs==max(LAs);
WA = abs(WA(:,maxLA));
VA = VA(:,maxLA);
lA = max(LAs);

%---------------------------

%---------------------------
% Simulations
X = 0.1; % extinction threshold
T = 75; % time horizon
CV = [0:0.1:1]; % CVs of variability
n = 1e1; % how many simulations to make


% Random variability in the parameters, using a normal distribution:
for c = 1:length(CV)
Rands(c,:,:) = normrnd(0,repmat(CV(c),[1,T,n])); 
end % end loop over CVs

% Do simulations (N is the stage structured, NA is age structured)
% Initial conditions:
N = zeros(5,T,n);
N(:,1,:) = repmat(W/sum(W),[1,1,n]);
NA = zeros(size(LA,1),T,n);
NA(:,1,:) = repmat(WA/sum(WA),[1,1,n]);

N(:,1,:) = repmat( N(:,1,1) .* (1+0.75*(rand(size(L,1),1,1)-0.5)), [1, 1, n]);
NA(:,1,:) = repmat( NA(:,1,1) .* (1+0.75*(rand(size(LA,1),1,1)-0.5)), [1, 1, n]);

N = max(0,N);
NA = max(0,NA);

N(:,1,:) = N(:,1,:)./repmat(sum(N(:,1,:)),[size(N,1),1,1]);
NA(:,1,:) = NA(:,1,:)./repmat(sum(NA(:,1,:)),[size(NA,1),1,1]);

Results = struct([]);

else % else if the runs already exist
    load loggerhead_ext_results.mat
end % 


for s = 1:3 % loop over stages
for c = 1:length(CV)
    
    if doruns
    
    for nn = 1:n
        for t = 2:T
            
        
            
            % stage
            Lt = L;
            try
            Tmp = min(1,max(0,L(Ind{s})+L(Ind{s})*Rands(c,t,nn)));
            catch
                keyboard
            end
            Lt(Ind{s}) = Tmp;
            LA(2:end,:) = min(1,LA(2:end,:));
            N(:,t,nn) = Lt*N(:,t-1,nn);
            
            % age
            LAt = LA;
            Tmp = min(1,max(0,LA(IndA{s})+LA(IndA{s})*Rands(c,t,nn)));
            LAt(IndA{s}) = Tmp;
            LAt(2:end,:) = min(1,LAt(2:end,:));
            NA(:,t,nn) = LAt*NA(:,t-1,nn);
            
            
        end % end T
    end % end n

    
    Results(s).CV(c).N = N;
    Results(s).CV(c).NA = NA;
    
    else % if not doruns
        
        N = Results(s).CV(c).N;
        NA = Results(s).CV(c).NA;
        
    end % end if doruns
    
% Quasi-exinctions
NN = squeeze(sum(N,1));
Ex = sum(NN<=X)>0;
Px(c) = mean(Ex);

NNA = squeeze(sum(NA,1));
ExA = sum(NNA<=X)>0;
PxA(c) = mean(ExA);
        

% Plotting:
if any(c == [1,5,10])
figure(2)

if c == 1; SPi=1; end
if c == 5; SPi=2; end
if c == 10; SPi=3; end

SP = reshape(1:12,[3,4])';
SP = SP(2:end,:);

subplot(4,3,SP(SPi,s))
hold on

for i = 1:1
    plot(sum(NA(:,:,i)),'k','linewidth',1)
    plot(sum(N(:,:,i)),'k--','linewidth',1)
    
end
ylim([0 1.2])
xlabel('Time (y)')
ylabel('Population density (N)')
set(gca,'tickdir','out','ticklength',[0.02 0.02],'xcolor','k','ycolor','k')
set(gca,'ytick',0:0.2:2,'xtick',0:20:100)
end % end if plot

end % end CV

%----------------------------------
% Plotting
figure(2)

subplot(4,3,s)
hold on
plot(CV,Px,'k')
plot(CV,PxA,'k--')
xlabel('CV of variability')
ylabel('Probability of extinction')
legend('Stage','Age','Location','best')

end % end loop over stages

if doruns
save loggerhead_ext_results.mat -v7.3
end


