function Taylor_model

% Discrete-time version of Taylor's (1979) model
% Used to illustrate effect of different demographic parameters on the
% importance of transient dynamics.

% Things to vary:
% - survival rate
% - overall fecundity
% - std of spawning age distribution
% - mode of spawning age distribution

% Figures: lambda1 & rho as a function of survival vs. a, overall fec vs.
% a, std vs. a

maxA = 20;
A = 1:maxA;
Fec1 = 30;

N = 20;
Survs = linspace(0.3,0.9,N);
Fecs = linspace(20,50,N);
Modes = linspace(5,maxA,N);
Vars = linspace(0,5,N);


% Setup figure
figure(1)
clf
set(gcf,'units','cent','position',[30 10 21 10])
colormap(1-gray)

% Left-hand: vary survival & mode of spawning distribution

for s = 1:N
    for a = 1:N
        
        Var = 5; % variance of spawning age distribution
        
        % Assemble matrix:
        Surv = ones(1,maxA-1).*Survs(s);
       % Fec1 = 10;
        
        % Hold variance constant:
        % mode = (k-1)q
        % var = kq^2
        % V = kq^2
        % mode = (Vq^-2 - 1)q = V/q - q
        % mq = V - q^2
        % q^2 - mq + V = 0
        % q = ( m +/- sqrt(m^2 - 4V))/2
        Th = max([ (Modes(a) + sqrt(Modes(a)^2 - 4*Var))/2, (Modes(a) - sqrt(Modes(a)^2 - 4*Var))/2]);
        k = Modes(a)/Th + 1;
        Fec = gampdf(A,k,Th);
        Fec = Fec/sum(Fec)*Fec1;
        
        L = diag(Surv);
        L = [L, zeros(length(Surv),1)];
        L = [Fec;L];
        
        Lam = eig(L);
        Lam = sort(Lam,'descend');
        L1(s,a) = Lam(1);
        Rho(s,a) = Lam(1)/abs(Lam(2));
        
    end % end a
end % end s

subplot(2,3,1)
contourf(Modes,Survs,(L1))
ylabel('Survival','fontsize',14)
%caxis([0 0.5])
%set(gca,'ytick',0:5)
format_axis(gca)
ch=colorbar;
set(ch,'tickdir','out','ticklength',0.02)
%set(ch,'ytick',log10(1:3),'yticklabel',1:3)

%ylabel(ch,'log10 Lambda')

subplot(2,3,4)
contourf(Modes,Survs,(Rho))
ylabel('Survival','fontsize',14)
%caxis([0.3 2])
format_axis(gca)
ch=colorbar;
set(ch,'tickdir','out','ticklength',0.02)
%set(ch,'ytick',log10([1:9 10:10:100]),'yticklabel',[1:9 10:10:100])

%----------------------------------------------------------
% Middle: vary fec & mode of spawning distribution

for s = 1:N
    for a = 1:N
        
        Var = 5; % variance of spawning age distribution
        Surv = 0.5;
        
        
        % Assemble matrix:
        Surv = ones(1,maxA-1).*Surv;
        
        % Hold variance constant:
        % mode = (k-1)q
        % var = kq^2
        % V = kq^2
        % mode = (Vq^-2 - 1)q = V/q - q
        % mq = V - q^2
        % q^2 - mq + V = 0
        % q = ( m +/- sqrt(m^2 - 4V))/2
        Th = max([ (Modes(a) + sqrt(Modes(a)^2 - 4*Var))/2, (Modes(a) - sqrt(Modes(a)^2 - 4*Var))/2]);
        k = Modes(a)/Th + 1;
        Fec = gampdf(A,k,Th);
        Fec = Fec/sum(Fec)*Fecs(s);
        
        L = diag(Surv);
        L = [L, zeros(length(Surv),1)];
        L = [Fec;L];
        
        Lam = eig(L);
        Lam = sort(Lam,'descend');
        L1(s,a) = Lam(1);
        Rho(s,a) = Lam(1)/abs(Lam(2));
        
    end % end a
end % end s

subplot(2,3,2)
contourf(Modes,Fecs,(L1))
xlabel('Modal age of spawning (y)','fontsize',14)
ylabel('Fecundity','fontsize',14)
title('log(\lambda)','fontsize',14)
%caxis([0 0.5])
%set(gca,'ytick',0:5)
format_axis(gca)
ch=colorbar;
set(ch,'tickdir','out','ticklength',0.02)
%set(ch,'ytick',log10(1:3),'yticklabel',1:3)


subplot(2,3,5)
contourf(Modes,Fecs,(Rho))
xlabel('Modal age of spawning (y)','fontsize',14)
ylabel('Fecundity','fontsize',14)
title('log(\rho)','fontsize',14)
format_axis(gca)
%caxis([0.3 2])
ch=colorbar;
set(ch,'tickdir','out','ticklength',0.02)
%set(ch,'ytick',log10([1:9 10:10:100]),'yticklabel',[1:9 10:10:100])



%----------------------------------------------------------
% Right: vary variance & mode of spawning distribution

for s = 1:N
    for a = 1:N
        
        Var = 10.^Vars(s); % variance of spawning age distribution
        Surv = 0.5;
        
        % Assemble matrix:
        Surv = ones(1,maxA-1).*Surv;
        
        % Hold variance constant:
        % mode = (k-1)q
        % var = kq^2
        % V = kq^2
        % mode = (Vq^-2 - 1)q = V/q - q
        % mq = V - q^2
        % q^2 - mq + V = 0
        % q = ( m +/- sqrt(m^2 - 4V))/2
        Th = max([ (Modes(a) + sqrt(Modes(a)^2 + 4*Var))/2, (Modes(a) - sqrt(Modes(a)^2 + 4*Var))/2]);
        Th = real(Th);
        k = Modes(a)/Th + 1;
        Fec = gampdf(A,k,Th);
        Fec = Fec/sum(Fec)*Fec1;
        
        L = diag(Surv);
        L = [L, zeros(length(Surv),1)];
        L = [Fec;L];
        
        Lam = eig(L);
        Lam = sort(Lam,'descend');
        L1(s,a) = Lam(1);
        Rho(s,a) = Lam(1)/abs(Lam(2));
        
    end % end a
end % end s


subplot(2,3,3)
contourf(Modes,sqrt(Vars),(L1))
ylabel('Standard deviation','fontsize',14)
%caxis([0 0.5])
set(gca,'ytick',0:5)
format_axis(gca)
ch=colorbar;
set(ch,'tickdir','out','ticklength',0.02)
%set(ch,'ytick',log10(1:0.1:2),'yticklabel',1:0.1:2)


subplot(2,3,6)
contourf(Modes,sqrt(Vars),(Rho))
ylabel('Standard deviation','fontsize',14)
format_axis(gca)
set(gca,'ytick',0:5)
%caxis([0.3 2])
ch=colorbar;
set(ch,'tickdir','out','ticklength',0.02)
%set(ch,'ytick',log10([1:9 10:10:100]),'yticklabel',[1:9 10:10:100])




function format_axis(A)
set(A,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k',...
    'xtick',5:5:20)
axis square



