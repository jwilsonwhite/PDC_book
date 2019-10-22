function MPA_persistence(doSims1,doSims2,doSims3)

% Recreate the isoclines of persistence from Botsford et al. 2001

Rw = [0.1:0.1:2.5]; % reserve width
Fr = [0.01 0.1:0.025:1]; % fraction in reserves


 % Scenario 1: No advection, vary CRT
Sig = 10;
Mu = 0;

CRTs = [0.195 0.345 0.595];
E = zeros(length(Fr),length(Rw),length(CRTs));
if doSims1
for f = 1:length(Fr)
    for r = 1:length(Rw) 
        for c = length(CRTs):-1:1
        
        % Assemble coastline
        [X, Sig_tmp] = assemble_coastline(Fr(f),Rw(r),Sig);
        
        % Assemble connectivity matrix:
        D = dispersal_matrix(X,Mu,Sig_tmp);
        
        % Do persistence calculation:
        if c < 3 && length(CRTs) > 1
            if E(f,r,c+1) == 1 % if it was persistent for the more stringent criterion (save computing time)
                E(f,r,c) = 1;
            else
        E(f,r,c) = persistence(D,X,CRTs(c));
            
            end
        else  
        try    
        E(f,r,c) = persistence(D,X,CRTs(c));
        catch
            keyboard
        end
        end
        
        
        
        end
    end % end Rw
end % end Fr
save MPA_persist_sims1.mat
else
    load MPA_persist_sims1.mat Rw Fr E
end

% Figure panel 1:
figure(1)
clf
sh(1) = subplot(3,1,1);
hold on

for c = 1:length(CRTs)
    contour(Rw,Fr,E(:,:,c),[1 Inf],'k-')
end

ylabel('Fraction of coastline in reserves')
xlabel('Reserve width (relative to dispersal distance')

%--------------------------------------
% LOOP over advection (Mu) holding CRT constant
Rw = [0.1:0.5:3]; % reserve width
Fr = [0.01 0.1:0.1:1]; % fraction in reserves

Sig = 10;
Mus = [0,5,10];

CRT = 0.345;
E = zeros(length(Fr),length(Rw),length(Mus));
if doSims2
for f = 1:length(Fr)
    for r = 1:length(Rw) 
        for c = length(Mus):-1:1
            
        % Assemble coastline
        [X, Sig_tmp,Fact] = assemble_coastline(Fr(f),Rw(r),Sig,'semi');
        
        % Assemble connectivity matrix:
        D = dispersal_matrix(X,Mus(c)*Fact,Sig_tmp,'semi');
        
        % Do persistence calculation:
        if c < 3 && length(Mus) > 1
            if E(f,r,c+1) == 1 % if it was persistent for the more stringent criterion (save computing time)
                E(f,r,c) = 1;
            else
        E(f,r,c) = persistence(D,X,CRT);
            
            end
        else  
        try    
        E(f,r,c) = persistence(D,X,CRT);
        catch
            keyboard
        end
        end
        
        end
    end % end Rw
end % end Fr
save MPA_persist_sims2.mat
else
    load MPA_persist_sims2.mat Rw Fr E
end

% Figure panel 2:
sh(2) = subplot(3,1,2);
hold on

for c = 1:length(Mus)
    contour(Rw,Fr,E(:,:,c),[1 Inf],'k-')
end

ylabel('Fraction of coastline in reserves')
xlabel('Reserve width (relative to dispersal distance')

%--------------------------------------
% LOOP over FLEP but Ladv = 0, CRT constant
Sig = 10;
Mu = 0;

CRT =  0.345;
FLEPs = [0 0.1 0.2];
E = zeros(length(Fr),length(Rw),length(FLEPs));
if doSims3
for f = 1:length(Fr)
    for r = 1:length(Rw) 
        for c = length(FLEPs):-1:1
        
        % Assemble coastline
        [X, Sig_tmp] = assemble_coastline(Fr(f),Rw(r),Sig);
        
        % Assemble connectivity matrix:
        D = dispersal_matrix(X,Mu,Sig_tmp);
        
        % Do persistence calculation:
        if c < 3 && length(FLEPs) > 1
            if E(f,r,c+1) == 1 % if it was persistent for the more stringent criterion (save computing time)
                E(f,r,c) = 1;
            else
        E(f,r,c) = persistence(D,X,CRT,FLEPs(c));
            
            end
        else  
        try    
        E(f,r,c) = persistence(D,X,CRT,FLEPs(c));
        catch
            keyboard
        end
        end
        
        
        
        end
    end % end Rw
end % end Fr
save MPA_persist_sims3.mat
else
    load MPA_persist_sims3.mat Rw Fr E
end

% Figure panel 1:
sh(3) = subplot(3,1,3);
hold on

for c = 1:length(CRTs)
    contour(Rw,Fr,E(:,:,c),[1 Inf],'k-')
end

ylabel('Fraction of coastline in reserves')
xlabel('Reserve width (relative to dispersal distance')

set(sh(:),'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015],...
    'xlim',[0 2.2],'xtick',0:0.5:3,'ylim',[0 1],'ytick',0:0.2:1)
axis square

% -------------------------------------
function [X,Sig_tmp,Fact] = assemble_coastline(Fr,Rw,Sig,type)

if ~exist('type','var')
    type= 'inf';
end


        if Fr > 0 && Rw > 0
            
            Res = Sig*Rw;
            Frac = 1/Fr;
            Facts = 1:20;
            Rem = rem(Frac*Facts,1); % remainder after multiplication by a factor
            Fact = Facts(find(Rem==min(Rem),1)); % find lowest factor that provides integer
            
            try
        X = [ones(1,round(Res*Fact)), zeros(1, round(Res*Fact/Fr - Res*Fact))];
            catch
                keyboard
            end
        Sig_tmp = Sig*Fact;
        else
        X = 0;
        Sig_tmp = Sig;
        end
        
        if length(X) == 0
            keyboard
        end
        
        switch type
            case 'semi'
                X = repmat(X,[1,round(50/Fact)]);
        end
        
% -------------------------------------
function D = dispersal_matrix(X,Mu,Sig,type)

if ~exist('type','var')
    type= 'inf';
end

x = length(X);
XX = repmat(1:x,[x,1]);
Dist = XX - XX';

D = zeros(x);

% Loop a few times to simulate infinite (or semi-infinite) coastline:
switch type
    case 'inf'
for i = -50:50
D = D + normcdf(Dist+0.5+x*i,Mu,Sig) - normcdf(Dist-0.5+x*i,Mu,Sig);    
end
    case 'semi'

D = normcdf(Dist+0.5,Mu,Sig) - normcdf(Dist-0.5,Mu,Sig);    
       
end

% Helper function to calculate persistence
%---------------------------------------
function E = persistence(D,X,CRT,FLEP)
if ~exist('FLEP','var'); FLEP = 0; end

        Xm = repmat(X,[length(X),1]);
        Xm(Xm==0) = FLEP; % non-reserves have egg production = FLEP. Reserves have egg prod = 1
        
        C = D.*Xm./CRT;
        

E = double(max(eig(C))>1);  
            