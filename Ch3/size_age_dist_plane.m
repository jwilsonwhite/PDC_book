function size_age_dist_plane

% Create figure illustrating the continuous-time, continuous-age approach
% to population modeling, and the idea of partial derivatives and
% characteristic lines

% Age and size vectors
Age = 0:0.1:12;
Time = 0:0.1:24;

% A toy population model that has some funky dynamics to give the surface a
% shape
N = nan(length(Age),length(Time));
N(:,1) = exp(-0.2*Age); % Initial age distribution
L1 = 0.02*(cos(Time/5)+1)+0.99; % Sinusoidal function (arbitrary)
L2 = 1; % growth factor, if desired
for t = 2:length(Time)
    N(:,t) = L1(t)*L2*N(:,t-1);
end

% Plotting
% Figure setup
figure(1)
clf
set(gcf,'unit','cent','position',[1,1,12,30])
% camera view angle
Az = -147; % Azimuth
El = 45; % Elevation

% Plot 1: plane at constant T
subplot(3,1,1)


%surf(fliplr(Time),Age,N,ones(size(N)),'facealpha',0.5)
surf(fliplr(Time),Age,N,'facealpha',0.5)
shading flat
hold on

% rotate 
view([Az,El])

% dN/dt (plane at constant A)
Vert = [0, 6, 0; 0, 6, 0.5*max(N(:)); max(Time), 6,0.5*max(N(:)); max(Time), 6, 0];
Fac = [1 2 3 4];
patch('Faces',Fac,'Vertices',Vert,'FaceColor',[0.8 0.8 0.8],'facealpha',1,'edgecolor','none')
plot3(fliplr(Time),ones(size(N,2))*6,N(Age==6,:),'k-')

% Set axis labels, etc
fig_cleanup(Age,Time)

% Plot 1: plane at constant A
subplot(3,1,2)

surf(fliplr(Time),Age,N,'facealpha',0.5)
shading flat
hold on
% rotate properly
view([Az,El])


% dN/dA (plane at constant T)
Vert = [5, 0, 0; 5, 0, max(N(:)); 5, max(Age),max(N(:)); 5, max(Age), 0];
Fac = [1 2 3 4];
patch('Faces',Fac,'Vertices',Vert,'FaceColor',[0.8 0.8 0.8],'facealpha',1,'edgecolor','none')
plot3(ones(size(N,1))*5,Age,N(:,Time==(max(Time)-5)),'k-')

% Set axis labels, etc
fig_cleanup(Age,Time)

% Plot 3: Characteristic lines
subplot(3,1,3)

%surf(fliplr(Time),Age,N,ones(size(N)),'facealpha',0.5)
surf(fliplr(Time),Age,N,'facealpha',0.5)
shading flat
hold on

% rotate properly
view([Az,El])

% Plot some characteristic lines (individual moving through age at a
% constant rate). This is more complicated than it looks because we are
% plotting time in reverse (so that the matrix looks right when plotted in
% 3D)
Is = [1, 5, 10];
for i = 1:length(Is)
    Index = find(Time==Is(i));
NN = N(:,Index:length(Age)+Index-1);
plot3(fliplr(Age)+(max(Time)-max(Age)-Is(i)),Age,diag(NN),'k-')
end

% Set axis labels, etc
fig_cleanup(Age,Time)

function fig_cleanup(Age,Time)
set(gca,'xtick',(0:6:24),'ytick',0:4:12,'ztick',0:2:10)
set(gca,'xticklabel',20:-5:0,'yticklabel',0:5:15,'zticklabel',1:4)
set(gca,'xlim',[-0.1,max(Time)],'ylim',[0,max(Age)])
zlabel('Age density n(a,t)')
xlabel('Time (t)')
ylabel('Age (a)')
cm = linspace(0.2,0.8,64);
colormap(repmat(cm(:),[1,3]))
