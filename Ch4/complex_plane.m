function complex_plane

% Plot examples of complex plane & de Moivre's formula

figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 16])

% unit circle
Theta = linspace(0,2*pi,100);
[x,y] = pol2cart(Theta,1);

% The complex number - from the eigenvalues of the Ch 3 Leslie matrix

% Assemble the Leslie matrix:
F = [0 0 0 0 0.8 1 1 1 0.1]; % fecundity
S = [0.7 0.7 0.9 0.9 0.9 0.9 0.3 0.3]; % survival
Smat = [diag(S), zeros(length(S),1)]; % lower half of matrix
L = [F;Smat];  % the full matrix
Lambda = sort(eig(L),1,'descend'); % eigenvalues
z = Lambda(2);
zz = Lambda(1);

subplot(2,1,1)
hold on

% plot unit circle
plot(x,y,'k-')
plot([-2 2],[0 0],'k--')
plot([0 0],[-2 2],'k--')
plot([0 real(z)],[0,imag(z)],'k-')

% plot the points on the circles
plot(z,'ko')
axis equal
xlim([-1.1 1.1])
ylim([-1.1 1.1])
ylabel('Im','fontsize',14)
xlabel('Re','fontsize',14)
set(gca,'tickdir','out','ticklength',[0.02 0.02])

subplot(2,1,2)
hold on

% Now plot the radii to each point
for t = 1:7
    plot(z^t,'ko') % imaginary eigenvalues
    plot(zz^t,0,'kd') % real (dominant) eigenvalue
end
axis equal


% plot unit circle
plot(x,y,'k-')
plot([-2 2],[0 0],'k--')
plot([0 0],[-2 2],'k--')
plot([0 real(z)],[0,imag(z)],'k-')

plot(z,'ko')
axis equal
xlim([-1.3 1.3])
ylim([-1.3 1.3])
ylabel('Im','fontsize',14)
xlabel('Re','fontsize',14)
set(gca,'tickdir','out','ticklength',[0.02 0.02])