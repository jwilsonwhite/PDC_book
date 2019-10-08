function surv_fec_curves

% Make function showing survival, fecundity, and their product as a
% function of age

Ages = 0:50;
Surv = 0.2;
MatM = 10;
MatSD = 2;
Linf = 10;
k = 0.1;
a = 1;
b = 3;

S = exp(-Surv.*Ages);
L = Linf.*(1-exp(-k.*Ages));
Mat = normcdf(Ages,MatM,MatSD);
E = a.*(L.*Mat).^b;
E = E./max(E);

figure(1)
clf
hold on
set(gcf,'units','cent','position',[20 20 9 8])
plot(Ages,S,'k-')
plot(Ages,E./max(E),'k:')
plot(Ages,S.*E*50,'k--')
set(gca,'tickdir','out','ticklength',[0.015 0.015],'fontsize',8)
set(gca,'xtick',0:10:100)
set(gca,'ytick',0:0.4:2)
set(gca,'fontsize',10)
ylabel('Age density, n(a,t)','fontsize',14)
xlabel('Age, a','fontsize',14)


