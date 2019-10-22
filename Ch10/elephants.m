function elephants

% Code to make Fig 10.2, showing results for Armbruster & Lande's elephant
% paper

% No need to redo model, just plotting results from Table 3

Area = [20 50 100 200 500 1000]*2.5899;
E_mean = [1, 0.976, 0.888, 0.707, 0.431, 0.272;...
          0.917, 0.485, 0.197, 0.061, 0.012, 0.003;...
          0.573, 0.098, 0.02, 0.003, 0.002, 0];
      
E_sd = [0, 0.005, 0.01, 0.014, 0.007, 0.006;...
        0.009, 0.016, 0.013, 0.003, 0.002, 0.0008;...
        0.049, 0.009, 0.004, 0.0008, 0.0006, 0];
    
 figure(1)
 clf
 set(gcf,'units','cent','position',[10 10 9 8])
 
 hold on
 for i = 1:3
     patch([Area,fliplr(Area)],[E_mean(i,:)-2*E_sd(i,:),fliplr(E_mean(i,:)+2*E_sd(i,:))],[0.5 0.5 0.5])
 end 
 
  set(gca,'tickdir','out','ticklength',[0.02 0.02],'xlim',[0 2.5e3],'fontsize',10)
  set(gca,'xcolor','k','ycolor','k')
  axis square

 xlabel('Reserve area (km^2)','fontsize',14)
 ylabel('Extinction probability','fontsize',14)
 
 