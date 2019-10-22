function otters

% Re-create Lubina & Levin's (1988) model of sea-otter range expansion. 
% Am Nat 131: 526-543

% The data: (total poplation, not split N-S)

%  Total range size
D = [11, 43, 74, 89, 108, 125, 137, 152, 158, 177,...
     192, 244, 255, 263, 279, 293, 293, 299, 312,...
     312, 353, 353];
 
% Extent of increase in N and S
dxN = [0, 11, 8, 2, 3, 11, 6, 5, 0 , 6, 0, 23,...
       6, 8, 10, 8, 0, 0, 0, 0, 26, 0];
   
dxS = [0, 21, 23, 13, 16, 6, 6, 10, 6, 13, 15, 29,...
       5, 0, 6, 6, 0, 6, 13, 0, 15, 0];

% Range size in each direction
DN = D(1) + cumsum(dxN);
DS = D(1) + cumsum(dxS);

 % Year:
Y = [1914, 1938, 1947, 1950, 1955, 1957, 1959,...
     1963, 1966, 1969, 1972:1980, 1982:1984];
 
 
 OKY1 = Y >= 1938 & Y <= 1972; % time period 1, 1938-1972
 OKY2 = Y > 1972 & Y < 1980; % time period 1, 1972 forward
 
 [bN1, ~, ~, ~, statsN1]=regress(DN(OKY1)',[ones(sum(OKY1),1),Y(OKY1)']);
 [bN2, ~, ~, ~, statsN2]=regress(DN(OKY2)',[ones(sum(OKY2),1),Y(OKY2)']);

 [bS1, ~, ~, ~, statsS1]=regress(DS(OKY1)',[ones(sum(OKY1),1),Y(OKY1)']);
 [bS2, ~, ~, ~, statsS2]=regress(DS(OKY2)',[ones(sum(OKY2),1),Y(OKY2)']);
 
 % Calculate Ds:
 r = 0.056;
 Dnorth = 13.5;
 Dsouth = 54.7;
 
 RN = 2*sqrt(r*Dnorth);
 RS = 2*sqrt(r*Dsouth);
 
 
%---------------------------------------------
% Plotting
figure(4)
clf
set(gcf,'units','cent','position',[2 30 9 9])

 hold on
 plot(Y,DS,'ko')
 plot(Y,DN,'kd')
 plot(Y(OKY1),bN1(1)+bN1(2)*Y(OKY1),'k-')
 plot(Y(OKY2),bN2(1)+bN2(2)*Y(OKY2),'k-')
 plot(Y(OKY1),bS1(1)+bS1(2)*Y(OKY1),'k-')
 plot(Y(OKY2),bS2(1)+bS2(2)*Y(OKY2),'k-')
 
 set(gca,'xlim',[1935 1985],'fontsize',12)
 set(gca,'ylim',[0 250])
 set(gca,'tickdir','out','ticklength',[0.02 0.02])
 set(gca,'xtick',1930:10:1990,'ytick',0:50:400)
 ylabel('Range size (km)','fontsize',14)
 
keyboard