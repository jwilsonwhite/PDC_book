function pinksalmon

% Plot R. Peterman's pink salmon dataset as an example:

% Load data
D = importdata('~/matlab_files/Botsford/Salmon/Peterman_pink.csv');

Year = D.data(:,1);
S = D.data(:,2);
R = D.data(:,3);
Labels = D.textdata(2:end,1:3);
Regions = unique(Labels(:,2));
Stocks = unique(Labels(:,1));

OK = strcmp(Labels(:,1),'Area 2E');


figure(1)
clf
set(gcf,'units','cent','position',[10 10 9,5])

plot(Year(OK),S(OK),'color','k','marker','o')
set(gca,'xlim',[1949 1996])
set(gca,'ylim',[0 1000]);
set(gca,'xtick',1950:10:2000,'ytick',0:250:1000)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
