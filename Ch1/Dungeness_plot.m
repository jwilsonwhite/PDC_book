function Dungeness_plot

% Botsford et al. Fig 1.1
% Plot of California Dungeness crab catch over time

% Load in data.
% Downloaded from NOAA's Fisheries Information System Program
% http://www.st.nmfs.noaa.gov/data/fis/
D = importdata('MF_ANNUAL_LANDINGS.RESULTS_edited.csv');
D = D.data(:,1:2); % Col 1 = Year; Col 2 = Annual landings (metric tons)
D(:,2) = D(:,2)./1e3; % convert to millions of kg

figure(1)
clf
set(gcf,'position',[500 500 500 300])

h = plot(D(:,1),D(:,2));
set(h,'color','k','linewidth',1.5)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xlim',[1949 2014],'xtick',1950:5:2015)
set(gca,'xticklabel',{'1950','','1960','','1970','','1980','','1990','','2000',''})
set(gca,'fontsize',12)
xlabel('Year','fontsize',16)
ylabel('Dungeness crab catch (10^6 kg)','fontsize',16)

print -depsc2 'Fig1.1_Dungeness_timeseries.eps'

