function extinction_timeseries

% Plot data from Dennis et al. 1991
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 25])

Datasets = {'whoopingcrane','parrot','warbler','grizzly','condor',...
            'palila','laysan'};

     
        
for i = 1:7
    subplot(7,1,i)

D = importdata(strcat('data/',Datasets{i},'.xlsx'));
T = D.data(:,1);
N = D.data(:,2);
plot(T,N,'ko','linestyle','-')

set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k')

switch Datasets{i}
    case 'whoopingcrane'
        xlim([1935 1995]); ylim([0 200])
        ylabel('Total population size','fontsize',10)
    case 'parrot'
        xlim([1968 1990]); ylim([0 40])
        ylabel('Maximum count, Jan-Apr','fontsize',10)
    case 'warbler'
        xlim([1950 1990]); ylim([0 450])
        ylabel('Singing male count','fontsize',10)
    case 'grizzly'
        xlim([1958 1987]); ylim([0 40])
        ylabel('Estimated no. adult females','fontsize',10)
    case 'condor'
        xlim([1963 1983]); ylim([0 55])
        ylabel('Established population size','fontsize',10)
    case 'palila'
        xlim([1973 1989]); ylim([0 100])
        ylabel('Estimated population size (x10^3)','fontsize',10)
    case 'laysan'
        xlim([1969 1993]); ylim([0 25e3])
        ylabel('Estimated population size (x10^3)','fontsize',10)
    
end


end

subplot(7,1,7)
xlabel('Year','fontsize',12)


