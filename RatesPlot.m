function RatesPlot

% Plots earthquake and eruption time-series generated in RatesCalc.m and
% fits a curve for visualisation

% Also calculates the Spearman cross-correllation between the earthquake
% and volcano time-series, for observed data and any monte carlo
% simulations (for the raw data, not the curves)

% Figure 1 - Earthquake rates
% Figure 2 - Eruption rates
% Figure 3 - Seismic moment release
% Figure 4 - Spearman cross-correlation between earthquake rates and
%            eruption rates time-series
% Figure 5 - Spearman cross-correlation between seismic moment release and
%            eruption rates time-series
% Figure 6 - Combination of Figures 1,2,3, and 5

% (Figure legends and axes can be manually formatted because this is
% complicated to auto format. Scroll down to find the plotting commands to
% do this and alter appropriately depending on what your time-series are
% (line 260 onwards))


% INPUT

input_file='rates_test.mat'; % Input .mat file with earthquake and eruption rates from RatesCalc.m

resamp=1; % Bin size for time-series; choose 1 to keep original calendar year bins, otherwise choose a factor of the time-series length
curvemode='LOESS'; % Curve fitting method ('poly' or 'LOESS' or 'movingmean' or 'none')
n_min=8; % If 'poly' curvemode, Minimum polynomial fit to try
n_max=8; % If 'poly' curvemode, Maximum polynomial fit to try
lspan=0.15; % If 'LOESS' curvemode, Span parameter of loess fitting (0 - 1). If 'movingmean' curvemode, window size over which to calculate sliding moving mean
autochoose='Y'; % Autoselect the best fit polynomial ('Y' = autoselect, 'N' = option to look at the data and manually choose)
tlag=15; % Maximum time shift for cross-correlation (years)
inp_percent=[1 5 32 50 68 95 99]; % Percentile values to plot for any Monte Carlo simulations
percent_color={[1 0 0] [1 0 0] [1 0 0] [1 0 0] [1 0 0] [1 0 0] [1 0 0] [1 0 0] [1 0 0]}; % Plotting colours for percentile data (as many specs as desired percentiles in inp_percent)
percent_style={':' '-' '--' '-' '--' '-' ':'}; % Plotting line styles for percentile data (as many specs as desired percentiles in inp_percent)
plot_moment='Y'; % Generate plots of seismic moment against eruption rate as well? ('Y' or 'N')

% Define the home directory for function
home_dir= 'C:\Users\ex18983\eqvolc github';


% CALCULATION

% Switch to home directory
cd(home_dir)

% Load the input file
load(input_file);

% Determine the number of yearly bins
no_years=size(eqs,2);

% Determine the number of monte carlo simulations
if exist('eqsMC','var')
    no_MC=size(eqsMC,1);
else
    % dummy values if no simulations
    no_MC=0;
    eqsMC=zeros(1,no_years);
    volcsMC=zeros(1,no_years);
    seismomentMC=zeros(1,no_years);
end

% Apply resampling if applicable
if resamp~=1
    xyears=xyears(1:resamp:end);
    eqsresamp=zeros(1,no_years/resamp);
    volcsresamp=zeros(1,no_years/resamp);
    seismomentresamp=zeros(1,no_years/resamp);
    eqsMCresamp=zeros(no_MC,no_years/resamp);
    volcsMCresamp=zeros(no_MC,no_years/resamp);
    seismomentMCresamp=zeros(no_MC,no_years/resamp);
    for i=1:no_years/resamp
        eqsresamp(i)=(sum(eqs((i-1)*resamp+1:(i-1)*resamp+resamp)))/resamp;
        volcsresamp(i)=(sum(volcs((i-1)*resamp+1:(i-1)*resamp+resamp)))/resamp;
        seismomentresamp(i)=(sum(seismoment((i-1)*resamp+1:(i-1)*resamp+resamp)))/resamp;
    end
    for i=1:no_MC
        for j=1:no_years/resamp
            eqsMCresamp(i,j)=(sum(eqsMC(i,(j-1)*resamp+1:(j-1)*resamp+resamp)))/resamp;
            volcsMCresamp(i,j)=(sum(volcsMC(i,(j-1)*resamp+1:(j-1)*resamp+resamp)))/resamp;
            seismomentMCresamp(i,j)=(sum(seismomentMC(i,(j-1)*resamp+1:(j-1)*resamp+resamp)))/resamp;
        end
    end
    eqs=eqsresamp;
    volcs=volcsresamp;
    seismoment=seismomentresamp;
    eqsMC=eqsMCresamp;
    volcsMC=volcsMCresamp;
    seismomentMC=seismomentMCresamp;
end

% (Re)determine the number of yearly bins
no_years=size(eqs,2);

% Convert sesimic moment to log values - might have to do this then convert
% back after LOESS, as LOESS can mess up with big exponential differences-
% bit of a fudge but works better than not doing it! SEE ALSO LINE 139
% moment=log10(seismoment);
% momentMC=log10(seismomentMC);
% Don't convert seismic moment to log values
moment=seismoment;
momentMC=seismomentMC;
% Fix to remove -Inf values from log of seismic moments (arising from years
% with 0 eruptions in 1976 in Aleutian - using the Mw >5.5 gives 10^18)
moment(moment==-Inf)=18.1;
momentMC(momentMC==-Inf)=18.1;
moment(moment==0)=1.2e18;
momentMC(momentMC==0)=1.2e18;

% Generate vector of timelages
tlagslabel=-tlag:resamp:tlag;
no_lags=size(tlagslabel,2);
tlags=-tlag/resamp:1:tlag/resamp;

% Fit the curves for the observed data
if strcmp(curvemode,'none')==0
    [eqspv,~,eqval,volcspv,~,volcval,momentpv,~,momentval]=fitpoly(curvemode,n_min,n_max,lspan,xyears,eqs,volcs,moment,autochoose);
end

% Generate array to store the curve fits for the monte carlo
% simulations
if no_MC>0
    eqspvMC=zeros(no_MC,no_years);
    volcspvMC=zeros(no_MC,no_years);
    momentpvMC=zeros(no_MC,no_years);
    
    % Generate array to store the correlation coefficients for the monte
    % carlo simulations
    Rmcraw=zeros(no_MC,no_lags);
    Rmcrawm=zeros(no_MC,no_lags);
    
    % Generate arrays to store the mean_SD and percentiles data for the
    % monte carlo simulations
    meanSDraw=zeros(2,no_lags);
    meanSDrawm=zeros(2,no_lags);
    pctsraw=zeros(size(inp_percent,2),no_lags);
    pctsrawm=zeros(size(inp_percent,2),no_lags);
    
    if strcmp(curvemode,'none')==0
        for i=1:no_MC
            % For each monte carlo simulation result, calculate the curve fit
            [eqspvMC(i,:),~,~,volcspvMC(i,:),~,~,momentpvMC(i,:),~,~]=fitpoly(curvemode,n_min,n_max,lspan,xyears,eqsMC(i,:),volcsMC(i,:),momentMC(i,:),autochoose);
        end
    end
end

% If previously converted moment to log values to prevent the LOESS messing
% up, switch them back now before analysis
% moment=10.^moment;
% momentpv=10.^momentpv;
% momentMC=10.^momentMC;
% momentpvMC=10.^momentpvMC;

% Create copies of timeseries for plotting purposes
eqs1=eqs;
volcs1=volcs;
moment1=moment;
if strcmp(curvemode,'none')==0
    eqspv1=eqspv;
    volcspv1=volcspv;
    momentpv1=momentpv;
end

% Generate array to store the observed correlation coefficients (arrays
% to store Monte Carlo previously generated)
Rraw=zeros(1,no_lags);
Rrawm=zeros(1,no_lags);

% Loop to calculate the correlation coefficients at each timestep
for i=1:no_lags
    % Negative timeshift
    if tlags(i)<0
        Rraw(i)=corr(eqs(1:end+tlags(i))',volcs(1-tlags(i):end)','type','Spearman');
        Rrawm(i)=corr(moment(1:end+tlags(i))',volcs(1-tlags(i):end)','type','Spearman');
        % Loop for Monte Carlo simulations
        for j=1:no_MC
            Rmcraw(j,i)=corr(eqsMC(j,1:end+tlags(i))',volcsMC(j,1-tlags(i):end)','type','Spearman');
            Rmcrawm(j,i)=corr(momentMC(j,1:end+tlags(i))',volcsMC(j,1-tlags(i):end)','type','Spearman');
        end
    % Positive timeshift
    else
        Rraw(i)=corr(eqs(1+tlags(i):end)',volcs(1:end-tlags(i))','type','Spearman');
        Rrawm(i)=corr(moment(1+tlags(i):end)',volcs(1:end-tlags(i))','type','Spearman');
        % Loop for Monte Carlo simulations
        for j=1:no_MC
            Rmcraw(j,i)=corr(eqsMC(j,1+tlags(i):end)',volcsMC(j,1:end-tlags(i))','type','Spearman');
            Rmcrawm(j,i)=corr(momentMC(j,1+tlags(i):end)',volcsMC(j,1:end-tlags(i))','type','Spearman');
        end
    end
end

% Find the highest (+ve) correlations of observed data 
[Rmaxobsraw,Rposobsraw]=max(Rraw);
[Rmaxobsrawm,Rposobsrawm]=max(Rrawm);
% Calculate time lags of highest correlation
tlagobsraw=tlags(Rposobsraw);
tlagobsrawm=tlags(Rposobsrawm);

% Run a loop for each timelag to determine the mean, SD and percentiles
% of the monte carlo data
if no_MC>0
    for i=1:size(tlags,2)
        meanSDraw(1,i)=mean(Rmcraw(:,i));
        meanSDraw(2,i)=std(Rmcraw(:,i),1);
        meanSDrawm(1,i)=mean(Rmcrawm(:,i));
        meanSDrawm(2,i)=std(Rmcrawm(:,i),1);
        % Run a loop for each desired percentile
        for j=1:size(inp_percent,2)
            pctsraw(j,i)=prctile(Rmcraw(:,i),inp_percent(j));
            pctsrawm(j,i)=prctile(Rmcrawm(:,i),inp_percent(j));
        end
    end
    
    % Calculate the percentile of the highest observed correlation relative to
    % monte carlo data at that timelag (big P values)
    Pobsraw=num2str((sum(Rmaxobsraw>=Rmcraw(:,Rposobsraw)))/no_MC*100);
    Pobsrawm=num2str((sum(Rmaxobsrawm>=Rmcrawm(:,Rposobsrawm)))/no_MC*100);
    % Change if P = 100
    if strcmp(Pobsraw,'100')
        Pobsraw='\geq99.99';
    end
    if strcmp(Pobsrawm,'100')
        Pobsrawm='\geq99.99';
    end
    
    % Calculate actual p-vales using the Davison 2003 formula
    pobsraw=num2str((sum(Rmcraw(:,Rposobsraw)>=Rmaxobsraw)+1)/(1+no_MC));
    pobsrawm=num2str((sum(Rmcrawm(:,Rposobsrawm)>=Rmaxobsrawm)+1)/(1+no_MC));
    % Change if no simulated values stronger than observed correlation
    pmodraw='';
    pmodrawm='';
    if sum(Rmcraw(:,Rposobsraw)>=Rmaxobsraw)==0
        pmodraw='\leq';
    end
    if sum(Rmcrawm(:,Rposobsrawm)>=Rmaxobsrawm)==0
        pmodrawm='\leq';
    end
    
    % Create p-values for simulated timeseries to see how many signfiicant
    % results there are for each
    pMCrawm=zeros(no_MC,no_lags);
    simsig=zeros(1,no_lags); % index 1 = no significant values, index 2 = significant value etc...
    for i=1:no_MC
        for j=1:no_lags
            pMCrawm(i,j)=(sum(Rmcrawm(:,j)>=Rmcrawm(i,j))+1)/(1+no_MC);
        end
        nosig=sum(pMCrawm(i,:)<=0.05);
        simsig(nosig+1)=simsig(nosig+1)+1;
    end
    % simsig=simsig/no_MC;
    
    % Create labels for plots of max correlation details
    Lobsraw={strcat('\leftarrow R=',num2str(Rmaxobsraw)),strcat('     P=',num2str(Pobsraw)),strcat('     p=',pmodraw,num2str(pobsraw))};
    Lobsrawm={strcat('\leftarrow R=',num2str(Rmaxobsrawm)),strcat('     P=',num2str(Pobsrawm)),strcat('     p=',pmodrawm,num2str(pobsrawm))};
end

% Plot observed rates and their best fitting polynomials
figure(1)
hold on
plot(xyears,eqs1,'-bo','LineWidth',0.5)
if strcmp(curvemode,'none')==0
    plot(xyears,eqspv1,'--b','LineWidth',1)
end
%plot(xyears,eqspolyvals+2*eqerror,'--k','LineWidth',0.1)
%plot(xyears,eqspolyvals-2*eqerror,'--k','LineWidth',0.1)
% Format plot
title('Earthquake rates')
xlabel('Year')
ylabel('Frequency')
if strcmp(curvemode,'poly')
    legend('Earthquakes Mw\geq7',strcat('Polyfit(',num2str(eqval),') eqs'),'Polyfit 95%','Polyfit 95%')
elseif strcmp(curvemode,'LOESS')
    legend('Earthquakes Mw\geq7',strcat('LOESS(',num2str(lspan),') eqs'),'Polyfit 95%','Polyfit 95%')
elseif strcmp(curvemode,'movingmean')
    legend('Earthquakes Mw\geq7',strcat('Moving mean(',num2str(lspan),') eqs'),'Polyfit 95%','Polyfit 95%')
elseif strcmp(curvemode,'none')
    legend('Earthquakes Mw\geq7')
end
    
% Plot cross correlation results (raw)
figure(4)
hold on
% Plot every monte carlo raw result
%if no_MC>0
%    for i=1:no_MC
%        plot(tlags,mcraw(i,:),'k','LineWidth',0.1)
%    end
%end
%Plot the mean, SD, and percentiles data for raw version
%plot(tlags,meanSDraw(1,:)+2*meanSDraw(2,:),'-.k','LineWidth',0.5) %mean+2sd
%plot(tlags,meanSDraw(1,:)+meanSDraw(2,:),'--k','LineWidth',0.5) %mean+1sd
%plot(tlags,meanSDraw(1,:),'-k','LineWidth',0.5) %mean
%plot(tlags,meanSDraw(1,:)-meanSDraw(2,:),'--k','LineWidth',0.5,'HandleVisibility','off') %mean-1sd
%plot(tlags,meanSDraw(1,:)-2*meanSDraw(2,:),'-.k','LineWidth',0.5,'HandleVisibility','off') %mean+2sd
if no_MC>0
    for i=1:size(inp_percent,2)
        if i==ceil(size(inp_percent,2)/2)
            plot(tlagslabel,pctsraw(i,:),'color',percent_color{i},'LineStyle',percent_style{i},'LineWidth',0.5) % percentiles loop
            text(tlagslabel(2),pctsraw(i,2),num2str(inp_percent(i)),'FontSize',8,'color','r')
            text(tlagslabel(end-1),pctsraw(i,end-1),num2str(inp_percent(i)),'FontSize',8,'color','r')
        else
            plot(tlagslabel,pctsraw(i,:),'color',percent_color{i},'LineStyle',percent_style{i},'LineWidth',0.5,'HandleVisibility','off')
            text(tlagslabel(2),pctsraw(i,2),num2str(inp_percent(i)),'FontSize',8,'color','r')
            text(tlagslabel(end-1),pctsraw(i,end-1),num2str(inp_percent(i)),'FontSize',8,'color','r')
        end
    end
end
% Plot the actual observed coefficients over the top of the Monte Carlo
plot(tlagslabel,Rraw,'-kx','LineWidth',1) 
% Add label
if no_MC>0
    text(tlagobsraw,Rmaxobsraw,Lobsraw)
end
% Format plot
ylim([-1 1])
xlabel('Time lag (years)')
title('Earthquake-eruption xcorr (raw)')
ylabel('Correlation')
if no_MC>0
    legend('Simulated Percentiles','Observed Data','Location','best')
else 
    legend('Observed Data','Location','best')
end

% Seismic moment plots
if plot_moment=='Y'   
    
    % Plot the eruptions by themselves (for seismic moment version)
    figure(2)
    hold on
    plot(xyears,volcs1,'-ro','LineWidth',0.5)
    if strcmp(curvemode,'none')==0
        plot(xyears,volcspv1,'--r','LineWidth',1)
    end
    % Format plot
    title('Eruption rates')
    xlabel('Year')
    ylabel('Frequency')
    if strcmp(curvemode,'poly')
        legend('Eruptions VEI\geq2',strcat('Polyfit(',num2str(volcval),') erups'))
    elseif strcmp(curvemode,'LOESS')
        legend('Eruptions VEI\geq2',strcat('LOESS(',num2str(lspan),') erups'))
    elseif strcmp(curvemode,'movingmean')
        legend('Eruptions VEI\geq2',strcat('Moving mean(',num2str(lspan),') erups'))
    elseif strcmp(curvemode,'none')
        legend('Eruptions VEI\geq2')
    end
    
    % Plot the log of moment by itself (for seismic moment version)
    figure(3)
    semilogy(xyears,moment1,'-bo','LineWidth',0.5)
    hold on
    if strcmp(curvemode,'none')==0
        plot(xyears,momentpv1,'--b','LineWidth',1)
    end
    % Format plot
    title('Cumulative yearly log(seismic moment)')
    xlabel('Year')
    ylabel('log(Cumulative yearly Mo)')
    if strcmp(curvemode,'poly')
        legend('log(Mo)',strcat('Polyfit(',num2str(momentval),') log(Mo)'))
    elseif strcmp(curvemode,'LOESS')
        legend('log(Mo)',strcat('LOESS(',num2str(lspan),') log(Mo)'))
    elseif strcmp(curvemode,'movingmean')
        legend('log(Mo)',strcat('Moving mean(',num2str(lspan),') log(Mo)'))
    elseif strcmp(curvemode,'none')
        legend('log(Mo)')
    end
    
    % Plot cross correlation results (raw) (for seismic moment version)
    figure(5)
    hold on
    % Plot every monte carlo raw result
    % for i=1:no_MC
    %     plot(tlags,mcrawmoment(i,:),'k','LineWidth',0.1)
    % end
    % Plot the mean, SD, and percentiles data for raw version
%     plot(tlags,meanSDrawm(1,:)+2*meanSDrawm(2,:),'-.k','LineWidth',0.5) %mean+2sd
%     plot(tlags,meanSDrawm(1,:)+meanSDrawm(2,:),'--k','LineWidth',0.5) %mean+1sd
%     plot(tlags,meanSDrawm(1,:),'-k','LineWidth',0.5) %mean
%     plot(tlags,meanSDrawm(1,:)-meanSDrawm(2,:),'--k','LineWidth',0.5,'HandleVisibility','off') %mean-1sd
%     plot(tlags,meanSDrawm(1,:)-2*meanSDrawm(2,:),'-.k','LineWidth',0.5,'HandleVisibility','off') %mean+2sd
    if no_MC>0
        for i=1:size(inp_percent,2)
            if i==ceil(size(inp_percent,2)/2)
            plot(tlagslabel,pctsrawm(i,:),'color',percent_color{i},'LineStyle',percent_style{i},'LineWidth',0.5) % percentiles loop
            text(tlagslabel(2),pctsrawm(i,2),num2str(inp_percent(i)),'FontSize',8,'color','r')
            text(tlagslabel(end-1),pctsrawm(i,end-1),num2str(inp_percent(i)),'FontSize',8,'color','r')
        else
            plot(tlagslabel,pctsrawm(i,:),'color',percent_color{i},'LineStyle',percent_style{i},'LineWidth',0.5,'HandleVisibility','off')
            text(tlagslabel(2),pctsrawm(i,2),num2str(inp_percent(i)),'FontSize',8,'color','r')
            text(tlagslabel(end-1),pctsrawm(i,end-1),num2str(inp_percent(i)),'FontSize',8,'color','r')
            end
        end
    end
    % Plot the actual observed noneicients over the top of the Monte Carlo
    plot(tlagslabel,Rrawm,'-kx','LineWidth',1)
    % Add label
    if no_MC>0
        text(tlagobsrawm,Rmaxobsrawm,Lobsrawm)
    end
    % Format plot
    ylim([-1 1])
    xlabel('Time lag (years)')
    title('Seismic moment-eruption xcorr (raw)')
    ylabel('Correlation')
    if no_MC>0
        legend('Simulated Percentiles','Observed Data','Location','best')
    else
        legend('Observed Data','Location','best')
    end
    
    % Combined subplot figure showing frequency, moment, and erups and corr
    figure(6)
    subplot('Position',[0.05 0.68 0.35 0.3])
    plot(xyears,eqs1,'bo','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')
    hold on
    if strcmp(curvemode,'none')==0
        plot(xyears,eqspv1,'-b','LineWidth',1)
    end
    a=gca;
    a.FontSize=12;
    a.XTick=ceil(xyears(1)/10)*10 : 10 : floor(xyears(end)/10)*10;
    a.XTickLabel={};
    a.LineWidth=0.75;
    %title('a) M_w \geq 7 Earthquake frequency', 'Units', 'normalized', 'Position', [0.5, 0.85, 0],'FontSize',14);
    ylabel({'M_w \geq 7 Earthquakes/', 'year'},'FontSize',14)
    xlim([xyears(1) xyears(end)])
    subplot('Position',[0.05 0.38 0.35 0.3])
    semilogy(xyears,moment1,'bo','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')
    hold on
    if strcmp(curvemode,'none')==0
        semilogy(xyears,momentpv1,'-b','LineWidth',1)
    end
    a=gca;
    a.FontSize=12;
    a.XTick=ceil(xyears(1)/10)*10 : 10 : floor(xyears(end)/10)*10;
    a.XTickLabel={};
    a.LineWidth=0.75;
    %title('b) M_w \geq 7 Seismic moment release', 'Units', 'normalized', 'Position', [0.5, 0.9, 0],'FontSize',14);
    ylabel('M_o/ year (Nm)','FontSize',14)
    xlim([xyears(1) xyears(end)])
    subplot('Position',[0.05 0.08 0.35 0.3])
    plot(xyears,volcs1,'ro','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','r')
    hold on
    if strcmp(curvemode,'none')==0
        plot(xyears,volcspv1,'-r','LineWidth',1)
    end
    a=gca;
    a.FontSize=12;
    a.XTick=ceil(xyears(1)/10)*10 : 10 : floor(xyears(end)/10)*10;
    a.LineWidth=0.75;
    %title('c) VEI \geq 2 Eruption frequency', 'Units', 'normalized', 'Position', [0.5, 0.9, 0],'FontSize',14);
    xlabel('Year','FontSize',14)
    ylabel({'VEI \geq 2 Eruptions/', 'year'},'FontSize',14)
    xlim([xyears(1) xyears(end)])
    subplot('Position',[0.45 0.1 0.5 0.8]) % 3 plot for global
%     subplot('Position',[0.45 0.1 0.36 0.55]) % 2 plot for regional
    if no_MC>0
        for i=1:size(inp_percent,2)
            if i==ceil(size(inp_percent,2)/2)
                plot(tlagslabel,pctsrawm(i,:),'color',percent_color{i},'LineStyle',percent_style{i},'LineWidth',1) % percentiles loop
                hold on
                text(tlagslabel(2),pctsrawm(i,2),num2str(inp_percent(i)),'FontSize',12,'color','r')
                text(tlagslabel(end-1),pctsrawm(i,end-1),num2str(inp_percent(i)),'FontSize',12,'color','r')
            else
                plot(tlagslabel,pctsrawm(i,:),'color',percent_color{i},'LineStyle',percent_style{i},'LineWidth',1,'HandleVisibility','off')
                hold on
                text(tlagslabel(2),pctsrawm(i,2),num2str(inp_percent(i)),'FontSize',12,'color','r')
                text(tlagslabel(end-1),pctsrawm(i,end-1),num2str(inp_percent(i)),'FontSize',12,'color','r')
            end
        end
    end
    plot(tlagslabel,Rrawm,'ko','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','k')
    %text(tlagobsrawm,Rmaxobsrawm,Lobsrawm)
    ylim([-1 1])
    xlabel('Timeshift (years)','FontSize',14)
    %title('d) Seismic moment - Eruptions cross-correlation', 'Units', 'normalized', 'Position', [0.5, 0.85, 0],'FontSize',14);
    ylabel('Correlation (M_o vs Eruptions)','FontSize',14)
    if no_MC>0
        a=legend('PERM','Observed','Location','best');
    else
        a=legend('Observed','Location','best');
    end
    a.FontSize=14;
    a=gca;
    a.LineWidth=0.75;
    a.FontSize=12;
    
end

end