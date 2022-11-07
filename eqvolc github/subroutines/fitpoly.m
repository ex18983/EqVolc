function [eqspolyvals,eqerror,eqval,volcspolyvals,volcerror,volcval,momentpolyvals,momenterror,momentval]=fitpoly(curvemode,n_min,n_max,lspan,xyears,eqs,volcs,seismoment,autochoose)

% Fits the desired curve to the input earthquake and eruption rate
% data. 
% This code uses MATLAB polyfit function and an internet-sourced equation 
% to select best fitting polynomial
% This code uses the downloaded LOESS fitting software from 
% https://uk.mathworks.com/matlabcentral/fileexchange/22470-lowess-locally-weighted-scatterplot-smoothing-for-linear-and-non-linear-data-enhanced

% Generate vector to store computed variance values
eqsvars=zeros(1,n_max-n_min+1);
volcsvars=zeros(1,n_max-n_min+1);
momentvars=zeros(1,n_max-n_min+1);

if strcmp(curvemode,'poly')==1
    
    % Generate vector of polynomial orders to try
    nints=n_min:1:n_max;

    % Loop to find the bestfit polynomial for observed data
    for i=n_min:n_max
        
        % Fit polynomial of order i to the input earthquake and eruption data
        % (using weird matlab error and scaling functions to improve it)
        [eqspoly,eqerror,eqmu]=polyfit(xyears,eqs,i);
        [volcspoly,volcerror,volcmu]=polyfit(xyears,volcs,i);
        [momentpoly,momenterror,momentmu]=polyfit(xyears,seismoment,i);
        
        % Generate the polynomial values at the input years (using weird matlab
        % error and scaling functions to improve it)
        [eqspolyvals,~]=polyval(eqspoly,xyears,eqerror,eqmu);
        [volcspolyvals,~]=polyval(volcspoly,xyears,volcerror,volcmu);
        [momentpolyvals,~]=polyval(momentpoly,xyears,momenterror,momentmu);
        
        % Calculate residuals
        eqsres=eqs-eqspolyvals;
        volcsres=volcs-volcspolyvals;
        momentres=seismoment-momentpolyvals;
        % Square residuals
        eqsres=eqsres.*eqsres;
        volcsres=volcsres.*volcsres;
        momentres=momentres.*momentres;
        % Sum squares of residuals
        eqsressum=sum(eqsres);
        volcsressum=sum(volcsres);
        momentressum=sum(momentres);
        
        % Calculate the variance foruma (from
        % https://autarkaw.org/2008/07/05/finding-the-optimum-polynomial-order-to-use-for-regression/)
        eqsvars(i)=eqsressum/(size(eqs,2)-i-1);
        volcsvars(i)=volcsressum/(size(eqs,2)-i-1);
        momentvars(i)=momentressum/(size(eqs,2)-i-1);
    end
    
    % Choose the best polynomial for observed data
    if autochoose=='N'
        disp('Examine the earthquakes polynomial variance (eqsvars) and choose the bestfit polynomial order. Continue code when decided then enter value.')
        eqval=input('Enter desired earthquake polynomial order:     ');
        disp('Examine the eruptions polynomial variance (volcsvars) and choose and enter the bestfit polynomial order. Continue code when decided then enter value.');
        volcval=input('Enter the desired eruption polynomial order:     ');
        disp('Examine the eruptions polynomial variance (momentvars) and choose and enter the bestfit polynomial order. Continue code when decided then enter value.');
        momentval=input('Enter the desired eruption polynomial order:     ');
        
    elseif autochoose=='Y'
        [~,eqval]=min(eqsvars);
        [~,volcval]=min(volcsvars);
        [~,momentval]=min(momentvars);
        eqval=nints(eqval);
        volcval=nints(volcval);
        momentval=nints(momentval);
    end
    
    % Recalculate the chosen bestfit polynomials for observed data
    [eqspoly,eqerror,eqmu]=polyfit(xyears,eqs,eqval);
    [volcspoly,volcerror,volcmu]=polyfit(xyears,volcs,volcval);
    [momentpoly,momenterror,momentmu]=polyfit(xyears,seismoment,momentval);
    
    % Generate the polynomial values at the input years (using weird matlab
    % error and scaling functions to improve it)
    [eqspolyvals,eqerror]=polyval(eqspoly,xyears,eqerror,eqmu);
    [volcspolyvals,volcerror]=polyval(volcspoly,xyears,volcerror,volcmu);
    [momentpolyvals,momenterror]=polyval(momentpoly,xyears,momenterror,momentmu);
    
elseif strcmp(curvemode,'LOESS')==1
    
    % Call the LOWESS function of Jeff Burkey to fit LOWESS curves
    [eqsout,~,~,~]=lowess([xyears' eqs'],lspan,0);
    [volcsout,~,~,~]=lowess([xyears' volcs'],lspan,0);
    [momentout,~,~,~]=lowess([xyears' seismoment'],lspan,0);
    
    % Assign the LOWESS values to output
    eqspolyvals=[eqsout(:,3)]';
    volcspolyvals=[volcsout(:,3)]';
    momentpolyvals=[momentout(:,3)]';
    % Assign null values for other outputs not caluculated for LOWESS
    eqerror=0;
    eqval=0;
    volcerror=0;
    volcval=0;
    momenterror=0;
    momentval=0;
    
elseif strcmp(curvemode,'movingmean')
    
    % Calculate the moving averages
    eqspolyvals=movmean(eqs,lspan);
    volcspolyvals=movmean(volcs,lspan);
    momentpolyvals=movmean(seismoment,lspan);
    
    % Assign null values for other outputs not caluculated for moving mean
    eqerror=0;
    eqval=0;
    volcerror=0;
    volcval=0;
    momenterror=0;
    momentval=0;
    
end

end