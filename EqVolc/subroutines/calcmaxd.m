function [maxd] = calcmaxd(maxd,Mw)

% Calculates a variable maximum search distance between earthquakes and
% eruptions, as a function of the input earthquake magnitude, using
% relations based on rupute length and the rate of decay of stress change
% with distance


% Extract the scaling number
loc = strfind(maxd,' ');
x = str2double(maxd(1:loc-1));

% Determine scaling method
if strcmp(maxd(loc+1:end),'WC')

    % Use Wells & Coppersmith 1994 to calculate rupture length
    dist = 10^((Mw-5.08)/1.16);
    
    % Find final maxd
    maxd = x * dist;
    maxd = round(maxd,1);
    
elseif strcmp(maxd(loc+1:end),'Mo')
    
    % Determine moment of earthquake (in Nm, not dyne cm)
    Mo = 10^(1.5*Mw+16.1);
    Mo = Mo * 10^(-7);
    
    % Scale by value * Mo^(1/3) (Note: a value of 4e-6 km gives a distance
    % of around 5 km for a Mw 6, and 150 km for a Mw 9.1 (Churchill et al.
    % in prep)
    maxd = x * Mo^(1/3);

end