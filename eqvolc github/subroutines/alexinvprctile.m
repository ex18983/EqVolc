function [Pmin,Pmed,Pmax]=alexinvprctile(data_in,x)

% Calculates probabilities of non-exceedance for an observed value, given
% an array of simualted values. This code is designed for when the
% simulated values are 'clustered', having multiple repeated values. This
% code then outputs several probabilities, to show the range of possible
% answers

% Pmin   - Lower bound, ie probability that observed > data
% Pmed   - Avergae value, = (Pmin+Pmax)/2
% Pmax   - Upper bound, ie probability that observed >= data
% data_in- Array of simulated data
% x      - Observed value


% Determine number of simulated values
no_sim = length(data_in);

% Sort simulated values ascending
data_in = sort(data_in,1);

% Assign percentile values to the sorted data
pp = 0 : (100/no_sim) : 100;

% Calculate the indeces of the relevant values
imin = find(x<=data_in,1);
imax = find(x<data_in,1);

% Output the relevant perctile values
if isempty(imin)
    % Special cases for when observed value exceeds range of simulated values
    Pmin = 100;
else
    Pmin = pp(imin);
end
if isempty(imax)
    % Special cases for when observed value exceeds range of simulated values
    Pmax = 100;
else
    Pmax = pp(imax);
end

% Calculate average probability
Pmed = (Pmin+Pmax)/2;

end