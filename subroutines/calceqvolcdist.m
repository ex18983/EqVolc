function [eqvolcdist]=calceqvolcdist(eqlat,eqlon,volcs,novolcs)

% Calculates the distance between an earthquake and each eruption in the
% input catalogue and outputs this as an array.


% Generate array to store the distances between the earthquake and
% eruptions
eqvolcdist=zeros(1,novolcs);

% Calculate distances and extract eruptions occuring within the specified
% distance from earthquake
for i=1:novolcs
    
    % Calculate the distance between the earthquake and eruption
    [eqvolcdist(i),~]=distance(eqlat,eqlon,volcs{i,5},volcs{i,6});
    eqvolcdist(i)=deg2km(eqvolcdist(i));
    
end

end