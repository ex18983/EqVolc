function [fashockresult]=foreaftershock(fashockdist,fashockcheck,EQdate,EQMw,EQlat,EQlon,inpeqs,fashockresult)

% Determines if the input earthquake is a foreshock or aftershock of a
% larger event within the time period defined. Fore/aftershocks are
% excluded from the input database

% Loop through the complete earthquake database to check if the event in
% question is a fore/aftershock
for i=2:size(inpeqs,1)
    % Earthquake falls within 'maxt' days from event in question
    % EDIT - comment out -fashockcheck to run aftershock only, not foreshock
    if inpeqs{i,1}<=EQdate+fashockcheck && inpeqs{i,1}>=EQdate-fashockcheck
        % Calculate the distance between the earthquake and event in question
        [distb,~]=distance(EQlat,EQlon,inpeqs{i,2},inpeqs{i,3});
        distb=deg2km(distb);
        % Earthquake falls within 'maxd' of event in question
        if distb<=fashockdist
            % Earthquake is a greater magnitude than event in question
            if inpeqs{i,5}>EQMw
                % Foreshock or aftershock identified, do not include
                % this event
                fashockresult='Y';
            end      
        end
    end
end
end