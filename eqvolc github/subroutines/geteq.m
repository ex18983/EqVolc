function [eqs,noeqs]=geteq(DB_inpeqs,lonlims,latlims,start_date,end_date,minEQ,maxEQ,mindepth,maxdepth,sliptype,fashockcheck,fashockdist)

% Generates database of input earthquakes using desired parameters

% Import the raw input database earthquake data
[~,~,inpeqs]=xlsread(DB_inpeqs);

% Format the input earthquakes dates
noinpeqs=size(inpeqs,1);
for i=2:noinpeqs
    inpeqs{i,1}=datenum(inpeqs{i,1},'dd/mm/yyyy');
end

% Preallcoate output array of input earthquakes for speed (probably...)
eqs=cell(8000,11);

% Start output line counter
j=1;

% Extract events which meet filter requirements from input database
for i=2:noinpeqs
    % Event within area filter
    [inarea]=checkarea(lonlims,latlims,inpeqs{i,3},inpeqs{i,2});
    if inarea=='Y'
        % Date filter
        if inpeqs{i,1}>=start_date && inpeqs{i,1}<=end_date
            % Magnitude filter
            if inpeqs{i,5}>=minEQ && inpeqs{i,5}<maxEQ
                % Depth filter
                if inpeqs{i,4}>=mindepth && inpeqs{i,4}<maxdepth
                    % Slip type filter
                    if (strcmp(sliptype,'all') || (strcmp(sliptype,inpeqs{i,8})))
                        % Fore/aftershock filter 
                        fashockresult='N'; % Reset result
                        % Foreshock and aftershock identification routine
                        if fashockcheck>0
                            [fashockresult]=foreaftershock(fashockdist,fashockcheck,inpeqs{i,1},inpeqs{i,5},inpeqs{i,2},inpeqs{i,3},inpeqs,fashockresult);
                        end
                        if fashockresult=='N'
                            % All filters met, add event to output
                            eqs{j,1}=inpeqs{i,1}; % Event date
                            eqs{j,2}=inpeqs{i,2}; % Event latitude
                            eqs{j,3}=inpeqs{i,3}; % Event longitude
                            eqs{j,4}=inpeqs{i,4}; % Event depth
                            eqs{j,5}=inpeqs{i,5}; % Event Mw
                            eqs{j,6}=inpeqs{i,6}; % Event Location Descriptor
                            eqs{j,7}=inpeqs{i,8}; % Event slip type
                            % Move to next output line
                            j=j+1;
                        end
                    end
                end
                
            end
        end
    end
end

% Trim empty values from output cell array
eqs=eqs(~cellfun('isempty',eqs(:,1)),:);

% Calculate number of earthquakes
noeqs=size(eqs,1);

end