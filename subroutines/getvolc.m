function [volcs,novolcs]=getvolc(lonlims,latlims,first_volc,last_volc,minVEI,DB_inperups,uncert_erups)
    
% Generates database of input eruptions using desired parameters
% Set the desired input databases in the following line of code:

% Import the raw input eruption data
[~,~,inpvolcs]=xlsread(DB_inperups);

% Preallcoate output array of input eruptions for speed (probably...)
volcs=cell(2000,7);

% Start line counters
rawsize=size(inpvolcs,1); % Number of input lines (eruptions in database)
j=1; % Output line counter
% Extract events which meet filter requirements
for i=3:rawsize
    % Event within area filter
    [inarea]=checkarea(lonlims,latlims,inpvolcs{i,24},inpvolcs{i,23});
    if inarea=='Y'
        % VEI filter
        if str2double(inpvolcs{i,6})>=minVEI
            % Uncertain start year filters
            if (isempty(inpvolcs{i,8})==1 && isempty(inpvolcs{i,10})==1) || uncert_erups=='A'
                % Uncertain start day filters - no uncertainty or accept
                % all eruptions
                if (isempty(inpvolcs{i,12})==1 && isempty(inpvolcs{i,14})==1 && isempty(inpvolcs{i,13})==0 && str2double(inpvolcs{i,13})~=0 && str2double(inpvolcs{i,11})~=0) || uncert_erups=='A'
                    % If accepting all eruptions, change zero dates to 1
                    if inpvolcs{i,13}=='0'
                        inpvolcs{i,13}='1';
                    end
                    if inpvolcs{i,11}=='0'
                        inpvolcs{i,11}='1';
                    end
                    % Date boundary filter
                    syear=strcat(inpvolcs{i,13},'/',inpvolcs{i,11},'/',inpvolcs{i,9});
                    syear=datenum(syear,'dd/mm/yyyy');
                    if syear>=first_volc && syear<=last_volc
                        % Filters met, add event to output
                        volcs{j,1}=inpvolcs{i,1}; % Volcano ID
                        volcs{j,2}=inpvolcs{i,2}; % Volcano name
                        volcs{j,3}=syear; % Eruption date
                        volcs{j,4}=inpvolcs{i,6}; % VEI
                        volcs{j,5}=inpvolcs{i,23}; % Volcano latitude
                        volcs{j,6}=inpvolcs{i,24}; % Volcano longitude
                        if (isempty(inpvolcs{i,8})==1 && isempty(inpvolcs{i,10})==1 && isempty(inpvolcs{i,12})==1 && isempty(inpvolcs{i,14})==1 && isempty(inpvolcs{i,13})==0 && str2double(inpvolcs{i,13})~=0 && str2double(inpvolcs{i,11})~=0)
                            volcs{j,7}='certain'; % Eruption date qualifier
                        else
                            volcs{j,7}='uncertain'; % Eruption date qualifier
                        end
                        % Move to next output line counter
                        j=j+1;
                    end
                % Uncertain start day filters - uncertainty quantified in days
                % before or after
                elseif uncert_erups=='Y' && isempty(inpvolcs{i,14})==0 && isempty(inpvolcs{i,12})==1 && isempty(inpvolcs{i,13})==0 && str2double(inpvolcs{i,13})~=0 && str2double(inpvolcs{i,11})~=0
                    % Work out if eruption falls for certain within designated year
                    syear=strcat(inpvolcs{i,13},'/',inpvolcs{i,11},'/',inpvolcs{i,9});
                    syear=datenum(syear,'dd/mm/yyyy'); % Eruption date
                    syearmin=syear-str2double(inpvolcs{i,14});
                    syearmax=syear+str2double(inpvolcs{i,14});
                    syearmin=datestr(syearmin,'dd/mm/yyyy');
                    syearmax=datestr(syearmax,'dd/mm/yyyy');
                    certyear=strcmp(syearmin(end-3:end),syearmax(end-3:end));
                    if certyear==1 && syear>=first_volc && syear<=last_volc
                        % Filters met, add event to output
                        volcs{j,1}=inpvolcs{i,1}; % Volcano ID
                        volcs{j,2}=inpvolcs{i,2}; % Volcano name
                        syear=strcat(inpvolcs{i,13},'/',inpvolcs{i,11},'/',inpvolcs{i,9});
                        volcs{j,3}=datenum(syear,'dd/mm/yyyy'); % Eruption date
                        volcs{j,4}=inpvolcs{i,6}; % VEI
                        volcs{j,5}=inpvolcs{i,23}; % Volcano latitude
                        volcs{j,6}=inpvolcs{i,24}; % Volcano longitude
                        if (isempty(inpvolcs{i,8})==1 && isempty(inpvolcs{i,10})==1 && isempty(inpvolcs{i,12})==1 && isempty(inpvolcs{i,14})==1 && isempty(inpvolcs{i,13})==0 && str2double(inpvolcs{i,13})~=0 && str2double(inpvolcs{i,11})~=0)
                            volcs{j,7}='certain'; % Eruption date qualifier
                        else
                            volcs{j,7}='uncertain'; % Eruption date qualifier
                        end
                        % Move to next output line counter
                        j=j+1;
                    end
                % Uncertain start day filters - uncertainty quantified in days
                % before only
                elseif uncert_erups=='Y' && isempty(inpvolcs{i,14})==0 && isempty(inpvolcs{i,12})==0 && inpvolcs{i,12}=='<' && isempty(inpvolcs{i,13})==0 && str2double(inpvolcs{i,13})~=0 && str2double(inpvolcs{i,11})~=0
                    % Work out if eruption falls for certain within designated year
                    syear=strcat(inpvolcs{i,13},'/',inpvolcs{i,11},'/',inpvolcs{i,9});
                    syear=datenum(syear,'dd/mm/yyyy'); % Eruption date
                    syearmin=syear-str2double(inpvolcs{i,14});
                    syearmax=syear;
                    syearmin=datestr(syearmin,'dd/mm/yyyy');
                    syearmax=datestr(syearmax,'dd/mm/yyyy');
                    certyear=strcmp(syearmin(end-3:end),syearmax(end-3:end));
                    if certyear==1 && syear>=first_volc && syear<=last_volc
                        % Filters met, add event to output
                        volcs{j,1}=inpvolcs{i,1}; % Volcano ID
                        volcs{j,2}=inpvolcs{i,2}; % Volcano name
                        syear=strcat(inpvolcs{i,13},'/',inpvolcs{i,11},'/',inpvolcs{i,9});
                        volcs{j,3}=datenum(syear,'dd/mm/yyyy'); % Eruption date
                        volcs{j,4}=inpvolcs{i,6}; % VEI
                        volcs{j,5}=inpvolcs{i,23}; % Volcano latitude
                        volcs{j,6}=inpvolcs{i,24}; % Volcano longitude
                        if (isempty(inpvolcs{i,8})==1 && isempty(inpvolcs{i,10})==1 && isempty(inpvolcs{i,12})==1 && isempty(inpvolcs{i,14})==1 && isempty(inpvolcs{i,13})==0 && str2double(inpvolcs{i,13})~=0 && str2double(inpvolcs{i,11})~=0)
                            volcs{j,7}='certain'; % Eruption date qualifier
                        else
                            volcs{j,7}='uncertain'; % Eruption date qualifier
                        end
                        % Move to next output line counter
                        j=j+1;
                    end
                end
            end
        end
    end
end

% Trim empty values from output cell array
volcs=volcs(~cellfun('isempty',volcs(:,1)),:);

% Calculate number of eruptions in database
novolcs=size(volcs,1);

end