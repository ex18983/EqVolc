function [eqs,volcs,indvolcs,no_indvolcs,alldates,eqvolcdist] = randomcat(eqs,noeqs,volcs,novolcs,first_eq,last_eq,first_volc,last_volc,MC_mode,erup_lim,eqvolcdist,indvolcs,no_indvolcs,alldates)

% Generates simulated earthquake and eruption catalogues using a variety of
% methods

% Randomise earthquake oreruption times
if MC_mode=='V' || MC_mode=='B'
    for j=1:novolcs
        volcs{j,3}=round(last_volc-rand*(last_volc-first_volc));
    end
end

% Randomise earthquake times
if MC_mode=='E' || MC_mode=='B'
    for j=1:noeqs
        eqs{j,1}=round(last_eq-rand*(last_eq-first_eq));
    end
end

% Permutate earthquake and eruption times
if MC_mode=='P'
    
    % Special case for first simulation - extract observed eruption dates
    if isempty(alldates)
        alldates=zeros(noeqs+novolcs,1);
        for j=1:noeqs
            alldates(j)=eqs{j,1};
        end
        for j=1:noinpvolcs
            alldates(noeqs+j)=volcs{j,3};
        end
    end
    
    % Permutate then assign dates back
    alldates=alldates(randperm(length(alldates)));
    for j=1:noeqs
        eqs{j,1}=alldates(j);
    end
    for j=1:novolcs
        volcs{j,3}=alldates(noeqs+j);
    end
end

% Permutate only eruption times
if MC_mode=='Q'
    
    % Special case for first simulation - extract observed eruption dates
    if isempty(alldates)
        alldates=zeros(novolcs,1);
        for j=1:novolcs
            alldates(j)=volcs{j,3};
        end
    end
    
    % Permutate dates then assign them back
    alldates=alldates(randperm(length(alldates)));
    for j=1:novolcs
        volcs{j,3}=alldates(j);
    end
end

% Randomise eruption times with recurrence limit
if MC_mode=='X'
    
    % If first simulation, calculate the number of individual volcanoes,
    % how many eruptions they have
    if isempty(no_indvolcs)
        % Make array to store extract individual volcano numbers, the
        % indeces of their eruptions in the original (current) volcs
        % file, and their total numbers of eruptions
        indvolcnos = zeros(300,1); % volcano numbers only
        volcindeces = zeros(novolcs,1); % order extracted from volcs array (needed for reordering eqvolcsit)
        indvolcs = cell(300,3); % volcano number, indeces of eruptions, total eruptions
        for i = 1 : size(indvolcs,1)
            indvolcs{i,2} = []; % Empty array to append erution indeces too
        end
        % Individual volcano counter
        no_indvolcs = 1;
        % Eruption counter
        no_erups = 1;
        
        % Loop through input file
        for i = 1 : novolcs
            % Test to see if this volcano already extracted
            cond = indvolcnos == volcs{i,1};
            if isempty(find(cond,1))
                % Extract volcano number and first eruption index
                indvolcnos(no_indvolcs) = volcs{i,1};                
                indvolcs{no_indvolcs,1} = volcs{i,1};
                indvolcs{no_indvolcs,2} = [indvolcs{no_indvolcs,2} i];
                volcindeces(no_erups) = i;
                no_erups = no_erups + 1;
                % Determine the number of total eruptions from this volcano
                % and their indeces
                indvolcs{no_indvolcs,3} = 1;
                for j = i+1 : novolcs
                    if volcs{j,1} == indvolcs{no_indvolcs,1}
                        indvolcs{no_indvolcs,2} = [indvolcs{no_indvolcs,2} j];
                        indvolcs{no_indvolcs,3} = indvolcs{no_indvolcs,3} + 1;
                        volcindeces(no_erups) = j;
                        no_erups = no_erups + 1;
                    end
                end
                % Add to counter
                no_indvolcs = no_indvolcs + 1;
            end
        end
        
        % Trim empty data
        no_indvolcs = no_indvolcs - 1;
        indvolcs = indvolcs(~cellfun('isempty',indvolcs(:,1)),:);
        
        % Make new volcs file with all the eruptions from one volcano next
        % to each other to allow easier randomising with recurrence limits
        volcsnew = cell(size(volcs));
        % Start counter 
        loc = 1;
        % Fill the new volcs file in order of volcano
        for i = 1 : no_indvolcs
            for j = 1 : indvolcs{i,3}
                volcsnew{loc,1} = volcs{indvolcs{i,2}(j),1};
                volcsnew{loc,2} = volcs{indvolcs{i,2}(j),2};
                volcsnew{loc,3} = volcs{indvolcs{i,2}(j),3};
                volcsnew{loc,4} = volcs{indvolcs{i,2}(j),4};
                volcsnew{loc,5} = volcs{indvolcs{i,2}(j),5};
                volcsnew{loc,6} = volcs{indvolcs{i,2}(j),6};
                volcsnew{loc,7} = volcs{indvolcs{i,2}(j),7};
                loc = loc + 1;
            end
        end
        % Overwrite
        volcs = volcsnew;
        % Also reorder the eqvolcdist array to match the new volcs array
        eqvolcdist = eqvolcdist(:,volcindeces);
    end   
    
    % Randomise eruption times, one volcano at a time to enforce recurrence
    % limit
    erup_no = 1; % Eruption counter for volcs file
    for i = 1 : no_indvolcs
        % Generate empty array of randomised dates for ith volcano
        dates = zeros(1,indvolcs{i,3});
        % Generate random times for each eruption at ith volcano
        for j = 1 : indvolcs{i,3}
            % First eruption
            if j == 1
                % Randomise
                dates(j) = round(last_volc-rand*(last_volc-first_volc));
                % Assign to volcs file
                volcs{erup_no,3} =  dates(j);
                erup_no = erup_no + 1;
            % Multiple eruptions
            else
                reclim = 1; % Recurrence limit test value
                % Try random dates until one is allowed
                while isempty(reclim) == 0
                    % Try random date
                    dates(j) = round(last_volc-rand*(last_volc-first_volc));
                    % Calculate recurrence times with other eruptions at
                    % ith volcano
                    rectest = zeros(1,j-1);
                    for k = 1 : j-1
                        rectest(k) = dates(j) - dates(k);
                    end
                    % Test to see if recurrence times okay
                    cond = abs(rectest) <= erup_lim;
                    reclim = find(cond,1);
                end
                % Assign eruption to volcs if recurrence times okay
                volcs{erup_no,3} =  dates(j);
                erup_no = erup_no + 1;
            end
        end
    end
end

% Permutate only eruption times with recurrence limit
if MC_mode=='Y'
    
    % Special case for first simulation - extract observed eruption dates
    if isempty(alldates)
        alldates=zeros(novolcs,1);
        for j=1:novolcs
            alldates(j)=volcs{j,3};
        end
    end
    
    % If first simulation, calculate the number of individual volcanoes,
    % how many eruptions they have, and reorder volcs file to allow for
    % easier eruption permutating with recurrence limits
    if isempty(no_indvolcs)
        % Make array to store extract individual volcano numbers, the
        % indeces of their eruptions in the original (current) volcs
        % file, and their total numbers of eruptions
        indvolcnos = zeros(300,1); % volcano numbers only
        volcindeces = zeros(novolcs,1); % order extracted from volcs array (needed for reordering eqvolcsit)
        indvolcs = cell(300,3); % volcano number, indeces of eruptions, total eruptions
        for i = 1 : size(indvolcs,1)
            indvolcs{i,2} = []; % Empty array to append erution indeces too
        end
        % Individual volcano counter
        no_indvolcs = 1;
        % Eruption counter
        no_erups = 1;
        
        % Loop through input file
        for i = 1 : novolcs
            % Test to see if this volcano already extracted
            cond = indvolcnos == volcs{i,1};
            if isempty(find(cond,1))
                % Extract volcano number and first eruption index
                indvolcnos(no_indvolcs) = volcs{i,1};
                indvolcs{no_indvolcs,1} = volcs{i,1};
                indvolcs{no_indvolcs,2} = [indvolcs{no_indvolcs,2} i];
                volcindeces(no_erups) = i;
                no_erups = no_erups + 1;
                % Determine the number of total eruptions from this volcano
                % and their indeces
                indvolcs{no_indvolcs,3} = 1;
                for j = i+1 : novolcs
                    if volcs{j,1} == indvolcs{no_indvolcs,1}
                        indvolcs{no_indvolcs,2} = [indvolcs{no_indvolcs,2} j];
                        indvolcs{no_indvolcs,3} = indvolcs{no_indvolcs,3} + 1;
                        volcindeces(no_erups) = j;
                        no_erups = no_erups + 1;
                    end
                end
                % Add to counter
                no_indvolcs = no_indvolcs + 1;
            end
        end
        
        % Trim empty data
        no_indvolcs = no_indvolcs - 1;
        indvolcs = indvolcs(~cellfun('isempty',indvolcs(:,1)),:);
        
        % Sort in order of decreasing number of eruptions to allow easier
        % permuating (i.e. permutate eruption dates to the volcanoes with
        % the most eruptions first to prevent them being rejected because
        % of recurrence limit)
        indvolcs = sortrows(indvolcs,3,'descend');
        
        % Make new volcs file with all the eruptions from one volcano next
        % to each other to allow easier permutating with recurrence limits
        volcsnew = cell(size(volcs));
        % Start counter 
        loc = 1;
        % Fill the new volcs file in order of volcano
        for i = 1 : no_indvolcs
            for j = 1 : indvolcs{i,3}
                volcsnew{loc,1} = volcs{indvolcs{i,2}(j),1};
                volcsnew{loc,2} = volcs{indvolcs{i,2}(j),2};
                volcsnew{loc,3} = volcs{indvolcs{i,2}(j),3};
                volcsnew{loc,4} = volcs{indvolcs{i,2}(j),4};
                volcsnew{loc,5} = volcs{indvolcs{i,2}(j),5};
                volcsnew{loc,6} = volcs{indvolcs{i,2}(j),6};
                volcsnew{loc,7} = volcs{indvolcs{i,2}(j),7};
                loc = loc + 1;
            end
        end
        % Overwrite
        volcs = volcsnew;
        % Also reorder the eqvolcdist array to match the new volcs array
        eqvolcdist = eqvolcdist(:,volcindeces);
    end
    
    % Permutate eruption times, one volcano at a time to enforce recurrence
    % limit
    erup_no = 1; % Eruption counter for volcs file
    alldates1 = alldates; % Pool of observed dates to be randomly selected from (permutated)
    for i = 1 : no_indvolcs
        % Generate empty array of randomly permutated dates for ith volcano
        dates = zeros(1,indvolcs{i,3});
        % Select random permutated dates for each eruption at ith volcano
        for j = 1 : indvolcs{i,3}
            % First eruption
            if j == 1
                % Select random date from array
                dateselect = randi(length(alldates1));
                dates(j) = alldates1(dateselect);
                % Assign to volcs file
                volcs{erup_no,3} =  dates(j);
                erup_no = erup_no + 1;
                % Remove selected date from the pool
                alldates1(dateselect) = [];
            % Multiple eruptions
            else
                reclim = 1; % Recurrence limit test value
                % Try random dates until one is allowed
                while isempty(reclim) == 0
                    % Try random date
                    dateselect = randi(length(alldates1));
                    dates(j) = alldates1(dateselect);
                    % Calculate recurrence times with other eruptions at
                    % ith volcano
                    rectest = zeros(1,j-1);
                    for k = 1 : j-1
                        rectest(k) = dates(j) - dates(k);
                    end
                    % Test to see if recurrence times okay
                    cond = abs(rectest) <= erup_lim;
                    reclim = find(cond,1);
                end
                % Assign eruption to volcs if recurrence times okay
                volcs{erup_no,3} =  dates(j);
                erup_no = erup_no + 1;
                % Remove selected date from the pool
                alldates1(dateselect) = [];
            end
        end
    end   
end

% Randomise eruption times at each individual volcano by permutating
% recurrence times from that volcano
if MC_mode == 'Z'
    
    % If first simulation, calculate the number of individual volcanoes,
    % how many eruptions they have, recurrence times, and reorder volcs
    % file to allow for easier eruption randomising with recurrence limits
    if isempty(no_indvolcs)
        % Make array to store extract individual volcano numbers, the
        % indeces of their eruptions in the original (current) volcs
        % file, and their total numbers of eruptions
        indvolcnos = zeros(300,1); % volcano numbers only
        volcindeces = zeros(novolcs,1); % order extracted from volcs array (needed for reordering eqvolcsit)
        indvolcs = cell(300,4); % volcano number, indeces of eruptions, total eruptions, recurrence times
        for i = 1 : size(indvolcs,1)
            indvolcs{i,2} = []; % Empty array to append erution indeces to
            indvolcs{i,4} = []; % Empty array to append recurrence times to
        end
        % Individual volcano counter
        no_indvolcs = 1;
        % Eruption counter
        no_erups = 1;
        
        % Loop through input file
        for i = 1 : novolcs
            % Test to see if this volcano already extracted
            cond = indvolcnos == volcs{i,1};
            if isempty(find(cond,1))
                % Extract volcano number and first eruption index
                indvolcnos(no_indvolcs) = volcs{i,1};
                indvolcs{no_indvolcs,1} = volcs{i,1};
                indvolcs{no_indvolcs,2} = [indvolcs{no_indvolcs,2} i];
                volcindeces(no_erups) = i;
                no_erups = no_erups + 1;
                % Determine the number of total eruptions from this volcano
                % and their indeces and recurrence times
                indvolcs{no_indvolcs,3} = 1;
                for j = i+1 : novolcs
                    if volcs{j,1} == indvolcs{no_indvolcs,1}
                        indvolcs{no_indvolcs,2} = [indvolcs{no_indvolcs,2} j];
                        indvolcs{no_indvolcs,4} = [indvolcs{no_indvolcs,4} (volcs{indvolcs{no_indvolcs,2}(indvolcs{no_indvolcs,3}),3})-volcs{j,3}];
                        indvolcs{no_indvolcs,3} = indvolcs{no_indvolcs,3} + 1;     
                        volcindeces(no_erups) = j;
                        no_erups = no_erups + 1;
                    end
                end
                % Add to counter
                no_indvolcs = no_indvolcs + 1;
            end
        end
        
        % Trim empty data
        no_indvolcs = no_indvolcs - 1;
        indvolcs = indvolcs(~cellfun('isempty',indvolcs(:,1)),:);
        
        % Make new volcs file with all the eruptions from one volcano next
        % to each other to allow easier randomising with recurrence limits
        volcsnew = cell(size(volcs));
        % Start counter 
        loc = 1;
        % Fill the new volcs file in order of volcano
        for i = 1 : no_indvolcs
            for j = 1 : indvolcs{i,3}
                volcsnew{loc,1} = volcs{indvolcs{i,2}(j),1};
                volcsnew{loc,2} = volcs{indvolcs{i,2}(j),2};
                volcsnew{loc,3} = volcs{indvolcs{i,2}(j),3};
                volcsnew{loc,4} = volcs{indvolcs{i,2}(j),4};
                volcsnew{loc,5} = volcs{indvolcs{i,2}(j),5};
                volcsnew{loc,6} = volcs{indvolcs{i,2}(j),6};
                volcsnew{loc,7} = volcs{indvolcs{i,2}(j),7};
                loc = loc + 1;
            end
        end
        % Overwrite
        volcs = volcsnew;
        % Also reorder the eqvolcdist array to match the new volcs array
        eqvolcdist = eqvolcdist(:,volcindeces);
        
    end
    
    % Randomise eruption times at each volcano, by assigning a random(ish)
    % first date, then permutating the recurrence times from that volcano
    erup_no = 1; % Eruption counter for volcs file
    for i = 1 : no_indvolcs
        % Generate empty array of random dates for ith volcano
        dates = zeros(1,indvolcs{i,3});
        % Select random permutated dates for each eruption at ith volcano
        for j = 1 : indvolcs{i,3}
            % First eruption
            if j == 1
                % Special case for only 1 eruption recorded, therefore any
                % random date allowed
                if indvolcs{i,3} == 1
                    dates(j) = round(last_volc-rand*(last_volc-first_volc));
                    volcs{erup_no,3} =  dates(j);
                    erup_no = erup_no + 1;
                % Volcano has multiple eruptions
                else
                    % Determine the allowable time range for the first eruption
                    % based on the total recurrence time
                    total_rec = sum(indvolcs{i,4});
                    max_firsterup = last_volc - total_rec;
                    % Randomise within these limits
                    dates(j) = round(max_firsterup-rand*(max_firsterup-first_volc));
                    volcs{erup_no,3} =  dates(j);
                    erup_no = erup_no + 1;
                    % Extract the subsequent recurrence times
                    rec_times = indvolcs{i,4};
                end
            % Multiple eruptions
            else
                % Select random date from array
                recselect = randi(length(rec_times));
                dates(j) = dates(j-1) + rec_times(recselect);
                % Assign to volcs file
                volcs{erup_no,3} =  dates(j);
                erup_no = erup_no + 1;
                % Remove selected date from the pool
                rec_times(recselect) = [];
            end
        end
    end
    
end

end