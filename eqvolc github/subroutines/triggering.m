function [E_before,E_after,E_average]=triggering(mind,maxd,mint,maxt,eqdate,last_eq,eqvolcdist,volcs,mult_erups,sameday_erups,novolcs,volc_length,avgEr,AvEx)

% Calculates number of eruptions before and after input earthquake for
% given filters and inputs. Also calculates the average eruption rate over
% the specified time period within the specified distance of input
% earthquake


% Generate the eruption date boundaries before and after the earthquake
t_before1=eqdate-maxt;
t_before2=eqdate-mint;
t_after1=eqdate+mint;
t_after2=eqdate+maxt;

% Generate cell to store the eruptions before and after and average 
E_before=cell(100,6);
E_after=cell(100,6);
E_average=cell(200,2);
% Start counters for before and after and average eruption numbers
n_before=1;
n_after=1;
n_total=1; % AER1

% Extract eruptions occuring within the specified distance from earthquake
for i=1:novolcs
    
    % If eruption occurs within the specified distance
    if eqvolcdist(i)>=mind && eqvolcdist(i)<maxd
        
        % Add to average eruption rate if appropriate:
        % Eruptions within certain time from earthquake excluded from
        % baseline rate calculation or not?
        if avgEr == 'A' || ( avgEr == 'E' && abs(eqdate-volcs{i,3})>AvEx )
            % Add to AER1 array
            % Multiple eruptions from same volcano allowed or not
            if mult_erups=='Y'
                E_average{n_total,1}=volcs{i,3}; % Eruption date
                E_average{n_total,2}=volcs{i,1}; % Volcano ID
                n_total=n_total+1;
            elseif mult_erups=='N'
                % Volcano in question does not already have an eruption in AER1
                if n_total==1 || max(volcs{i,1}==cell2mat(E_average(:,2)))==0
                    E_average{n_total,1}=volcs{i,3}; % Eruption date
                    E_average{n_total,2}=volcs{i,1}; % Volcano ID
                    n_total=n_total+1;
                    % Volcano already has an eruption in AER1
                else
                    % Calculate time between the eruptions
                    cond=find(cell2mat(E_average(:,2))==volcs{i,1});
                    timebtw=zeros(1,length(cond));
                    for j=1:length(cond)
                        timebtw(j)=abs(E_average{cond(j),1}-volcs{i,3});
                    end
                    % Add to AER1 if sufficent time gap between eruptions
                    if timebtw>(maxt-mint)
                        E_average{n_total,1}=volcs{i,3}; % Eruption date
                        E_average{n_total,2}=volcs{i,1}; % Volcano ID
                        n_total=n_total+1;
                    end
                end
            end
        end
        
        % Eruption located in the specified time period before earthquake
        if volcs{i,3}>=t_before1 && volcs{i,3}<t_before2
            % Multiple eruptions allowed 
            if mult_erups=='Y'
                E_before{n_before,1}=volcs{i,3}; % Eruption date
                E_before{n_before,2}=volcs{i,5}; % Volcano lat
                E_before{n_before,3}=volcs{i,6}; % Volcano lon
                E_before{n_before,4}=volcs{i,1}; % Volcano ID
                E_before{n_before,5}=volcs{i,2}; % Volcano name
                E_before{n_before,6}=volcs{i,4}; % VEI
                n_before=n_before+1;
            % Multiple eruptions not allowed
            elseif mult_erups=='N'
                % Does volcano already have a before eruption 
                if n_before==1 || max(volcs{i,1}==cell2mat(E_before(:,4)))==0
                    E_before{n_before,1}=volcs{i,3}; % Eruption date
                    E_before{n_before,2}=volcs{i,5}; % Volcano lat
                    E_before{n_before,3}=volcs{i,6}; % Volcano lon
                    E_before{n_before,4}=volcs{i,1}; % Volcano ID
                    E_before{n_before,5}=volcs{i,2}; % Volcano name
                    E_before{n_before,6}=volcs{i,4}; % VEI
                    n_before=n_before+1;
                end
            end
            
        % Eruption located in the specified time period after earthquake
        elseif (volcs{i,3}>t_after1 && volcs{i,3}<=t_after2) || (sameday_erups=='A' && mint==0 && (volcs{i,3}>=t_after1 && volcs{i,3}<=t_after2))
            % Multiple eruptions allowed 
            if mult_erups=='Y'
                E_after{n_after,1}=volcs{i,3}; % Eruption date
                E_after{n_after,2}=volcs{i,5}; % Volcano lat
                E_after{n_after,3}=volcs{i,6}; % Volcano lon
                E_after{n_after,4}=volcs{i,1}; % Volcano ID
                E_after{n_after,5}=volcs{i,2}; % Volcano name
                E_after{n_after,6}=volcs{i,4}; % VEI
                n_after=n_after+1;
            % Multiple eruptions not allowed
            elseif mult_erups=='N'
                % Does volcano already have a before eruption 
                if n_after==1 || max(volcs{i,1}==cell2mat(E_after(:,4)))==0
                    E_after{n_after,1}=volcs{i,3}; % Eruption date
                    E_after{n_after,2}=volcs{i,5}; % Volcano lat
                    E_after{n_after,3}=volcs{i,6}; % Volcano lon
                    E_after{n_after,4}=volcs{i,1}; % Volcano ID
                    E_after{n_after,5}=volcs{i,2}; % Volcano name
                    E_after{n_after,6}=volcs{i,4}; % VEI
                    n_after=n_after+1;
                end
            end
        end
    end        
end

% Trim empty before and after cell values
E_before=E_before(~cellfun('isempty',E_before(:,1)),:);
E_after=E_after(~cellfun('isempty',E_after(:,1)),:);
E_average=E_average(~cellfun('isempty',E_average(:,1)),:);

% Calculate AER1
if avgEr == 'A'
    E_average=(maxt-mint)*((size(E_average,1))/(volc_length)); % AER1 over entire eruption catalogue
elseif avgEr == 'E'
    % Determine whether earthquake occurs within AvEx of end of
    % time-serires
    endgap = last_eq - eqdate;
    if endgap < AvEx
        eqbuffer = AvEx + endgap;
    else
        eqbuffer = 2 * AvEx;
    end
    E_average=(maxt-mint)*((size(E_average,1))/(volc_length-eqbuffer)); % AER1 over entire eruption catalogue, excluding period around earthquake
end

end