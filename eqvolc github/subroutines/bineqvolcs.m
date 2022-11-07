function [eqs,volcs,seismoment]=bineqvolcs(no_years,datebounds,inpeqs,noinpeqs,inpvolcs,noinpvolcs,mult_erups,vmode)

% Assigns earthquakes and eruptions to their yearly bins

% Generate vectors to store the binned earthquakes and eruptions
eqs=zeros(1,no_years);
volcs=zeros(1,no_years);
% Generate vector to store cumulative seismic moment each year
seismoment=zeros(1,no_years);
% Generate vector to store cumulative VEI each year
vei=zeros(1,no_years);
% 1000 eruptions in each year bin preallocated - alter if more needed!
volcdeets=cell(1,no_years);
for i=1:size(volcdeets,2)
    volcdeets{1,i}=cell(100,2);
end

% Loop through earthquakes, assigning each to the correct bin
for i=1:noinpeqs
    % Find the correct bin
    bin_no=find(inpeqs{i,1}<datebounds,1)-1;
    % Add to counter for bin
    eqs(bin_no)=eqs(bin_no)+1;
    % Add to seismic moment for bin
    seismoment(bin_no)=seismoment(bin_no)+10^(1.5*(inpeqs{i,5}+6.06));
end

% Loop through eruptions, assigning each to correct bin
for i=1:noinpvolcs
    % Find the correct bin
    bin_no=find(inpvolcs{i,3}<datebounds,1)-1;
    % Determine the eruption number in the timeseries bin
    row_no=find(cellfun('isempty',volcdeets{1,bin_no}(:,1))==1,1);
    
    % Multiple eruptions allowed or not
    if mult_erups=='Y'
        volcdeets{bin_no}{row_no,1}=datestr(inpvolcs{i,3},'dd/mm/yyyy'); % Eruption date
        volcdeets{bin_no}{row_no,2}=inpvolcs{i,1}; % Volcano ID
        if vmode=='A'
            vei(bin_no)=vei(bin_no)+inpvolcs{i,4};
        elseif vmode=='E'
            vei(bin_no)=vei(bin_no)+10^(inpvolcs{i,4}-5); % VEI 5 erup =1, 4 =0.1, 3=0.01 etc
        end
    elseif mult_erups=='N'
        % Does volcano already have an eruption in timebin
        if row_no==1 || max(inpvolcs{i,1}==cell2mat(volcdeets{1,bin_no}(:,2)))==0
            volcdeets{bin_no}{row_no,1}=datestr(inpvolcs{i,3},'dd/mm/yyyy'); % Eruption date
            volcdeets{bin_no}{row_no,2}=inpvolcs{i,1}; % Volcano ID
            if vmode=='A'
                vei(bin_no)=vei(bin_no)+inpvolcs{i,4};
            elseif vmode=='E'
                vei(bin_no)=vei(bin_no)+10^(inpvolcs{i,4}-5); % VEI 5 erup =1, 4 =0.1, 3=0.01 etc
            end
        end
    end
end

% Trim empty volcano detail cells and add up the eruptions
for i=1:no_years
    volcdeets{1,i}=volcdeets{1,i}(~cellfun('isempty',volcdeets{1,i}(:,1)),:);
    volcs(i)=size(volcdeets{1,i},1);
end

% If using VEI scale for output
if vmode=='A'
    volcs=vei;
elseif vmode=='E'
    volcs=log(vei);
end

end