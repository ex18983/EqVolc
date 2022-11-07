function cmt2input

% Converts the gCMT .ndk database format to a usable input spreadsheet file

% INPUT

input_file = 'gCMTall2.txt'; % Name of input CMT data as text file (ndk format)
output_name = 'gCMTmw7_2.xlsx'; % Desired name of output earthquake events file for EqVolc
min_mw = 0; % Minimum desired Mw of output events
min_depth = 0; % Minumum depth of output events (km?)
max_depth = 9999; % Maximum depth of output events (km?)

% Define the home directory for function
home_dir= 'C:\Users\ex18983\OneDrive - University of Bristol\Documents\MATLAB\eqvolc github';
cd(home_dir)

% CALCULATION

% Open the CMT input file
fid_min=fopen(input_file,'r');

% Initiate an event counter
event_no=2;
% Initiate break point
nexteq='Y';

% Create an array to store the 5 lines CMT detail for each event
event=cell(5,1);
% Preallocate array to store selected output events (max 10000 - can change)
output=cell(10000,13);
% Write output headers
output{1,1}='Date';
output{1,2}='Lat';
output{1,3}='Lon';
output{1,4}='Depth';
output{1,5}='Mw';
output{1,6}='Earthquake';
output{1,7}='Ranking';
output{1,8}='SlipType';
output{1,9}='MagType';

% % Multiple slip type category checker
% counts=0;
% countsmw=zeros(55,1);
% countsd=zeros(55,1);

% Loop through events
while nexteq=='Y'
    
    % Read in the next event
    for j=1:5
        event{j}=fgetl(fid_min);
    end
    
    % Make sure event exists
    if isa(event{1},'char')==1
        % Event exists - calculate event magnitude
        exp=str2double(event{4}(1:2));
        moment=str2double(event{5}(50:56));
        moment=moment*10^exp;
        mw=(2/3)*(log10(moment)-16.1); % See cmt webpage for equation. Can differ by 0.1 prior to 2006
        mw=round(mw,2);
        
        % Moment magntiude filter
        if mw>=min_mw
            % Depth filter 
            depth=str2double(event{3}(49:53));
            if depth>min_depth && depth<=max_depth
                
                % Determine slip type and remove multiple slip type events
                rake1=str2double(event{5}(65:68));
                rake2=str2double(event{5}(77:80));
                if rake1<=-160 || rake2<=-160
                    type='strike-slip'; % Probably a strike-slip eq
                    % Check for multiple categories
                    if (rake1>=-110 && rake1<=-70) || (rake2>=-110 && rake2<=-70) || (rake1>=70 && rake1<=110) || (rake2>=70 && rake2<=110)
                        type='unknown'; % Falls into multiple categories
%                         counts=counts+1;
%                         countsmw(counts)=mw;
%                         countsd(counts)=depth;
                    end
                elseif (rake1>=-110 && rake1<=-70) || (rake2>=-110 && rake2<=-70)
                    type='normal'; % Probably a normal eq
                    % Check for multiple categories
                    if (rake1>=70 && rake1<=110) || (rake2>=70 && rake2<=110)
                        type='unknown'; % Falls into multiple categories
%                         counts=counts+1;
%                         countsmw(counts)=mw;
%                         countsd(counts)=depth;
                    end
                elseif (rake1>=-20 && rake1<=20) || (rake2>=-20 && rake2<=20)
                    type='strike-slip'; % Probably a strike-slip eq
                    % Check for multiple categories
                    if (rake1>=-110 && rake1<=-70) || (rake2>=-110 && rake2<=-70) || (rake1>=70 && rake1<=110) || (rake2>=70 && rake2<=110)
                        type='unknown'; % Falls into multiple categories
%                         counts=counts+1;
%                         countsmw(counts)=mw;
%                         countsd(counts)=depth;
                    end
                elseif (rake1>=70 && rake1<=110) || (rake2>=70 && rake2<=110)
                    type='reverse'; % Probably a reverse eq
                elseif rake1>=160 || rake2>=160
                    type='strike-slip'; % Probably a strike-slip eq
                    % Check for multiple categories
                    if (rake1>=-110 && rake1<=-70) || (rake2>=-110 && rake2<=-70) || (rake1>=70 && rake1<=110) || (rake2>=70 && rake2<=110)
                        type='unknown'; % Falls into multiple categories
%                         counts=counts+1;
%                         countsmw(counts)=mw;
%                         countsd(counts)=depth;
                    end
                elseif rake1<0 && rake2<0
                    type='oblique-normal'; % Probably oblique-normal eq
                elseif rake1>0 && rake2>0
                    type='oblique-reverse'; % Probably oblique-reverse eq
                end
                
                % Filters satisfied - fill output event details
                day=event{1}(14:15);
                month=event{1}(11:12);
                year=event{1}(6:9);
                output{event_no,1}=strcat(day,'/',month,'/',year); %date
                output{event_no,2}=str2double(event{3}(24:29)); %lat
                output{event_no,3}=str2double(event{3}(36:42)); %lon
                output{event_no,4}=depth; %depth
                output{event_no,5}=mw; %Mw
                output{event_no,6}=event{1}(57:80); %location
                output{event_no,7}=10; %ranking: (10=from gCMT input file)      
                output{event_no,8}=type; %type
                output{event_no,9}='mw'; %Mw
                % Add to event counter
                event_no=event_no+1;
            end
        end
        
    else % Event doesn't exist - break while loop
        nexteq='N';
    end
    
end

% Trim empty output cell values
output=output(~cellfun('isempty',output(:,1)),:);

% Write output file, close all
overwr='Y';
cd('eq databases')
[overwr]=overwrite_check(output_name,overwr);
if overwr=='Y'
    xlswrite(output_name,output)
end
fclose('all');

% Return to home
cd(home_dir)

end