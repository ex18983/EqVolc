function [overwr]=overwrite_check(file_name,overwr)

% Generic check to see if a file already exists in the output folder with a
% name matching the output file name
% If file exists, script will present user with options 'yes, no, all' for
% overwriting files

% Check if a file already exists with the output name
if isfile(file_name)==1
    % If so, ask permission to overwrite from user
    disp(['WARNING! File already exists with name "',file_name,'". Data will be over written.'])
    overwr=input('Continue? Y = yes, N = no, A = continue for all files: ','s');
    % If yes, continue
    if overwr=='Y' || overwr=='y'
        overwr='Y';
    % If no, stop and return to invoking function
    elseif overwr=='N' || overwr=='n'
        overwr='N';
        return
    % If all, continue and set overwr to A to prevent further checks
    elseif overwr=='A' || overwr=='a'
       overwr='A';
    % If invalid answer selected, throw error message and stop script
    elseif overwr~='A' && overwr~='a' && overwr~='Y' && overwr~='y' && overwr ~='N' && overwr~='n'
        disp('Invalid response, stopping script')
        overwr='N';
        return
    end
end
end