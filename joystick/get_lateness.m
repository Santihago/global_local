%% Get lateness

addpath('D:\Dropbox\MyScience\MyPostdoc\exp\global_local\1-exp');
addpath('D:\Dropbox\MyScience\MyPostdoc\exp\global_local\3-analysis\eeg\processing\0-yijie-preprocessing-scripts\yijie-preprocessing')

% specify participant number
partnum=input('Please insert the number of participant(s): ');
if ~ismember(partnum,1:21)
    partnum=input('Please insert the number of participant(s): ');
end

% load all trials
load('saloglo_input.mat');

for id = 1:length(partnum)
    % add latenesses to the trials
    % this function adds trial codes to trials.
    % first digit: block type;
    % second digit: salience:
    % global standard is always 0;
    % in block type 1 (glostd=locdev):
    % X1=1 - change in ear;
    % Y2=2 - change in pitch;
    % third digit: count of lateness:
    % test block (block 1) and block with no deviance (block 12) are
    % excluded;
    % first 15 trials of each block is given 0;
    % each global deviant has value 0.
    
    % correct for id in input matrix : id=21 was given stim for id=1, 
    % and id=14 was given stim for id=13)
    if partnum(id)==14
        j = 13;
    elseif partnum(id)==21
        j = 1;
    else
        j = id;
    end
    
    % Obtain trialcodes with the custom function
    trial_info = trialcode_addition(saloglo_input, partnum(j), [7 8 10 11]);
    
    %     %reshape table with each cell's value
    
    for blkn = [7:8,10:11]
        
        trial_n = length(trial_info{partnum(j), blkn}.stim);
        
        out_table = cell2table(cell(trial_n,9), 'VariableNames', ...
        {'id', 'block', 'trial', 'blocktype', 'deviance', 'laterality', ...
        'stim', 'trialtype', 'trialcode'});
        
        % 1. Create empty output table with id in first column
        out_table.id(1:trial_n) = num2cell(repelem(id, trial_n)');
        % 2. Add trial number
        out_table.trial(1:trial_n) = num2cell(1:trial_n);
        % 3. Add block number
        out_table.block(1:trial_n) = num2cell(repelem(blkn, trial_n)');
        % 3. Blocktype
        blocktype = trial_info{partnum(j), blkn}.blocktype;
        out_table.blocktype(1:trial_n) = num2cell(repelem(blocktype, trial_n)');
        % 4. Deviance
        deviance = trial_info{partnum(j), blkn}.deviance;
        out_table.deviance(1:trial_n) = num2cell(repelem(deviance, trial_n)');
        % 5. Laterality
        laterality = trial_info{partnum(j), blkn}.laterality;
        out_table.laterality(1:trial_n) = repelem({laterality}, trial_n)';
        % 6. Stim name
        out_table.stim(1:trial_n) = {trial_info{partnum(j), blkn}.stim.trialans};
        % 7. Trial type
        out_table.trialtype(1:trial_n) = {trial_info{partnum(j), blkn}.stim.trialname};
        % 8. Trialcode (3 digits code, to decompose later in R)
        out_table.trialcode(1:trial_n) = {trial_info{partnum(j), blkn}.stim.trialcode};
        
        
        %
        % ?why step of 'reduce the number of data points for earlier latenesses' in Yijie's code?
        %
        
        % Save as csv
        out_folder = 'D:\Dropbox\MyScience\MyPostdoc\exp\global_local\2-data\trial_info\';
        out_fname = strcat('trial_info_id_', num2str(id), '_block_', num2str(blkn), '.csv');
        fprintf(1, 'Now saving %s\n', out_fname);
        out_fullfilename = fullfile(out_folder, out_fname);
        writetable(out_table, out_fullfilename, 'Delimiter', ',')
        
    end
    
end
