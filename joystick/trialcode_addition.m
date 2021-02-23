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

% only add trial codes for the current participant

function saloglo_input=trialcode_addition(saloglo_input,partnum,blocks)
    % loop through every block for each participant, except for test block
    % or block with no deivance, which is block 12. 
    for blocknum=blocks
        % give the first digit of trialcode - blocktype
        if saloglo_input{partnum,blocknum}.blocktype==0 
            % alternatively, I can count to determine the block type
            for trlnum=1:length(saloglo_input{partnum,blocknum}.stim)
                saloglo_input{partnum,blocknum}.stim(trlnum).trialcode='0';
            end
        elseif saloglo_input{partnum,blocknum}.blocktype==1
            for trlnum=1:length(saloglo_input{partnum,blocknum}.stim)
                saloglo_input{partnum,blocknum}.stim(trlnum).trialcode='1';
            end
        else
            error('Wrong blocktype.')
        end
        
        % give the second digit of trialcode - salience 
        % you can deduct the actual audio played by blocktype * salience 
        % loop through every trial 
        for trialnum=1:length(saloglo_input{partnum,blocknum}.stim)
            switch saloglo_input{partnum,blocknum}.stim(trialnum).trialname
                case 'glostd'
                    saloglo_input{partnum,blocknum}.stim(trialnum).trialcode(2)='0';
                case 'glodev'
                    if ismember(saloglo_input{partnum,blocknum}.stim(trialnum).trialans(3:4),['X1' 'X2'])
                        saloglo_input{partnum,blocknum}.stim(trialnum).trialcode(2)='1';
                    elseif strcmp(saloglo_input{partnum,blocknum}.stim(trialnum).trialans(3:4),'Y2')
                        saloglo_input{partnum,blocknum}.stim(trialnum).trialcode(2)='2';
                    else 
                        error('Wrong trialans.')
                    end
                otherwise
                    error('Wrong trialname.')
            end
        end
        
        % give the third digit of trialcode - lateness 
        % give 0 to the first 15 trials
        if ismember(blocknum,[1 12])
            for trialnum=1:length(saloglo_input{partnum,blocknum}.stim)
                saloglo_input{partnum,blocknum}.stim(trialnum).trialcode(3)='0';
            end
        else
            for trialnum=1:15
                saloglo_input{partnum,blocknum}.stim(trialnum).trialcode(3)='0';
            end

            % for the following trials, count lateness in standards from 1
            count=0;
            for trialnum=16:length(saloglo_input{partnum,blocknum}.stim)
                if strcmp(saloglo_input{partnum,blocknum}.stim(trialnum).trialname,'glostd')
                    count=count+1;
                elseif strcmp(saloglo_input{partnum,blocknum}.stim(trialnum).trialname,'glodev')
                    count=0;
                else 
                    error('Wrong trialname.')
                end
                saloglo_input{partnum,blocknum}.stim(trialnum).trialcode(3)=num2str(count);

                % throw error is count is too big. 
                if count>5
                    error('Your lateness count is too big.')
                end
            end
        end
    end

end

