%% Open joystick files and output .csv
%author: S. Munoz Moldes
% Opens the .mat files and outputs them as csv with slight modification
% Matrices are inverted (columns become rows) and an ID column is added

clear all
close all

% Set folder containing files
in_folder = 'D:\Dropbox\MyScience\MyPostdoc\exp\global_local\2-data\joystick\';
% List of .mat files inside the folder
file_list = dir(fullfile(in_folder, '*.mat'));

%% Loop through each file
for k = 1:length(file_list)

    % File basename
    name = file_list(k).name;
    % Number of ID from the basename
    id = str2num(name(1:2));
    % Full filename
    fullfilename = fullfile(in_folder, name);
    % Open file
    fprintf(1, 'Now reading %s\n', name);
    this_mat = load(fullfilename);
    this_joystick_data = this_mat.pos';  %reverse matrix
    % Add ID number as column
    id_column = repelem(id, length(this_joystick_data))';
    final = [ id_column this_joystick_data];
    % Save as csv
    out_folder = 'D:\Dropbox\MyScience\MyPostdoc\exp\global_local\3-analysis\joystick\1-csv files\';
    out_fname = strcat(name(1:end-4), '.csv');
    out_fullfilename = fullfile(out_folder, out_fname);
    writematrix(final , out_fullfilename, 'Delimiter', ',')
    
    clear this_mat
end