clc         % clear command window
close all   % close all figures
clear       % clear workspace

%% settings 
settings.startfolder        = 'C:\Users\paulk\Documents\Programmieren\Matlab'; % input root folder
settings.fastmode           = true;    % improves speed by omitting real/imaginary fit plots
settings.flipdirection      = false;    % sweep through files in flipped direction
settings.autoanalysisfolder = true;    % automatically creats an analysis folder, select if false
settings.multimodes         = true;    % set true if you want to fit multiple modes within a file
settings.oldmode            = false;
settings.repickworkingrange = false;    % select new workingrange after finished fitting round
settings.multifit           = 70;   % X fits before selecting a new fit range
% Unused
settings.useactualtemp  = false;    % not used for anything so far

%% Do work
[settings,datainfo,loopinfo] = initialize(settings);
[settings,datainfo] = get_filenames_and_detect_experiment(settings,datainfo);
datainfo = sort_filenames(settings,datainfo);
if settings.oldmode
    [datainfo,loopinfo] = oldmode(datainfo,loopinfo);
else
    datainfo = extract_filename_content(settings,datainfo);
    loopinfo = select_modes(settings,datainfo);
end

disp('Job done')

%% Functions

% Initialize parameters
function [settings,datainfo,loopinfo] = initialize(settings)
    % Multifit names
    if settings.multifit == 1
        settings.multifitname = 'fit';
    else
        settings.multifitname = 'fits';
    end

    % Datapath and parentfolder
    datainfo.datapath=uigetdir(settings.startfolder,'Select data folder');
    if datainfo.datapath == 0
        disp('Canceled data folder selection.')
        return
    end
    datainfo.previousdatapath = fileparts(datainfo.datapath);

    % Analysis folder
    if ~settings.autoanalysisfolder % Manually select analysis folder
        datainfo.analysispath=uigetdir(datainfo.previousdatapath,'Select analysis folder');
        if datainfo.analysispath == 0
            disp('Canceled analysis folder selection.')
            return
        end
    else
        analysisFolderPath = fullfile(datainfo.previousdatapath, 'Analysis');
        if ~isfolder(analysisFolderPath)
            mkdir(analysisFolderPath);
        end
        datainfo.analysispath = analysisFolderPath;
    end

    % Initialize loopinfo parameters
    loopinfo.usename    = 1;
    loopinfo.usesweeper = 1;
    loopinfo.sweeper    = 1;
    loopinfo.fixedtemp  = 1;
    loopinfo.fixedfield = 1;
    loopinfo.fixedpower = 1;
end

% Add filenames & detect experiment from filenames
function [settings,datainfo] = get_filenames_and_detect_experiment(settings,datainfo)
    % Get filenames
    datainfo.fileinfo = dir([datainfo.datapath '\' '*.dat']);
    datainfo.filename = {datainfo.fileinfo.name}';
    
    % experiment detection
    settings.vti = contains(datainfo.filename(1),'_setTemp=') && contains(datainfo.filename(1),'K_Field=') ...
        && contains(datainfo.filename(1),'T_Power=') && contains(datainfo.filename(1),'_actTemp=') ...
        && contains(datainfo.filename(1),'percent.dat');
    settings.dilutionfridge = contains(datainfo.filename(1),'_setTmc=') && contains(datainfo.filename(1),'K_actTsample=')...
        && contains(datainfo.filename(1),'K_B1(setB=') && contains(datainfo.filename(1),'dBm.dat');
    % case: errorexit
    if ~settings.vti && ~settings.dilutionfridge && ~settings.oldmode
        errordlg({'Tool is not yet adapted for the used data.','Please use the settings.oldmode.'},'Data Error');
        return
    end
end

% Sort filenames regarding datetime, also apply flipdirection
function datainfo = sort_filenames(settings,datainfo)
    % Check if "_" was used for something else than a delimiter
    try % This throws the error, if the output of split has varying dimensions, ADD: compare with expected # of delimiters
        split(datainfo.filename,'_');
    catch
        errordlg({'Apparently a "_" was used incorrectly in a filename.',['Please replace all "_" that are not used ...' ...
            'as delimiters with "-".'],['You can use the program "fix_underscore_in_oldfridge_filenames.m" in ...' ...
            ' "\\S4\Datenpool\Microwave Group\Programs".']},'Filename Error');
    end

    % Sort filenames regarding their creation date             
    [~, sortorder] = sort([datainfo.fileinfo.datenum]); % Sort based on creation dates
    datainfo.filename = datainfo.filename(sortorder);   % Update sorted filenames
    datainfo = rmfield(datainfo, 'fileinfo');           % Remove fileinfo substruct    
    
    % Flip file order according to settings
    if settings.flipdirection 
        datainfo.filename = flip(datainfo.filename);
    end
end

% Prepare Parameters for the oldmode case
function [datainfo, loopinfo] = oldmode(datainfo, loopinfo)
    loopinfo.currentfiles = datainfo.filename;  % Set all files to current files
    % Select sweeper file
    [datainfo.sweeperfile,datainfo.sweeperpath] = uigetfile('*.*','Choose sweeper file',datainfo.previousdatapath);
    loopinfo.analysissweeper = importdata([datainfo.sweeperpath '\' datainfo.sweeperfile]);
    loopinfo.analysissweeper = num2str(loopinfo.analysissweeper);
    loopinfo.analysissweeper = cellstr(loopinfo.analysissweeper);
    if settings.flipdirection % Flip sweeperfile
        loopinfo.analysissweeper = flip(loopinfo.analysissweeper);
    end
    datafoldername = split(datainfo.datapath,'\');
    datainfo.datafoldername = datafoldername(end);
    loopinfo.analysisfile   = ['analysis_XXX_' char(datainfo.datafoldername) '_step=' num2str(str2double(char(...
        loopinfo.analysissweeper(2)))-str2double(char(loopinfo.analysissweeper(1))))  '_range=' sprintf('%g',...
        str2double(loopinfo.analysissweeper(1))) '-' sprintf('%g',str2double(loopinfo.analysissweeper(end))) '.dat'];
end

% Extract information from filenames
function datainfo = extract_filename_content(settings,datainfo)
    split_file_names = split(datainfo.filename,'_');

    % Name
    if settings.vti
        datainfo.name = split_file_names(:,2);
    elseif settings.dilutionfridge
        datainfo.name = split_file_names(:,6); % Also adds the center frequency to the name
    end

    % Settemp
    if settings.vti || settings.dilutionfridge
        settemp_split = split(split_file_names(:,3),'=');
        settemp_split = split(settemp_split(:,2),'K');
    end
    datainfo.temp = settemp_split(:,1);

    % Field
    if settings.vti
        field_split = split(split_file_names(:,4),'=');
        field_split = split(field_split(:,2),'T');
    elseif settings.dilutionfridge
        field_split = split(split_file_names(:,5),'=');
        field_split = split(field_split(:,2),'T)');
    end
    datainfo.field = field_split(:,1);

    % Power
    if settings.vti
        power_split = split(split_file_names(:,5),'=');
        power_split = split(power_split(:,2),'dBm');
    elseif settings.dilutionfridge
        power_split = split(split_file_names(:,7),'dBm.dat');
        power_split_string = strings([length(power_split),1]);
        for l = 1:length(power_split)
            if contains(power_split(l),'min')
                power_split_string(l) = replace(string(power_split(l)),'min','-'); % negative values
            else
                power_split_string(l) = power_split(l); % positive values
            end
        end    
        power_split = cellstr(power_split_string);
    end
    datainfo.power = power_split(:,1);

    % Actualtemp
    if settings.useactualtemp
        if settings.vti
            actualtemp_split = split(split_file_names(:,6),'=');
            actualtemp_split = split(actualtemp_split(:,2),'K');
        elseif settings.dilutionfridge
            actualtemp_split = split(split_file_names(:,4),'=');
            actualtemp_split = split(actualtemp_split(:,2),'K');
        end
        datainfo.acttemp = actualtemp_split(:,1);
    end
end

% Select modes/spectra and multimodes
function loopinfo = select_modes(settings,datainfo)
    loopinfo.name = unique(datainfo.name,'stable'); % Unique modes/spectra

    % Select modes/spectra to sweep
    if length(loopinfo.name) > 1
        [usenameindex,tf] = listdlg('PromptString',{'Select all modes you want to fit.'},'ListString',loopinfo.name);
        if ~tf % Exit if nothing is selected
            disp('No mode selected.')
            return
        end
        loopinfo.usename = loopinfo.name(usenameindex);
    else
        loopinfo.usename = loopinfo.name;
    end

    % Select Multimodes
    if settings.multimodes 
        if length(loopinfo.usename) > 1 
            [usemultiindex,tf] = listdlg('PromptString',{'Select all modes with multiple modes.',''},'ListString',loopinfo.usename);
            if ~tf % exit
                disp('No mode selected.')
                return
            end
            loopinfo.multimodes = loopinfo.usename(usemultiindex);
        else
            loopinfo.multimodes = loopinfo.usename;
        end
        for z = 1:length(loopinfo.multimodes) % Loop selected multimodes
            for y = 1:length(datainfo.name) % Plot spectra and ask for amount of modes
                if isequal(loopinfo.multimodes(z),datainfo.name(y))
                    selectfile = importdata(fullfile(datainfo.datapath, char(datainfo.filename(y))));
                    figure(3)
                        plot(selectfile.data(:,1),selectfile.data(:,4)),title(['Spectra for ' char(loopinfo.multimodes(z)) ':'])
                        answer = str2double(inputdlg('How many modes do you want to fit within this spectra?','Multimodes'));
                    close all
                    for w = 1:answer-1
                        index = find(ismember(loopinfo.usename,loopinfo.multimodes(z)),1);
                        loopinfo.usename(index+1:end+1) = loopinfo.usename(index:end);
                        loopinfo.usename(index) = loopinfo.multimodes(z);
                    end
                    loopinfo.multimodenumber(z) = answer;
                    break
                end
            end
            loopinfo.multimodecounter(z) = 0;
        end
    end
end
