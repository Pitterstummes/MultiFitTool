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
settings.fitdelay           = 0;        % Time delay inbetween fits
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
dowork(settings,datainfo,loopinfo);

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

% Loop modes, do fitting and save
function dowork(settings,datainfo,loopinfo)
    % Loop modes
    for usename = loopinfo.usename
        % Initialize multimodes
        if settings.multimodes && ismember(usename,loopinfo.multimodes) && ~settings.oldmode
            multimodeindex = ismember(loopinfo.multimodes,usename);
            loopinfo.multimodecounter(multimodeindex) = loopinfo.multimodecounter(multimodeindex) + 1;
            currentmultimodemax = loopinfo.multimodenumber(multimodeindex);
            currentmultimodestep = loopinfo.multimodecounter(multimodeindex);
        end 

        % Select relevant data & sweeps
        if ~settings.oldmode
            loopinfo = select_data(settings,datainfo,loopinfo,usename,currentmultimodestep,currentmultimodemax);
        end

        % Loop over sweeps
        for usesweeper = loopinfo.usesweeper
            % Loop over fixed parameters for each sweep
            for sweeper = loopinfo.sweeper
                if ~settings.oldmode
                    % Specify parameters that are fixed while sweeping
                    loopinfo = specify_sweeps(settings,loopinfo,usename,usesweeper,sweeper,currentmultimodestep,currentmultimodemax);

                    % Set values for fitting
                    if isequal(usesweeper,{'temp'}) 
                        fixedtemp  = 1;
                        fixedfield = length(loopinfo.usefield);
                        fixedpower = length(loopinfo.usepower);
                    elseif isequal(usesweeper,{'field'})
                        fixedtemp  = length(loopinfo.usetemp);
                        fixedfield = 1;
                        fixedpower = length(loopinfo.usepower);
                    elseif isequal(usesweeper,{'power'})
                        fixedtemp  = length(loopinfo.usetemp);
                        fixedfield = length(loopinfo.usefield);
                        fixedpower = 1;
                    end
                end

                % Sweep trough selected cases and get related files
                counter = 0; 
                for usefixedtemp = 1:fixedtemp
                    for usefixedfield = 1:fixedfield
                        for usefixedpower = 1:fixedpower
                            [loopinfo,shortsweep,shortsweepunit,counter,maxcounter,currentfiles] = more_stuff(settings,loopinfo,...
                                usename,usesweeper,counter,usefixedtemp,usefixedfield,usefixedpower,fixedtemp,fixedfield,fixedpower,currentmultimodestep);

                            % Maxmultifit
                            if settings.multifit >= length(currentfiles)
                                settings.multifit = length(currentfiles);
                            end
                            file.end = importdata([datainfo.datapath,'\',char(currentfiles(end))]);

                            % Init for the final loop trough files & fitmenu
                            i = 1;
                            settings.usemultifit     = settings.multifit;       % use*values since they could be changed via fitmenu
                            settings.usemultifitname = settings.multifitname;
                            didthisalready           = [false,false];           % workrange, breakmenu
                            breakswitch              = 1;

                            while i <= length(currentfiles)
                                file.file = importdata([datainfo.datapath,'\',char(currentfiles(i))]);
                                if i == 1 && ~didthisalready(1)
                                    didthisalready(1) = true;
                                    if ~settings.oldmode
                                        % get sweeper value array for analysisfile
                                        if isequal(usesweeper,{'temp'}) 
                                            loopinfo.analysissweeper = loopinfo.usetemp;
                                        elseif isequal(usesweeper,{'field'})
                                            loopinfo.analysissweeper = loopinfo.usefield;
                                        elseif isequal(usesweeper,{'power'})
                                            loopinfo.analysissweeper = loopinfo.usepower;
                                        end
                                    end
                                    % pick working range
                                    if settings.repickworkingrange || ~settings.repickworkingrange && counter == 1 || settings.oldmode
                                        figure(1)
                                            plot(file.file.data(:,1),file.file.data(:,4),file.end.data(:,1),file.end.data(:,4)), ...
                                                if ~settings.oldmode
                                                    if settings.multimodes && ismember(usename,loopinfo.multimodes)
                                                        title({...
                                                        [char(usename) ', ' char(usesweeper) 'sweep (' num2str(counter) '/' num2str(maxcounter) '), multimode (' num2str(currentmultimodestep) '/' num2str(currentmultimodemax) '):' ] ...
                                                        'pick working range'})
                                                    else
                                                        title({ ...
                                                        [char(usename) ', ' char(usesweeper) 'sweep (' num2str(counter) '/' num2str(maxcounter) '):'] ...
                                                        'pick working range'})
                                                    end
                                                else
                                                    title('pick working range')
                                                end
                                            [lbgrob,ubgrob,~] = getbounds(file.file.data);    
                                    end
                                end
                                if i == 1 
                                    % phase fit range
                                    figure(1)
                                        plot(file.file.data(lbgrob:ubgrob,1),file.file.data(lbgrob:ubgrob,2),file.file.data(lbgrob:ubgrob,1),file.file.data(lbgrob:ubgrob,3),file.file.data(lbgrob:ubgrob,1),file.file.data(lbgrob:ubgrob,4)),...
                                            if ~settings.oldmode
                                                if settings.multimodes && ismember(usename,loopinfo.multimodes)
                                                    title({...
                                                    [char(usename) ', ' char(usesweeper) 'sweep (' num2str(counter) '/' num2str(maxcounter) '), multimode (' num2str(currentmultimodestep) '/' num2str(currentmultimodemax) '):' ] ...
                                                    'pick range for phase fit'})
                                                else
                                                    title({ ...
                                                    [char(usename) ', ' char(usesweeper) 'sweep (' num2str(counter) '/' num2str(maxcounter) '):'] ...
                                                    'pick range for phase fit'})
                                                end
                                            else
                                                title('pick range for phase fit')
                                            end
                                        [lbtau,ubtau,~] = getbounds(file.file.data);
                                    % phase fit
                                    linfun = @(v,xdata) v(1)*xdata+v(2);
                                    best   = lsqcurvefit(linfun,[-193,0],file.file.data(lbtau:ubtau,1)./1E9,unwrap(angle(complex(file.file.data(lbtau:ubtau,2),file.file.data(lbtau:ubtau,3)))));
                                    f2     = figure(2);
                                    plot(file.file.data(lbtau:ubtau,1)./1E9,unwrap(angle(complex(file.file.data(lbtau:ubtau,2),file.file.data(lbtau:ubtau,3)))),file.file.data(lbtau:ubtau,1)./1E9,linfun(best,file.file.data(lbtau:ubtau,1)./1E9)),title('this is your phase fit');
                                    tau    = -best(1);
                                end
                                % load, plot and pick plotintervall, fitmenu
                                if mod(i-1,settings.usemultifit) == 0 || multifitrepick
                                    multifitrepick = false;
                                    if i-1+settings.usemultifit < length(currentfiles)
                                        file.step=importdata([datainfo.datapath,'\',char(currentfiles(i-1+settings.usemultifit))]);
                                    else
                                        file.step = file.end;
                                    end
                                    figure(1);
                                        plot(file.file.data(lbgrob:ubgrob,1),file.file.data(lbgrob:ubgrob,4),file.step.data(lbgrob:ubgrob,1),file.step.data(lbgrob:ubgrob,4)),...
                                            if settings.usemultifit == 1
                                                if ~settings.oldmode
                                                    title(['Pick range for ' num2str(settings.usemultifit) ' ' settings.usemultifitname ', ' shortsweep ' = ' sprintf('%g',str2double(loopinfo.analysissweeper(i))) shortsweepunit ])
                                                else
                                                    title(['Pick range for ' num2str(settings.usemultifit) ' ' settings.usemultifitname ', sweepparameter = ' sprintf('%g',str2double(loopinfo.analysissweeper(i))) ])
                                                end
                                            else
                                                if i-1+settings.usemultifit <= length(currentfiles)
                                                    if ~settings.oldmode
                                                        title(['Pick range for ' num2str(settings.usemultifit) ' ' settings.usemultifitname ', ' shortsweep ' = ' sprintf('%g',str2double(loopinfo.analysissweeper(i))) '-' sprintf('%g',str2double(loopinfo.analysissweeper(i-1+settings.usemultifit))) shortsweepunit ])
                                                    else
                                                        title(['Pick range for ' num2str(settings.usemultifit) ' ' settings.usemultifitname ', sweepparameter = ' sprintf('%g',str2double(loopinfo.analysissweeper(i))) '-' sprintf('%g',str2double(loopinfo.analysissweeper(i-1+settings.usemultifit))) ])
                                                    end
                                                else
                                                    if ~settings.oldmode
                                                        title(['Pick range for ' num2str(length(currentfiles) - i + 1) ' ' settings.usemultifitname ', ' shortsweep ' = ' sprintf('%g',str2double(loopinfo.analysissweeper(i))) '-' sprintf('%g',str2double(loopinfo.analysissweeper(length(currentfiles)))) shortsweepunit ])
                                                    else
                                                        title(['Pick range for ' num2str(length(currentfiles) - i + 1) ' ' settings.usemultifitname ', sweeperparameter = ' sprintf('%g',str2double(loopinfo.analysissweeper(i))) '-' sprintf('%g',str2double(loopinfo.analysissweeper(length(currentfiles)))) ])
                                                    end
                                                end
                                            end                         
                                        [lb,ub,x] = getbounds(file.file.data);
                                        if x(1) == x(2) % open fitting menu via picking the same point 2 times
                                            fitmenuanswer = questdlg({'What would you like to do?','Exit menu by pressing X'}, 'Fitting Menu', 'Repick phase', 'End fit sequence', 'Change multifits', 'Change multifits');
                                            % Handle response
                                            switch fitmenuanswer
                                                case 'End fit sequence'
                                                    disp('breaking the current loop')
                                                    disp('fitting finished')
                                                    close all
                                                    break
                                                case 'Change multifits'
                                                    newmultifit = inputdlg('Type in new multifit number.','New Multifit');
                                                    if ~isempty(newmultifit)
                                                        if isa(str2double(newmultifit),'double') && isfinite(str2double(newmultifit)) && str2double(newmultifit) > 0 && mod(str2double(newmultifit),1) == 0 %check if natural
                                                            settings.usemultifit = str2double(newmultifit); 
                                                            if settings.usemultifit == 1
                                                                settings.usemultifitname = 'fit';
                                                            else
                                                                settings.usemultifitname = 'fits';
                                                            end
                                                        end
                                                    end
                                                    multifitrepick = true;
                                                    continue
                                                case 'Repick phase'
                                                    i = 1;
                                                    continue
                                            end
                                            if isempty(fitmenuanswer)
                                                continue
                                            end
                                        end
                                        if i == 1 && settings.repickworkingrange % close phase fit fig
                                            close(f2)
                                        end
                                end
                                % fit cases according to the normal multifittool, still optimizable
                                if i == 1 % once again, do this only one time
                                    % analysis file header
                                    analysis_ID = fopen([datainfo.analysispath '\' loopinfo.analysisfile],'w');
                                    if ~settings.oldmode
                                        fprintf(analysis_ID,[shortsweep '\t f0 \t fB \t Q \r\n' shortsweepunit '\t GHz \t GHz \t \r\n']);
                                    elseif settings.oldmode
                                        fprintf(analysis_ID,'X \t f0 \t fB \t Q \r\n \t GHz \t GHz \t \r\n');
                                    end
                                    fclose(analysis_ID);
                                    % get initial fit parameters
                                    offset_init         = (mean(file.file.data(lb:lb+9,4))+mean(file.file.data(ub-9:ub,4)))/2;
                                    f_v                 = file.file.data(lb:ub,1).*1e-9;
                                    [bandh_init,f0_ind] = max(abs(file.file.data(lb:ub,4)-offset_init));
                                    [~,fB_ind]          = min(abs(file.file.data(lb:ub,4)-bandh_init/2));
                                    f0_init             = f_v(f0_ind);
                                    fB_init             = max([abs(f_v(fB_ind)-f0_init)*2;0.001]);
                                    initial_guess       = [tau/(2*pi);bandh_init*fB_init*1i;f0_init+1i*fB_init/2;0;offset_init];
                                    % fit with initial_guess parameters
                                    fitparameter = fit_complex(file.file.data,initial_guess,lb,ub,0,1,tau);
                                    fB           = imag(fitparameter(3))*2;
                                    f0           = real(fitparameter(3));
                                    Q            = f0/fB;                            
                                    tau = fitparameter(1);
                                    fitparameter = fitparameter(2:end);
                                else % fit with optained fit parameters, keep tau fixed.
                                    fitparameter = fit_complex(file.file.data,fitparameter,lb,ub,settings.fastmode,0,tau);
                                    fB           = imag(fitparameter(2))*2;
                                    f0           = real(fitparameter(2));
                                    Q            = f0/fB;
                                end
                                analysis_ID = fopen([datainfo.analysispath '\' loopinfo.analysisfile],'a'); % write/save data
                                    fprintf(analysis_ID,'%f \t %.9e \t %.9e \t %.9e \r\n',...
                                    str2double(loopinfo.analysissweeper(i)),f0,fB,Q);
                                fclose(analysis_ID);
                                % Break menu: open by typing 'b' while fitprocess
                                currentchar1 = get(gcf, 'CurrentCharacter'); % get last pressed Character (keyboard)
                                if currentchar1 == 'b'
                                    didthisalready(2) = true;
                                end
                                if didthisalready(2) && breakswitch == 1
                                    breakswitch = 2;
                                    breakmenuanswer = questdlg({'What would you like to do?','Exit menu by pressing X or Esc.','If this opens unwanted, pick (confirm) ranges with "n".'}, 'Break Menu', 'Repick workingrange', 'End fit sequence', 'Repick phase', 'Repick phase');
                                    if isempty(breakmenuanswer)
                                        i = 1;
                                        continue
                                    end
                                    switch breakmenuanswer
                                        case 'End fit sequence'
                                            disp('breaking the current loop')
                                            disp('fitting finished')
                                            close all
                                            break
                                        case 'Repick workingrange'
                                            i = 1;
                                            didthisalready(1) = false;
                                            continue
                                        case 'Repick phase'
                                            i = 1;
                                            continue
                                    end
                                elseif breakswitch == 2
                                    breakswitch = 1;
                                end
                                didthisalready(2) = false;
                                if i == length(currentfiles) % current fit succession ended
                                    disp('fitting finished')
                                    pause(1)
                                    close all
                                end 
                                i = i + 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

% Select sweeps
function loopinfo = select_data(settings,datainfo,loopinfo,usename,currentmultimodestep,currentmultimodemax)
    keep = false(length(datainfo.name),1);
    for i = 1:length(datainfo.name)
        keep(i) = isequal(datainfo.name(i),usename);
    end
    loopinfo.usefile    = datainfo.filename(keep);
    loopinfo.temp       = unique(datainfo.temp(keep),'stable');
    loopinfo.field      = unique(datainfo.field(keep),'stable');
    loopinfo.power      = unique(datainfo.power(keep),'stable');

    % Initialize Sweepcases
    loopinfo.sweepername = {'temp', 'field', 'power'}'; 
    loopinfo.sweepertest = [length(loopinfo.temp) length(loopinfo.field) length(loopinfo.power)];
    loopinfo.sweeper = {};
    for i = 1:3 % generate possible sweeper list
        if loopinfo.sweepertest(i) > 1 % more than 1 value -> sweeper
            loopinfo.sweeper(length(loopinfo.sweeper)+1) = loopinfo.sweepername(i); % array of used sweepernames
        end
    end
    loopinfo.usesweeper = loopinfo.sweeper;
    loopinfo.usetemp    = loopinfo.temp;
    loopinfo.usefield   = loopinfo.field;
    loopinfo.usepower   = loopinfo.power;

    % Two or three sweeps
    if length(loopinfo.sweeper) >= 1
        if settings.multimodes && ismember(usename,loopinfo.multimodes)
            [sweeperindex,tf] = listdlg('PromptString',{['Regarding ' char(usename) ','],... % select sweeps
            ['multimode (' num2str(currentmultimodestep) '/' num2str(currentmultimodemax) '):'],...
            'What would you like to sweep?'},'ListString',loopinfo.sweeper);
        else
            [sweeperindex,tf] = listdlg('PromptString',{['Regarding ' char(usename) ':'],... % select sweeps
            'What would you like to sweep?'},'ListString',loopinfo.sweeper);
        end
        if ~tf % exit
            return
        end
        loopinfo.usesweeper = loopinfo.sweeper(sweeperindex);
    end
end

% Specify sweeps
function loopinfo = specify_sweeps(settings,loopinfo,usename,usesweeper,sweeper,currentmultimodestep,currentmultimodemax)
    if ~isequal(sweeper,usesweeper) % if there is something to choose, determine the case
        if isequal(sweeper,{'temp'})
            parameterarray  = loopinfo.temp;
            paramcase       = 1;
        elseif isequal(sweeper,{'field'})
            parameterarray  = loopinfo.field;
            paramcase       = 2;
        elseif isequal(sweeper,{'power'})
            parameterarray  = loopinfo.power;
            paramcase       = 3;
        end
        % choose fixed parameters
        if settings.multimodes && ismember(usename,loopinfo.multimodes)
            [fixparameterindex,tf] = listdlg('PromptString',{['At which ' char(sweeper) ...
            ' do you want your ' char(usesweeper) 'sweep for ' char(usename) ','...
            'multimode (' num2str(currentmultimodestep) '/' num2str(currentmultimodemax) ')?'],'',''},'ListString',parameterarray);
        else
            [fixparameterindex,tf] = listdlg('PromptString',{['At which ' char(sweeper) ...
            ' do you want your ' char(usesweeper) 'sweep for ' char(usename) '?'],''},'ListString',parameterarray);
        end
        if ~tf % exit
            return
        end
        % fill in the correct useparameter (since not done yet)
        if paramcase == 1
            loopinfo.usetemp  = loopinfo.temp(fixparameterindex);
        elseif paramcase == 2
            loopinfo.usefield = loopinfo.field(fixparameterindex);
        elseif paramcase == 3
            loopinfo.usepower = loopinfo.power(fixparameterindex);
        end                   
    else % trivial case, maybe allready done before?
        if isequal(sweeper,{'temp'})
            loopinfo.usetemp  = loopinfo.temp;
        elseif isequal(sweeper,{'field'})
            loopinfo.usefield = loopinfo.field;
        elseif isequal(sweeper,{'power'})
            loopinfo.usepower = loopinfo.power;
        end
    end
end

% More stuff LOL
function [loopinfo,shortsweep,shortsweepunit,counter,maxcounter,currentfiles] = more_stuff(settings,loopinfo,usename,usesweeper,...
    counter,usefixedtemp,usefixedfield,usefixedpower,fixedtemp,fixedfield,fixedpower,currentmultimodestep)
    if ~settings.oldmode
        currentfilesindex = false(length(loopinfo.usefile),1);
        for q = 1:length(loopinfo.usefile) % get index of all needed files for certain cases + analysis file
            if isequal(usesweeper,{'temp'})
                if settings.vti
                    currentfilesindex(q) = contains(loopinfo.usefile(q),['K_Field=' char(loopinfo.usefield(usefixedfield))]) && contains(loopinfo.usefile(q),['T_Power=' char(loopinfo.usepower(usefixedpower))]);
                elseif settings.dilutionfridge
                    if contains(loopinfo.usepower(usefixedpower),'-')
                        usepower_minfix = replace(string(loopinfo.usepower(usefixedpower)),'-','min');
                    end
                    currentfilesindex(q) = contains(loopinfo.usefile(q),['K_B1(setB=' char(loopinfo.usefield(usefixedfield))]) && contains(loopinfo.usefile(q),[char(usepower_minfix) 'dBm.dat']);
                end
                if q == length(loopinfo.usefile)
                    if settings.multimodes && ismember(usename,loopinfo.multimodes)
                        loopinfo.analysisfile = ['analysis_' char(usename) '-mode' num2str(currentmultimodestep) '_tempsweepstep=' num2str(str2double(char(loopinfo.usetemp(2)))-str2double(char(loopinfo.usetemp(1)))) 'K_temp=' sprintf('%g',str2double(loopinfo.usetemp(1))) '-' sprintf('%g',str2double(loopinfo.usetemp(end))) 'K_field=' sprintf('%g',str2double(loopinfo.usefield(usefixedfield))) 'T_power=' sprintf('%g',str2double(loopinfo.usepower(usefixedpower))) 'dBm.dat'];
                    else
                        loopinfo.analysisfile = ['analysis_' char(usename) '_tempsweepstep=' num2str(str2double(char(loopinfo.usetemp(2)))-str2double(char(loopinfo.usetemp(1)))) 'K_temp=' sprintf('%g',str2double(loopinfo.usetemp(1))) '-' sprintf('%g',str2double(loopinfo.usetemp(end))) 'K_field=' sprintf('%g',str2double(loopinfo.usefield(usefixedfield))) 'T_power=' sprintf('%g',str2double(loopinfo.usepower(usefixedpower))) 'dBm.dat'];
                    end
                end
                shortsweep     = 'T';
                shortsweepunit = 'K';
            elseif isequal(usesweeper,{'field'})
                if settings.vti
                    currentfilesindex(q) = contains(loopinfo.usefile(q),['_setTemp=' char(loopinfo.usetemp(usefixedtemp))]) && contains(loopinfo.usefile(q),['T_Power=' char(loopinfo.usepower(usefixedpower))]);
                elseif settings.dilutionfridge
                    if contains(loopinfo.usepower(usefixedpower),'-')
                        usepower_minfix = replace(string(loopinfo.usepower(usefixedpower)),'-','min');
                    end
                    currentfilesindex(q) = contains(loopinfo.usefile(q),['_setTmc=' char(loopinfo.usetemp(usefixedtemp))]) && contains(loopinfo.usefile(q),[char(usepower_minfix) 'dBm.dat']);
                end
                if q == length(loopinfo.usefile)
                    if settings.multimodes && ismember(usename,loopinfo.multimodes)
                        loopinfo.analysisfile = ['analysis_' char(usename) '-mode' num2str(currentmultimodestep) '_fieldsweepstep=' num2str(str2double(char(loopinfo.usefield(2)))-str2double(char(loopinfo.usefield(1))))  'T_temp=' sprintf('%g',str2double(loopinfo.usetemp(usefixedtemp))) 'K_field=' sprintf('%g',str2double(loopinfo.usefield(1))) '-' sprintf('%g',str2double(loopinfo.usefield(end))) 'T_power=' sprintf('%g',str2double(loopinfo.usepower(usefixedpower))) 'dBm.dat'];
                    else
                        loopinfo.analysisfile = ['analysis_' char(usename) '_fieldsweepstep=' num2str(str2double(char(loopinfo.usefield(2)))-str2double(char(loopinfo.usefield(1))))  'T_temp=' sprintf('%g',str2double(loopinfo.usetemp(usefixedtemp))) 'K_field=' sprintf('%g',str2double(loopinfo.usefield(1))) '-' sprintf('%g',str2double(loopinfo.usefield(end))) 'T_power=' sprintf('%g',str2double(loopinfo.usepower(usefixedpower))) 'dBm.dat'];
                    end
                end
                shortsweep     = 'B';
                shortsweepunit = 'T';
            elseif isequal(usesweeper,{'power'})
                if settings.vti
                    currentfilesindex(q) = contains(loopinfo.usefile(q),['_setTemp=' char(loopinfo.usetemp(usefixedtemp))]) && contains(loopinfo.usefile(q),['K_Field=' char(loopinfo.usefield(usefixedfield))]);
                elseif settings.dilutionfridge
                    currentfilesindex(q) = contains(loopinfo.usefile(q),['_setTmc=' char(loopinfo.usetemp(usefixedtemp))]) && contains(loopinfo.usefile(q),['K_B1(setB=' char(loopinfo.usefield(usefixedfield))]);
                end
                if q == length(loopinfo.usefile)
                    if settings.multimodes && ismember(usename,loopinfo.multimodes)
                        loopinfo.analysisfile = ['analysis_' char(usename) '-mode' num2str(currentmultimodestep) '_powersweepstep=' num2str(str2double(char(loopinfo.usepower(2)))-str2double(char(loopinfo.usepower(1))))  'dBm_temp=' sprintf('%g',str2double(loopinfo.usetemp(usefixedtemp))) 'K_field=' sprintf('%g',str2double(loopinfo.usefield(usefixedfield))) 'T_power=' sprintf('%g',str2double(loopinfo.usepower(1))) '-' sprintf('%g',str2double(loopinfo.usepower(end))) 'dBm.dat'];
                    else
                        loopinfo.analysisfile = ['analysis_' char(usename) '_powersweepstep=' num2str(str2double(char(loopinfo.usepower(2)))-str2double(char(loopinfo.usepower(1))))  'dBm_temp=' sprintf('%g',str2double(loopinfo.usetemp(usefixedtemp))) 'K_field=' sprintf('%g',str2double(loopinfo.usefield(usefixedfield))) 'T_power=' sprintf('%g',str2double(loopinfo.usepower(1))) '-' sprintf('%g',str2double(loopinfo.usepower(end))) 'dBm.dat'];
                    end
                end
                shortsweep     = 'P';
                shortsweepunit = 'dBm';
            end
        end
        currentfiles = loopinfo.usefile(currentfilesindex); 
        if ~settings.oldmode || ~settings.repickworkingrange
            counter = counter + 1; 
            if counter == 1
                maxcounter = fixedtemp*fixedfield*fixedpower;
            end
        end
    end
end

% Pick range
function [lb,ub,x] = getbounds(data)
    [x,~] = ginput(2);
    if x(1) < x(2)
        lb=find(data(:,1)>x(1),1,'first');
        ub=find(data(:,1)<x(2),1,'last');
    else
        lb=find(data(:,1)>x(2),1,'first');
        ub=find(data(:,1)<x(1),1,'last');
    end
end

% Fitfunction with fixed tau after first fit
function fitparameter = fit_complex(data,initial_guess,lb,ub,fastmode,initial_run,tau) 
    fitfreq = data(lb:ub,1).*1e-9;
    complexdata = complex(data(lb:ub,2),data(lb:ub,3));
    if initial_run
        fitfunction = @(p,f)exp(-p(1).*2.*pi.*1i.*f).*(p(2)./(f-p(3))+p(4)+p(5)*(f-real(p(3))));
    else
        fitfunction = @(p,f)exp(-tau.*2.*pi.*1i.*f).*(p(1)./(f-p(2))+p(3)+p(4)*(f-real(p(2))));
    end
    options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','MaxFunEvals',100000);
    fitparameter = lsqcurvefit(fitfunction,initial_guess,fitfreq,complexdata,[],[],options);
    fitcurve = fitfunction(fitparameter,fitfreq);
    if ~fastmode
        figure(2)
        plot(fitfreq,real(complexdata),fitfreq,imag(complexdata),fitfreq,real(fitcurve),fitfreq,imag(fitcurve))
        drawnow
    end
    figure(1)
    plot(fitfreq,abs(complexdata).^2,fitfreq,abs(fitcurve).^2)
    drawnow
end
