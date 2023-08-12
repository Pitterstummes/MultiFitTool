clc         % clear command window
close all   % closes all figures
clearvars   % remove all variables from current active workspace

% settings (false = 0, true = 1, both works)
settings.startfolder        = 'C'; % input root folder
settings.oldmode            = 0;    % just select datafolder, sweepfile and fit it. tells you when to activate
settings.debug              = 0;    % dont clear any variables
settings.fastmode           = 1;    % double fitting speed, may be too fast! (no more slowing down if mouse is on figure)
settings.flipdirection      = 0;    % sweep through files in flipped direction
settings.autoanalysisfolder = 1;    % automatically creats an analysis folder
settings.multimodes         = 0;    % set true if you want to fit multiple modes within a file
settings.repickworkingrange = 0;    % select new workingrange after finished fitting round
settings.multifit           = 70;   % fit X files before selecting a new fit range
% fixed
settings.useactualtemp  = false;    % not used for anything so far
settings.useheliumlevel = false;    % not used for anything so far
settings.programpath    = pwd;      % current folder has to contain the program. otherwise add to PATH, do not change

% multifit: initialize right names
if settings.multifit == 1
    settings.multifitname = 'fit';
else
    settings.multifitname = 'fits';
end

% get data path & path above
datainfo.datapath=uigetdir(settings.startfolder,'Select data folder');
if datainfo.datapath == 0
    disp('Canceled data folder selection.')
    return
end
cd(datainfo.datapath)
cd ..
datainfo.previousdatapath = pwd;

% analysis folder
if ~settings.autoanalysisfolder
    datainfo.analysispath=uigetdir(datainfo.previousdatapath,'Select analysis folder');
    cd(settings.programpath)
    if datainfo.analysispath == 0
        disp('Canceled analysis folder selection.')
        return
    end
else
    cd(datainfo.previousdatapath)
    if ~isfolder('Analysis')
      mkdir Analysis;    
    end
    datainfo.analysispath = strcat(pwd,'\Analysis');
    cd(settings.programpath)
end

% get file names
datainfo.fileinfo = dir([datainfo.datapath '\' '*.dat']); %(alphabetical sorted)
datainfo.filename = {datainfo.fileinfo.name}';

% setup detection, add more cases!
% vti
settings.vti = contains(datainfo.filename(1),'_setTemp=') && contains(datainfo.filename(1),'K_Field=') ...
    && contains(datainfo.filename(1),'T_Power=') && contains(datainfo.filename(1),'_actTemp=') ...
    && contains(datainfo.filename(1),'percent.dat');
% dilution refridgerator
settings.dilutionfridge = contains(datainfo.filename(1),'_setTmc=') && contains(datainfo.filename(1),'K_actTsample=')...
    && contains(datainfo.filename(1),'K_B') && contains(datainfo.filename(1),'(setB=') && contains(datainfo.filename(1),'dBm.dat');
% case: errorexit
if ~settings.vti && ~settings.dilutionfridge && ~settings.oldmode
    f = errordlg({'Tool is not yet adapted for the used data.','Please use the settings.oldmode.'},'Data Error');
    return
end

% sort filenames numerically, flip if selected
try %check if "_" was used accidently 
    samedelimiter_test = split(datainfo.filename,'_');
catch
    f = errordlg({'Apparently a "_" was used incorrectly in a filename.','Please replace all "_" that are not used as delimiters with "-".',...
        'You can use the program "fix_underscore_in_oldfridge_filenames.m" in "\\S4\Datenpool\Microwave Group\Programs".'},'Filename Error');
end
split_file_names  = split(datainfo.filename,'_');       % split in parts (alphabetical sorted)
doublenumbers     = str2double(split_file_names(:,1));  % type to double
[~, sortorder]    = sort(doublenumbers);                % get right sortorder
datainfo.filename = datainfo.filename(sortorder);       % sorted file names
if ~settings.debug
    clear doublenumbers sortorder
end
if settings.flipdirection % flip direction if true
    datainfo.filename = flip(datainfo.filename);
end

if ~settings.oldmode % skip when oldmode is active
    split_file_names = split(datainfo.filename,'_'); %split in parts (numerical sorted)

    % extract all data from sorted filename -> create datainfo struct
    % number
    datainfo.number = split_file_names(:,1);
    % name
    if settings.vti
        datainfo.name = split_file_names(:,2);
    elseif settings.dilutionfridge
        % name_split = split(split_file_names(:,6),'(');
        % datainfo.name = name_split(:,1); 
        datainfo.name = split_file_names(:,6); % also adds the center frequency to the name
    end
    % settemp
    if settings.vti || settings.dilutionfridge
        settemp_split = split(split_file_names(:,3),'=');
        settemp_split = split(settemp_split(:,2),'K');
    end
    datainfo.temp = settemp_split(:,1);
    % field
    if settings.vti
        field_split = split(split_file_names(:,4),'=');
        field_split = split(field_split(:,2),'T');
    elseif settings.dilutionfridge
        field_split = split(split_file_names(:,5),'=');
        field_split = split(field_split(:,2),'T)');
    end
    datainfo.field = field_split(:,1);
    % datainfo.field = field_split(:,1);
    % power (add case of positive powers at the dilution fridge?)
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
    % actualtemp
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
    % heliumlevel
    if settings.vti && settings.useheliumlevel
        heliumlevel_split = split(split_file_names(:,7),'percent.dat');
        datainfo.helium   = heliumlevel_split(:,1);
    end
    if ~settings.debug % clear no longer used variables
        clear split_file_names settemp_split power_split heliumlevel_split field_split actualtemp_split l name_split power_split_string
    end
    
    loopinfo.name = sort(unique(datainfo.name));
    % select modes (names) to sweep
    if length(loopinfo.name) > 1
        [usenameindex,tf] = listdlg('PromptString',{'Select all modes you want to fit.'},'ListString',loopinfo.name);
        if ~tf % exit
            return
        end
        loopinfo.usename = loopinfo.name(usenameindex);
        if ~settings.debug
            clear usenameindex tf
        end
    else
        loopinfo.usename = loopinfo.name;
    end
    
    % multimodes (select modenames(spectras) with multimodes, init, add multimodes to usename
    if settings.multimodes 
        if length(loopinfo.usename) > 1 % select  multimodes
            [usemultiindex,tf] = listdlg('PromptString',{'Select all modes with multiple modes.',''},'ListString',loopinfo.usename);
            if ~tf % exit
                return
            end
            loopinfo.multimodes = loopinfo.usename(usemultiindex);
        else
            loopinfo.multimodes = loopinfo.usename;
        end
        for z = 1:length(loopinfo.multimodes) % loop selected multimodes
            for y = 1:length(datainfo.name) % plot spectra and ask for multimodenumber
                if isequal(loopinfo.multimodes(z),datainfo.name(y))
                    file.choose = importdata([datainfo.datapath,'\',char(datainfo.filename(y))]);
                    figure(3)
                        plot(file.choose.data(:,1),file.choose.data(:,4)),title(['Spectra for ' char(loopinfo.multimodes(z)) ':'])
                        answer = str2double(inputdlg('How many modes do you want to fit within this spectra?','Multimodes'));
                    close all
                    for w = 1:answer-1
                        loopinfo.usename(end+1) = loopinfo.multimodes(z);
                    end
                    loopinfo.multimodenumber(z) = answer;
                    break
                end
            end
            loopinfo.multimodecounter(z) = 0;
        end
        loopinfo.usename = sort(loopinfo.usename);
    end
else % init parameters for oldmode to skip loops,..
    loopinfo.usename    = 1;
    loopinfo.usesweeper = 1;
    loopinfo.sweeper    = 1;
    fixedtemp           = 1;
    fixedfield          = 1;
    fixedpower          = 1;
    % get files, sweeperfile and analysisfilename
    currentfiles = datainfo.filename;
    [datainfo.sweeperfile,datainfo.sweeperpath] = uigetfile('*.*','Choose sweeper file',datainfo.previousdatapath);
    loopinfo.analysissweeper = importdata([datainfo.sweeperpath '\' datainfo.sweeperfile]);
    loopinfo.analysissweeper = num2str(loopinfo.analysissweeper);
    loopinfo.analysissweeper = cellstr(loopinfo.analysissweeper);
    if settings.flipdirection
        loopinfo.analysissweeper = flip(loopinfo.analysissweeper);
    end
    datainfo.datafoldername = split(datainfo.datapath,'\');
    datainfo.datafoldername = datainfo.datafoldername(end);
    loopinfo.analysisfile   = ['analysis_XXX_' char(datainfo.datafoldername) '_step=' num2str(str2double(char(loopinfo.analysissweeper(2)))-str2double(char(loopinfo.analysissweeper(1))))  '_range=' sprintf('%g',str2double(loopinfo.analysissweeper(1))) '-' sprintf('%g',str2double(loopinfo.analysissweeper(end))) '.dat'];
end

% work
for usenameind = 1:length(loopinfo.usename) % loop over selected modes
    usename = loopinfo.usename(usenameind);
    % multimodes init
    if settings.multimodes && ismember(usename,loopinfo.multimodes) && ~settings.oldmode
        multimodeindex = ismember(loopinfo.multimodes,usename);
        loopinfo.multimodecounter(multimodeindex) = loopinfo.multimodecounter(multimodeindex) + 1;
        currentmultimodemax = loopinfo.multimodenumber(multimodeindex);
        currentmultimodestep = loopinfo.multimodecounter(multimodeindex);
    end   
    if ~settings.oldmode % once again, skip for oldmode
        % select mode relevant data
        keep = false(length(datainfo.name),1);
        for i = 1:length(datainfo.name)
            keep(i) = isequal(datainfo.name(i),usename);
        end
        loopinfo.usefile    = datainfo.filename(keep);
        loopinfo.temp       = unique(datainfo.temp(keep),'stable');
        loopinfo.field      = unique(datainfo.field(keep),'stable');
        loopinfo.power      = unique(datainfo.power(keep),'stable');
        if ~settings.debug
            clear i keep
        end
        % setup for different sweepcases
        loopinfo.sweepername = {'temp', 'field', 'power'}'; 
        loopinfo.sweepertest = [length(loopinfo.temp) length(loopinfo.field) length(loopinfo.power)];
        loopinfo.sweeper = {};
        for i = 1:3 % generate possible sweeper list
            if loopinfo.sweepertest(i) > 1
                loopinfo.sweeper(length(loopinfo.sweeper)+1) = loopinfo.sweepername(i);
            end
        end
        if ~settings.debug
            clear i 
            loopinfo = rmfield(loopinfo,{'sweepername','sweepertest'});
        end
        % one sweep parameter
        if length(loopinfo.sweeper) == 1
            loopinfo.usesweeper = loopinfo.sweeper;
            loopinfo.usetemp    = loopinfo.temp;
            loopinfo.usefield   = loopinfo.field;
            loopinfo.usepower   = loopinfo.power;
        % two or tree sweep parameters
        else
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
        if length(loopinfo.power) == 1
            loopinfo.usepower = loopinfo.power;    
        end
        if length(loopinfo.temp) == 1
            loopinfo.usetemp = loopinfo.temp;
        end
        if length(loopinfo.field) == 1
            loopinfo.usefield = loopinfo.field;
        end
    end % end ~settings.oldmode
    for usesweeperind = 1:length(loopinfo.usesweeper) % loop over selected sweeps
        usesweeper = loopinfo.usesweeper(usesweeperind);
        for k = 1:length(loopinfo.sweeper) % loop for fixed parameters
            if ~settings.oldmode
                if ~isequal(loopinfo.sweeper(k),usesweeper) % if there is something to choose, determine the case
                    if isequal(loopinfo.sweeper(k),{'temp'})
                        parameterarray  = loopinfo.temp;
                        paramcase       = 1;
                    elseif isequal(loopinfo.sweeper(k),{'field'})
                        parameterarray  = loopinfo.field;
                        paramcase       = 2;
                    elseif isequal(loopinfo.sweeper(k),{'power'})
                        parameterarray  = loopinfo.power;
                        paramcase       = 3;
                    end
                    % choose fixed parameters
                    if settings.multimodes && ismember(usename,loopinfo.multimodes)
                        [fixparameterindex,tf] = listdlg('PromptString',{['At which ' char(loopinfo.sweeper(k)) ...
                        ' do you want your ' char(usesweeper) 'sweep for ' char(usename) ','...
                        'multimode (' num2str(currentmultimodestep) '/' num2str(currentmultimodemax) ')?'],'',''},'ListString',parameterarray);
                    else
                        [fixparameterindex,tf] = listdlg('PromptString',{['At which ' char(loopinfo.sweeper(k)) ...
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
                    if isequal(loopinfo.sweeper(k),{'temp'})
                        loopinfo.usetemp  = loopinfo.temp;
                    elseif isequal(loopinfo.sweeper(k),{'field'})
                        loopinfo.usefield = loopinfo.field;
                    elseif isequal(loopinfo.sweeper(k),{'power'})
                        loopinfo.usepower = loopinfo.power;
                    end
                end
            end % end ~settings.oldmode
            if k == length(loopinfo.sweeper) % when previous loop went trough
                if ~settings.oldmode
                    if isequal(usesweeper,{'temp'}) % again casedetection with assigning
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
                counter = 0; % detailedprocessinfo
                % loop trough selected cases and get corresponding files
                for usefixedtemp = 1:fixedtemp
                    for usefixedfield = 1:fixedfield
                        for usefixedpower = 1:fixedpower
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
                            if settings.multifit >= length(currentfiles) % maxmultifit
                                settings.multifit = length(currentfiles);
                            end
                            file.end = importdata([datainfo.datapath,'\',char(currentfiles(end))]);
                            % init for the final loop trough files & fitmenu
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
                                    f1 = figure(1);
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
if ~settings.debug % clear variables - clearing directly after the last use could improve performance slightly
    clearvars -except datainfo loopinfo settings
end
disp('job done :-)')


% functions
function [lb,ub,x] = getbounds(data) % select range
    [x,~] = ginput(2);
    if x(1) < x(2)
        lb=find(data(:,1)>x(1),1,'first');
        ub=find(data(:,1)<x(2),1,'last');
    else
        lb=find(data(:,1)>x(2),1,'first');
        ub=find(data(:,1)<x(1),1,'last');
    end
end

function fitparameter = fit_complex(data,initial_guess,lb,ub,fastmode,initial_run,tau) % fitfunction with fixed tau after first fit
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