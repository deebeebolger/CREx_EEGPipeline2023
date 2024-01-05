function CREx_EEGPipeline()

%% Preprocessing Pipeline for EEG.
% Date: 19/12/2023
% Programmed by: D. Bolger
%
% PREP pipeline reference: Bigdely-Shamlo N, Mullen T, Kothe C, Su K-M and Robbins KA (2015)
% The PREP pipeline: standardized preprocessing for large-scale EEG analysis
% Front. Neuroinform. 9:16. doi: 10.3389/fninf.2015.0001
%****************************************************************************************************************

%% Select all the data files that you wish to process.
[filename, filepath] = uigetfile('*.bdf', 'Choose a BDF data file ', 'multiselect', 'on');
[parentdir, childdir, ext] = fileparts(filepath);
main_dir = fileparts(parentdir);

%% Load in the main parameters file.

param_path = fullfile(main_dir, "PipelineParameters.txt");
opts = detectImportOptions(param_path);
opts = setvartype(opts, 'Value', 'string');
PIn = readtable(param_path, opts);          % Load in the parameters in table form.

%% Select the Bin Descriptor File to use when creating the EventLists.

[fileBDF, pathBDF] = uigetfile('*.txt','Select the Bin Descriptor File (.txt) to use when adding EventList.');
BinDFile_fullpath   = fullpath(pathBDF, fileBDF);

%% Loop through the participant data files.

for fcount = 1:length(string(filename))

    if ischar(filename)
        [EEG, ~, dat] = pop_biosig(fullfile(filepath,filename));
        EEG.filename = filename;
        currfname       = filename;
    else

        [EEG, ~, dat] = pop_biosig(fullfile(filepath,filename{1,fcount}));
        EEG.filename = filename{1,fcount};
        currfname       = filename{1, fcount};
    end

    EEG.filepath   = filepath;
    fname_split = split(currfname, '.');
    EEG.setname = fname_split{1};

    %% Add the title of the current dataset to the parameters table.

    addSujInfo = {{'Participant'}, {'Title'}, {EEG.setname}; {'Participant'}, {'File'}, {EEG.filename}};
    PIn = [addSujInfo; PIn];

    %% Use the EEG structure corresponding to the first raw file to create json file and structure with parameters

    params = [];                                                                                                  %Initialize the structure.
    EEG.data = double(EEG.data);                                                                    % Set eeg data to double precision.
    UserParam = make_parameters_json(PIn, params, main_dir, EEG);         % Call of function to establish the parameters and create a json file with the parameters. It saves a json file of the user parameters.

    %% Define scalp, external and auxiliary channels. Define channels to remove.

    scalpNum       =  str2double(UserParam.channels.scalpChannelNumber);
    channelScalp = 1: scalpNum;

    exgNum = str2double(UserParam.channels.exgChannelNumber);
    if ~isempty(exgNum) | exgNum > 0
        channelEXG = channelScalp(end)+(1 : exgNum); end
    if ~isempty(UserParam.channels.auxChannels)
        channelAux    = UserParam.channels.auxChannels; end

    ChanRej = UserParam.channels.noChannel; % Get channels to reject.
    fprintf('Removing the following channels from the dataset: %s', string(ChanRej));
    EEG = pop_select(EEG, 'nochannel', cellstr(ChanRej));

    %%  Add channel info to the EEG structure.

    if isempty(cell2mat({EEG.chanlocs.theta})) == 1
        chlocpath = fullfile(main_dir,'Chanlocs_files/standard-10-5-cap385.elp');
        EEG         = pop_chanedit(EEG, 'lookup',chlocpath);            % Load channel path information
        strparts    = split(chlocpath, '/');

        fprintf('\n--------------------------\n')
        fprintf("Adding channel information based on file %s.\n", strparts{length(strparts),1})
    else
        fprintf('\n--------------------------\n')
        fprintf("Channel locations appear to have already been added.\n");
    end

    

    %% Call of function Add_EventList() to add EventList to current EEG structure.
    %    This function uses the binlister() from ERPLAB which takes in a
    %    Bin Descriptor File (BDF)...not to be confused with Biosemi Data
    %    Format. Event codes are added at this point before potential data
    %    interval rejection, which could upset the structure of the trigger
    %    coding.

    trigs_forbidden = [];
    trigs_ignore      = [];

    EEG = Add_EventList(EEG, BinDFile_fullpath, savepath, trigs_forbidden, trigs_ignore);

    %% Carry out downsampling if defined.

    if str2double(UserParam.resample.resampleOff) == 1
        fprintf('Do not carry out downsampling.\n')
    elseif str2double(UserParam.resample.resampleOff) == 0
        newFS = str2double(UserParam.resample.resamplingFrequency);
        fprintf('Downsampling from %d to %d Hz', EEG.srate, newFS);
        EEG = pop_resample(EEG, newFS);
    end

    %% Detrending:  High-pass filter the signals to facilitate bad electrode detection.
    %    This is only necessary if we are applying the CleanLine algorithm
    %    to correct line noise artefact.

    fprintf('Preliminary detrend before applying the line-noise correction.\n');

    [EEGNew, detrend] = removeTrend(EEG, UserParam);
    EEG.etc.noiseDetection.detrend = detrend;

    %% Show the frequency spectrum of current EEG version

    [frefs, pows, BadChans] = showSpectrum(EEG, string({EEG.chanlocs(1:scalpNum).labels}), 1:scalpNum, 1:scalpNum,...
        'Channel Spectra before linenoise correction', scalpNum/2);

    %% Remove line noise using the Cleanline algorithm.

    fprintf('Line noise removal\n');

    [EEGClean, lineNoise] = removeLineNoise(EEGNew, UserParam);
    EEG.etc.noiseDetection.lineNoise = lineNoise;

    % Subtracting HP filtrage effect from line-noise corrected EEG data.
    % This avoids having to commit to a filtering strategy.
    EEG.data(lineNoise.lineNoiseChannels, :) = EEG.data(lineNoise.lineNoiseChannels, :) - EEGNew.data(lineNoise.lineNoiseChannels, :)...
        + EEGClean.data(lineNoise.lineNoiseChannels, :);

    %% Save the clean lined data here.

    if ~exist(fullfile(main_dir, 'Processed_data'), 'dir')
        mkdir(fullfile(main_dir, 'Processed_data'))
    else
        fprintf('Processed_data folder already exists. \n')
    end

    currfname_split = split(currfname, '.');
    saveFname_lnoise = strcat(currfname_split{1,1},'-linenoise');
    EEG.setname = saveFname_lnoise;
    savepath      = fullfile(main_dir, 'Processed_data');
    EEG = pop_saveset( EEG, 'filename',saveFname_lnoise, 'filepath',savepath);

     %% Show the frequency spectrum of current EEG version

    [frefs, pows, BadChans] = showSpectrum(EEG, string({EEG.chanlocs(1:scalpNum).labels}), 1:scalpNum, 1:scalpNum,...
        'Channel Spectra after linenoise correction', scalpNum/2);

    %% Compute the average reference before bad electrode detection.

    norefindx = find(~ismember({EEG.chanlocs.type}, 'EEG'));   % Only include the scalp elements in the

    if strcmp(UserParam.reference.type, 'average')
        EEG = pop_reref(EEG, [], 'exclude', norefindx);
    else
        refchans = cell2mat(eval(UserParam.reference.refchannels));
        EEG = pop_reref(EEG, refchans, 'keepref', 'on');                          % Keeping the original reference channels. 
    end

    saveFname_ref = strcat(EEG.setname,'-ref1');
    EEG.setname = saveFname_ref;
    savepath      = fullfile(main_dir, 'Processed_data');
    EEG = pop_saveset( EEG, 'filename',saveFname_ref, 'filepath',savepath);

    %% Apply the clean_rawdata plugin to detect bad electrodes. This requires computing the reference before-hand.
    OldEEG = EEG;
    chansPreclean = OldEEG.chanlocs;

    EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',5,'ChannelCriterion',0.87, 'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
        'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian', 'WindowCriterionTolerances',[-Inf 7] ,'fusechanrej',1);

    % EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',str2double(UserParam.cleancriteria.FlatlineCriterion),'ChannelCriterion',str2double(UserParam.cleancriteria.ChannelCriterion), ...
    %     'LineNoiseCriterion',str2double(UserParam.cleancriteria.LineNoiseCriterion),'Highpass',[0.25 0.75],'BurstCriterion',str2double(UserParam.cleancriteria.BurstCriterion), ...
    %     'WindowCriterion',str2double(UserParam.cleancriteria.WindowCriterion),'BurstRejection',UserParam.cleancriteria.BurstRejection,'Distance',UserParam.cleancriteria.Distance, ...
    %     'WindowCriterionTolerances', [-Inf 7] ,'fusechanrej',str2double(UserParam.cleancriteria.fusechanrej));

    noisyChannelsIndx = find([EEG.etc.clean_channel_mask == 0]);
    noisySamplesIndx  = find([EEG.etc.clean_sample_mask == 1]);
    timeOld = OldEEG.times;
    EEG.etc.rejected_samples_time = timeOld(noisySamplesIndx);

    saveFname_clean = strcat(EEG.setname,'-cleanraw');
    EEG.setname = saveFname_clean;
    savepath      = fullfile(main_dir, 'Processed_data');
    EEG = pop_saveset( EEG, 'filename',saveFname_clean, 'filepath',savepath);

    %% Recompute the average reference, interpolaing the bad electrode and removing them.
    
    if strcmp(UserParam.reference.type, 'average')
        norefindx1 = find(~ismember({EEG.chanlocs.type}, 'EEG'));   % Only include the scalp elements in the
        EEG = pop_reref(EEG, [], 'exclude', norefindx1);                     % do not interpolate here. Interpolate after ICA rejection.
    
        saveFname_ref2 = strcat(saveFname_clean,'-reref');
        EEG                      = pop_saveset( EEG, 'filename',saveFname_ref2, 'filepath',savepath);
    else
        fprintf('No need to recalculate the reference as the %s has already been applied. \n', UserParam.reference.type);
    end 

    %% Carry out filtering: high-pass (if defined).
    %    The data set will be highpass filtered at 1Hz to facilitate ICA.
    %    This yields 2 versions of the dataset: one highpass filtered at
    %    1Hz (for ICA) and one highpass filtered at 0.25Hz (for ERP analysis).

    if isempty(UserParam.cleancriteria.Highpass)
        fprtintf('Data was not high-pass filtered during data cleaning.\n Carry out high-pass filtering.\n');
        EEG = pop_eegfiltnew(EEG, str2double(postprocessParams.filter_highpass_cutoff), []);

        EEG.etc.postprocess.detrend.method  = 'high-pass';
        EEG.etc.postprocess.detrend_type      = 'Hamming windowed sinc FIR filter';
        EEG.etc.postprocess.detrend_cutoff    = str2double(UserParam.postprocess.filter_highpass_cutoff);

        saveFname_filterHP = strcat(saveFname_ref2,'-filterHP');
        EEG = pop_saveset( EEG, 'filename',saveFname_filterHP, 'filepath',savepath);

    elseif ~isempty(UserParam.cleancriteria.Highpass)
        fprintf('Data already highpass filtered during data cleaning.\n')
        hpcut = str2double(UserParam.cleancriteria.Highpass);              % Fetch highpass cutoff applied during data cleaning.

        EEG.etc.postprocess.detrend.method  = 'high-pass applied during data cleaning';
        EEG.etc.postprocess.detrend_type      = 'Hamming windowed sinc FIR filter';
        EEG.etc.postprocess.detrend_cutoff    = hpcut(1);
    end

    EEG_filterHP  = pop_eegfiltnew(EEG, 1, []);   % Apply agressive highpass filtering at 1Hz to facilitate ICA.

    %% Carry out filtering: lowpass (if defined)

    if str2double(UserParam.postprocess.filter_lowpass) == 1
        EEG           =  pop_eegfiltnew(EEG, [], str2double(postprocessParams.filter_lowpass_cutoff));
        EEG_filterLP =  pop_eegfiltnew(EEG_filterHO, [], str2double(UserParam.postprocess.filter_lowpass_cutoff));

        EEG.etc.postprocess.LPfilter                = 'Low pass';
        EEG.etc.postprocess.LPfilter_type       = 'Hamming windowed sinc FIR filter';
        EEG.etc.postprocess.LPfilter_cutoff     = str2double(UserParam.postprocess.filter_lowpass_cutoff);

        saveFname_filterLP = strcat(saveFname_filterHP,'-filterLP');
        EEG = pop_saveset( EEG, 'filename',saveFname_filterLP, 'filepath',savepath);
    else
        fprintf('No low pass filter applied. \n');
        EEG_filterLP = EEG_filterHP;
    end

    %% Run Independent Components Analysis (ICA)

    if str2double(UserParam.postprocess.doICA) == 1

        if strcmp(UserParam.reference.type, 'average')
            fprintf('The average reference reduces the data rank by one leading to rank deficiency. Need to remove one channel before ICA.\n')
            chans2remove = find(contains({EEG.chanlocs.labels}, 'T'));
            removeChan = randsample(chans2remove, 1);
            fprintf('Removing channel %s before ICA to render data full ranked.\n', EEG.chanlocs(removeChan).labels);

            EEG_filterLP = pop_select(EEG_filterLP, 'rmchannel',EEG.chanlocs(removeChan).labels);
            iseeg = find(ismember({EEG_filterLP.chanlocs.type}, 'EEG' ));   % Find the number of scalp channels.

        else
            fprintf('Reference type applied is %s.\nNo need to remove the channels. ', UserParam.reference.type)  % Question here of whether to include the mastoid reference channels in ICA.
            find(ismember({EEG_filterLP.chanlocs.type}, 'EEG' ));   % Find the number of scalp channels.
        end

        icaType = UserParam.postprocess.icaType;

        switch icaType
            case 'infomax'
                fprintf('Carry out Independent Components Analysis (ICA) by applying the %s algorithm', UserParam.postprocess.icaType)
                ncomps = iseeg;   % Define the number ICs
                [weights, sphere] = runica(EEG_filterLP.data(iseeg,:), 'ncomps', ncomps);

                % Copy the IC weigths and sphere information to EEG
                % dataset.
                EEG.icaweights = weights;
                EEG.icasphere  = sphere;
                EEG.icachansind = iseeg;
                EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(iseeg,:); % Calculate the ICA activations.
                EEG.icawinv = inv(weights*sphere);

                saveFname_ICA = strcat(saveFname_filterLP,'-ica');
                EEG = pop_saveset( EEG, 'filename',saveFname_ICA, 'filepath',savepath);



                % Use the ICLabel plugin to automatically identify component corresponding to artifacts.
                [EEG, icOut]   = pop_iclabel(EEG, 'default');
                icThreshold    = [0 0;0 0; 0.8 1; 0 0; 0 0; 0 0; 0 0];

                EEG = pop_icflag(EEG, icThreshold);
                ic2Rej = find(EEG.reject.gcompreject);        % Reject component/s

                saveFname_icaRej = strcat(saveFname_ICA,'-icarej');
                EEG = pop_saveset( EEG, 'filename',saveFname_icaRej, 'filepath',savepath);

            case 'amica'
                fprintf('Carry out Independent Components Analysis (ICA) by applying the %s algorithm', UserParam.postprocess.icaType)
                % to be continued...

        end % end of switch.

    end % if doICA loop.

    %% Carry out interpolation: interpolate the bad channels rejected during data cleaning.
    %    Spherical spline interpolation is applied.

    EEG = pop_interp(EEG, chansPreclean, UserParam.interpolation.type);

    saveFname_ssInterp = strcat(saveFname_icaRej,'-ssinterp');
    EEG = pop_saveset( EEG, 'filename',saveFname_ssInterp, 'filepath',savepath);

    %% Segment the continuous data.

    poststim = [-1000 800];
    prestim  = [-200 0];
    EEG_epoched = pop_epochbin( EEG , poststim,  prestim);
    namefig = EEG_epoched.setname;
    pop_plotepoch4erp(EEG_epoched, namefig)
    
    save_fname = strcat(EEG_epoched.setname,'_ToneLDT_1st_epoched');
    EEG = pop_saveset( EEG_epoched, 'filename',save_fname, 'filepath',Dirsave);

end









