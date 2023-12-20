%% Preprocessing Pipeline based on PREP Pipeline
% Date: 06/11/2023
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

    %% Show the frequency spectrum of current EEG version

    [frefs, pows, BadChans] = showSpectrum(EEG, string({EEG.chanlocs(1:scalpNum).labels}), 1:scalpNum, 1:scalpNum, 'Spectra All Channels (raw)', scalpNum/2);

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

    %% Remove line noise using the Cleanline algorithm.

    fprintf('Line noise removal\n');

    [EEGClean, lineNoise] = removeLineNoise(EEGNew, UserParam);
    EEG.etc.noiseDetection.lineNoise = lineNoise;

    % Subtracting HP filtrage effect from line-noise corrected EEG data.
    % This avoids having to commit to a filtering strategy.
    EEG.data(lineNoise.lineNoiseChannels, :) = EEG.data(lineNoise.lineNoiseChannels, :) - EEGNew.data(lineNoise.lineNoiseChannels, :)...
        + EEGClean.data(lineNoise.lineNoiseChannels, :);

    %% save the clean lined data here.

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

    %% Compute the average reference before bad electrode detection.

    norefindx = find(~ismember({EEG.chanlocs.type}, 'EEG'));   % Only include the scalp elements in the

    if strcmp(UserParam.reference.type, 'average')
        EEG = pop_reref(EEG, [], 'exclude', norefindx);
    else
        refchans = cell2mat(eval(UserParam.reference.refchannels));
        EEG = pop_reref(EEG, refchans, 'keepref', 'on');
    end

    saveFname_ref = strcat(EEG.setname,'-ref1');
    EEG.setname = saveFname_ref;
    savepath      = fullfile(main_dir, 'Processed_data');
    EEG = pop_saveset( EEG, 'filename',saveFname_ref, 'filepath',savepath);

    %% Apply the clean_rawdata plugin to detect bad electrodes. This requires computing the reference before-hand.
    OldEEG = EEG;

    EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',5,'ChannelCriterion',0.87, 'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
        'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian', 'WindowCriterionTolerances',[-Inf 7] ,'fusechanrej',1);

    % EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',str2double(UserParam.cleancriteria.FlatlineCriterion),'ChannelCriterion',str2double(UserParam.cleancriteria.ChannelCriterion), ...
    %     'LineNoiseCriterion',str2double(UserParam.cleancriteria.LineNoiseCriterion),'Highpass',[0.25 0.75],'BurstCriterion',str2double(UserParam.cleancriteria.BurstCriterion), ...
    %     'WindowCriterion',str2double(UserParam.cleancriteria.WindowCriterion),'BurstRejection',UserParam.cleancriteria.BurstRejection,'Distance',UserParam.cleancriteria.Distance, ...
    %     'WindowCriterionTolerances', [-Inf 7] ,'fusechanrej',str2double(UserParam.cleancriteria.fusechanrej));

    noisyChannelsIndx = find([EEG.etc.clean_channel_mask == 0]);
    noisySamplesIndx = find([EEG.etc.clean_sample_mask == 1]);
    timeOld = OldEEG.times;
    EEG.etc.rejected_samples_time = timeOld(noisySamplesIndx);

    saveFname_clean = strcat(EEG.setname,'-cleanraw');
    EEG.setname = saveFname_clean;
    savepath      = fullfile(main_dir, 'Processed_data');
    EEG = pop_saveset( EEG, 'filename',saveFname_clean, 'filepath',savepath);

    %% Recompute the average reference, interpolaing the bad electrode and removing them.

    norefindx1 = find(~ismember({EEG.chanlocs.type}, 'EEG'));   % Only include the scalp elements in the
    EEG = pop_reref(EEG, [], 'exclude', norefindx1);                    % do not interpolate here. Interpolate after ICA rejection.

    saveFname_ref2 = strcat(saveFname_clean,'-reref');
    EEG                      = pop_saveset( EEG, 'filename',saveFname_ref2, 'filepath',savepath);

    %% Carry out filtering: detrend and low-pass.
    %    Need to create two versions of the data:
    %    - One version that is highpass filtered at 0.1Hz; the version that
    %    will be used for ERP analysis.
    %   - A second version that is highpass filtered at 1Hz, the version that
    %   will be used for ICA.
    global EEG_filt_lp;

    if str2double(postprocessParams.filter_highpass) == 1
        EEG = pop_eegfiltnew(EEG, str2double(postprocessParams.filter_highpass_cutoff), []);   % This applies the default eeglab filter based on windowed-sinc filter with a hamming window applied. Version for ERP analysis.
        EEG_filt_hp  = pop_eegfiltnew(EEG, 1, []);                % Version that will be used for ICA.
    else
        fprintf('No highpass filtering \n')
    end

    EEG.etc.postprocess.detrend.method  = 'high-pass';
    EEG.etc.postprocess.detrend_type      = 'Hamming windowed sinc FIR filter';
    EEG.etc.postprocess.detrend_cutoff    = 1;

    EEG_filt_hp.etc.postprocess.detrend.method = 'high-pass';
    EEG_filt_hp.etc.postprocess.detrend_type     = 'Hamming windowed sinc FIR filter';
    EEG_filt_hp.etc.postprocess.detrend_cutoff   = 1;

    if str2double(postprocessParams.filter_lowpass) == 1
        EEG           =  pop_eegfiltnew(EEG, [], str2double(postprocessParams.filter_lowpass_cutoff));
        EEG_filt_lp =  pop_eegfiltnew(EEG_filt_hp, [], str2double(postprocessParams.filter_lowpass_cutoff));
    else
        fprintf('No low pass filter applied. \n');
        EEG_filt_lp = EEG_filt_hp;
    end

    EEG.etc.postprocess.LPfilter                = 'Low pass';
    EEG.etc.postprocess.LPfilter_type       = 'Hamming windowed sinc FIR filter';
    EEG.etc.postprocess.LPfilter_cutoff     = str2double(postprocessParams.filter_lowpass_cutoff);
    EEG_filt_lp.etc.postprocess.LPfilter                = 'Low pass';
    EEG_filt_lp.etc.postprocess.LPfilter_type       = 'Hamming windowed sinc FIR filter';
    EEG_filt_lp.etc.postprocess.LPfilter_cutoff     = str2double(postprocessParams.filter_lowpass_cutoff);

    %%  Add code to save filtered data here.

    fname = EEG.filename;
    fname_split = split(fname, '.');
    saveFname_filt = strcat(fname_split{1,1},'-filtered');
    EEG.setname = saveFname_filt;
    EEG = pop_saveset( EEG, 'filename',saveFname_filt, 'filepath',Dirsave);

    %% Call of function to carry out postprocess steps.

    PREPpipeline_crex_prepost_process(EEG, UserParam.postprocess, savepath);
end









