function PREPpipeline_crex_prepost_process(EEG, postprocessParams, Dirsave)
%% *********************************************************************************************************************
% Date: November 2023        Programmed by: D. Bolger
% This function carries two post-pre-processing steps:
% Filtering: HP or band-pass.
% If bandpass filtering is applied, high-pass and lowpass filtering is
% applied seperately.
% Removal of interpolated bad channels in preparation for ICA calculation.
% Input Data:
% Output Data:
%************************************************************************************************************************

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

%% Automatic identification of noisy time intervals on dataset 1 (conservative low-pass filtered data).

% notrig = 0; % 0 if you don't mind marking regions with triggers.
% trigs2avoid = [];
% 
% EEG_filt_lp = PREPpipeline_crex_rejcontinuous(EEG_filt_lp, notrig,trigs2avoid);  % Remove detected noisy data intervals from continuous temp EEG dataset.

% EEG = eeg_eegrej(EEG, TMPREJ(:,1:2));                                                               % Remove the same data intervals from the EEG dataset.
% 
% EEG_filt_lp.etc.postprocess.badtimewindows = TMPREJ;  % The first two columns detail the time windows marked for rejection (in samples).
% EEG_filt_lp.etc.postprocess.avoidtrigs = notrig;
% EEG_filt_lp.etc.postprocess.trigsavoid = trigs2avoid;
% 
% % Apply the noisy time intervals to the EEG data HP filtered for ERP
% % analysis.
% EEG.etc.postprocess.badtimewindows = TMPREJ;
% EEG.etc.postprocess.avoidtrigs = notrig;
% EEG.etc.postprocess.trigs2avoid = trigs2avoid;

%% Run ICA on data that has undergone the following processing:
%    - Detrended with high-pass cutoff of 1Hz.
%    - Interpolated electrodes removed.
%    - Noisy time intervals removed.

ica_type = 'infomax';   % Other options informax

iseeg = find(ismember({EEG_filt_lp.chanlocs.type}, 'EEG' ));

if strcmp(postprocessParams.icaType, 'amica')  % If the chosen ICA method is Adaptive Mixture ICA (AMICA)
    % Define parameters to run amica
    numprocs  = 1;            % Define the number of nodes
    threadmax = 1;            % Define the number of threads
    modelnum = 1;             % Define the number of ICA mixture models.
    itermax      = 1000;       % Define the number of learning steps.

    % Can calculate the rank with "pcakeep' option.

    [weights, sphere, mods] = runamica15(EEG_filt_lp.data, 'num_models', modelnum, 'outdir', savedir, 'numprocs', numprocs, 'max_threads', threadmax, 'max_iter', itermax);
elseif strcmp(postprocessParams.icaType, 'infomax') % If the chosen ICA method is Infomax ICA algorithm.
    % Define parameters to run infomax.
    ncomps = length(iseeg);

    % Can calculate the rank of the data using the pca option.
    [weights, sphere] = runica(EEG_filt_lp.data(iseeg,:), 'ncomps', ncomps);
end

% Copy the ICA weights and sphere information from EEG_filt_lp to the
% EEG highpass filtered for ERP analysis.

% chansnot_interp = find(~ismember(1:EEG.nbchan, interpolatedChannels));
% chansall = chansnot_interp(iseeg);
EEG.icaweights = weights;
EEG.icasphere  = sphere;
EEG.icachansind = iseeg;
EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(iseeg,:); % Calculate the ICA activations.
EEG.icawinv = inv(weights*sphere);

% Can use ICLabel plugin to automatically detect and l-mark the ICs to
% reject.
[EEG, icOut]   = pop_iclabel(EEG, 'default');
icThreshold    = [0 0;0 0; 0.8 1; 0 0; 0 0; 0 0; 0 0];

EEG               = pop_icflag(EEG, icThreshold);
ic2Rej            = find(EEG.reject.gcompreject);

%% Save dataset with ICA components.
saveFname_ica = strcat(saveFname_filt,'-ica');
EEG.setname = saveFname_ica;
EEG = pop_saveset( EEG, 'filename',saveFname_ica, 'filepath',Dirsave);

%%
% Visualise this component marked for rejection
% Put the code here...
% Reject the identified bad components. By defaul it identifies
% eye-movement related ICs.

EEG = pop_subcomp(EEG,ic2Rej,1);

%% Interpolate here.

%% Save ICA corrected data to file.

saveFname_icarej = strcat(saveFname_filt,'-icarej');
EEG.setname = saveFname_icarej;
EEG = pop_saveset( EEG, 'filename',saveFname_icarej, 'filepath', Dirsave);


%% Call of erplab function to add EVENTLIST field to the current EEG structure.

event_dir = fullfile(filesep, 'Users', 'bolger', 'matlab', 'Projects', 'PREPpipeline_crex', 'EventData');
bdffile_name = 'ToneLDT_CAT_v2_1st.txt';
bdffile_path = fullfile(event_dir, bdffile_name);

eventlist_path = fullfile(Dirsave, strcat(EEG.setname, '_eventslist.txt'));
EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Eventlist', eventlist_path );

forbiddenCodeArray = [];
ignoreCodeArray = [];
[EEG, EVENTLIST, ~, ~] = binlister(EEG, bdffile_path, 'none', 'none', forbiddenCodeArray, ignoreCodeArray, 1);
EEG. EVENTLIST = [];
EEG.EVENTLIST = EVENTLIST;

%% Segment the continuous data.

poststim = [-1000 800];
prestim  = [-200 0];
EEG_epoched = pop_epochbin( EEG , poststim,  prestim);
namefig = EEG_epoched.setname;
pop_plotepoch4erp(EEG_epoched, namefig)

save_fname = strcat(EEG_epoched.setname,'_ToneLDT_1st_epoched');
EEG = pop_saveset( EEG_epoched, 'filename',save_fname, 'filepath',Dirsave);














