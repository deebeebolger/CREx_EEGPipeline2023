function EEG_filt_lp = PREPpipeline_crex_rejcontinuous(EEG_filt_lp,notrig,triggers)
% Function that automatically detects time intervals of noisy data and
% marks them for rejection.
% Input variables
% 1. EEG structure
% 2. notrig: 1- if you wish to mark regions that do not include certain
% triggers, 0 - if you don't mind marking regions including triggers.
% 3. trig: matrix of trigger values to avoid when marking regions to reject.
% Output variables:
% 1. matrixrej_new : matrix of time intervals marked for rejection after
% visual verification.
%******************************************************************************

%% Need to check that triggers = [] if notrig == 0, if not correct.
%declare as global
shortisi= 50;
shortisisam  = floor(shortisi*EEG_filt_lp.srate/1000);  % to samples

[winrej,chanrej]      = basicrap_prepcrex(EEG_filt_lp,1:EEG_filt_lp.nbchan,50,2000,100,0,[],[]);
assignin('base', 'winrej', winrej)
[winrej2, chanrej2] = joinclosesegments_prepcrex(winrej, chanrej, shortisisam);

if size(chanrej2, 1)~=size(winrej2, 1)
    winrej2 = winrej2(1:length(chanrej2),:);
end

%% CHECK THAT MARKED REGIONS DO NOT COINCIDE WITH DEFINED TRIGGERS

if notrig == 1
    y=zeros(length(winrej2),length(triggers));
    for counter=1:length(triggers)

        x=([EEG_filt_lp.event.type]==triggers(counter));
        lats=[EEG_filt_lp.event(x).latency];
        x2=zeros(length(lats),length(winrej2));

        for counter1=1:length(winrej2)                             %for each rejection interval, determine if it coincides with trigger code.
            for counter2=1:length(lats)
                x2(counter2,counter1)=(lats(counter2)>=winrej2(counter1,1) & lats(counter2)<=winrej2(counter1,2));
            end
            y(counter1,counter)=isempty(find(x2(:,counter1)));

        end
    end

    %y tells us if a trigger has fallen within a certain window and which
    %trigger and the index of which window
    Y=(sum(y,2)==length(triggers));
    winrej2=winrej2(Y,:);
    chanrej2=chanrej(Y,:);

elseif notrig ==0
    fprintf('Do not take into account event triggers when marking noisy time intervals.\n');
end

%% Plot the EEG highlighting those time-intervals marked as Bad.

colorseg     = [1.0000    0.9765    0.5294];
colormatrej = repmat(colorseg, size(winrej2,1),1);
commrej     = sprintf('EEG_filt_lp = eeg_eegrej(EEG_filt_lp, TMPREJ(:,1:2)); eeglab redraw;');
matrixrej     = [winrej2 colormatrej chanrej2];
try
    waitforbuttonpress
    eegplot(EEG_filt_lp.data, 'winrej', matrixrej, 'srate', EEG_filt_lp.srate,'butlabel','REJECT','command', commrej,'events', EEG_filt_lp.event,'winlength', 50);
catch
    disp(TMPREJ)
end
end

