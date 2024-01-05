function EEG = Add_EventList(EEG, bdffile_path, dirsave, code_forbidden, code_ignore)

    %% Call of erplab function to add EVENTLIST field to the current EEG structure.

    
    eventlist_path = fullfile(dirsave, strcat(EEG.setname, '_eventslist.txt'));
    EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Eventlist', eventlist_path );
    
    forbiddenCodeArray = code_forbidden;
    ignoreCodeArray = code_ignore;
    [EEG, EVENTLIST, ~, ~] = binlister(EEG, bdffile_path, 'none', 'none', forbiddenCodeArray, ignoreCodeArray, 1);
    EEG. EVENTLIST = [];
    EEG.EVENTLIST = EVENTLIST;

end