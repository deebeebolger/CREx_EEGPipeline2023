function bigStruct = make_parameters_json(paramFile, structIn, savedir, signal)
% takes a cell array in which the starting parameters are defined by the
% user.
% Input: N X 2 cell array
%*************************************************************************************

params_titles = paramFile(:, 1);                                                                                      % Extract the parameter titles

%% Extract the actions, parameters and values from the parameter table.

actionCol = paramFile{:,1};
paramCol = paramFile{:,2};
valueCol = paramFile{:,3};

%%  Identify the actions possible and create a parameters structure based on table contents. 

actionsAll = unique(actionCol);
bigStruct = structIn; 

for counter = 1:length(actionsAll)

    I = cellfun(@(s) contains(s, actionsAll{counter, 1}), actionCol, 'UniformOutput',false);
    indxCurr = find(cell2mat(I));
    paramCurr = paramFile{indxCurr, 2};
    valueCurr  = paramFile{indxCurr, 3};
    
    % Ensure that the parameters and value columns are in string format
    % (and not cell).
    if ~strcmp(class(paramCurr), 'string')
        fprintf('The parameter column is of class %s \n', class(paramCurr));
        paramCurr = string(paramCurr);
    elseif strcmp(class(paramCurr), 'string')
        fprintf('The parameter column in an array of strings');
    end

    if  ~strcmp(class(valueCurr), 'string')
        fprintf('The value column is of class %s \n', class(valueCurr));
        valueCurr = string(valueCurr);
     elseif strcmp(class(valueCurr), 'string')
        fprintf('The value column in an array of strings');
    end

    % Create structure to house the parameters.
    for pvcount = 1:length(paramCurr)

        if contains(valueCurr{pvcount, 1}, 'signal')
            currVal = eval(valueCurr{pvcount, 1});
            bigStruct.(actionsAll{counter, 1}).(paramCurr{pvcount, 1}) = currVal;
        else
            bigStruct.(actionsAll{counter, 1}).(paramCurr{pvcount, 1}) = valueCurr{pvcount,1};
        end

    end 

end

%% Create a json file to store the parameters.

json_params = jsonencode(bigStruct, PrettyPrint=true);
json_title = 'userParams.json'
fid = fopen(fullfile(savedir, json_title), 'w');
fprintf(fid, '%s', json_params);
fclose(fid);
































