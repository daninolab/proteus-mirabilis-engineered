function allExptData = fillInIptgAra(allExptData)
    % For allExptData, the new combo gene flattening script does not fill
    % in iptgs & aras; fill this in for any expts in allExptData that need
    % it
    
    % Look for expts in allexptdata that do not have iptgs & aras filled
    % out
    if isfield(allExptData, 'iptgs')
        exptinds = find(cellfun('isempty', {allExptData.iptgs}));
    else
        exptinds = [];
    end
    
    if ~any(exptinds)
        % all filled in--ask user if they want to redo?
        chk = input('All iptgs & aras already calculated. Redo all? Type y/n: ', 's');
        if strcmpi(chk, 'y')|strcmpi(chk, 'yes')
           % If result was okay, do the rest automatically
           disp('Redoing all iptgs & aras in allExptData. ');
           exptinds = 1:length(allExptData);
        end
    end

    for i = 1:length(exptinds)
        % Load the img names
        disp(allExptData(exptinds(i)).date);
        imfiles = allExptData(exptinds(i)).imfiles;
        disp(imfiles');

        % Use function getConcsfromFileNames
        [iptgs, aras] = getConcsFromFileNames(imfiles);

        % Display results
        fprintf('Results obtained are: \n');
        disp([iptgs', aras']);

        % Display result; ask user to confirm or do it manually

        chk = input("Results look okay? y/n: ", 's');

        if strcmpi(chk, 'y')|strcmpi(chk, 'yes')
           % If result was okay, do the rest automatically
           disp('Using these iptgs & aras in allExptData. ');
        else
            % if not okay, iterate over images, ask user to type in iptg/ara for
            % each
            disp('Type iptgs & aras in manually');
            for k = 1:length(imfiles)
               tempimname = imfiles{k};
               disp(tempimname);
               iptgs{k} = input("Type iptg: ");
               aras{k} = input("Type ara: ");
            end
        end
        % Store results in allExptData
        allExptData(exptinds(i)).iptgs = iptgs;
        allExptData(exptinds(i)).aras = aras;
    end


end

function [iptg, ara] = getNums(tempimname)
    matchStr = regexpi(tempimname, '\d*\.*\d*_i', 'match');
    matchStr = matchStr{1};
    iptg = str2double(erase(matchStr, '_i'));
    matchStr = regexpi(tempimname, '\d*\.*\d*_a', 'match');
    matchStr = matchStr{1};
    ara = str2double(erase(matchStr, '_a'));

end