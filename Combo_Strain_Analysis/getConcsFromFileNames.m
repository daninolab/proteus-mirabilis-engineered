function [IPTGs, aras] = getConcsFromFileNames(filenames)
    % Filenames should be cell array of strings representing image file
    % names; get rid of capitals
    filenames = lower(filenames);
   
    
    disp(strcat("Example file name: ", filenames{1}));
    
    iptg_conc_str = input("Type string for iptg conc, eg '0pt1_iptg': ", 's');
    ara_conc_str = input("Type string for ara conc, eg '0pt1_ara': ", 's');

    iptg_splitchar = getSplitChar(iptg_conc_str);
    ara_splitchar = getSplitChar(ara_conc_str);
    
    % Get indices of concentration, iptg string
    iptg_concind = getConcInd(iptg_conc_str, iptg_splitchar);
    ara_concind = getConcInd(ara_conc_str, ara_splitchar);
    
    % Get the string matching iptg, ara
    iptg_str = split(iptg_conc_str, iptg_splitchar);
    iptg_str = iptg_str{contains(iptg_str, 'i', 'IgnoreCase', 1)};
    ara_str = split(ara_conc_str, ara_splitchar);
    ara_str = ara_str{contains(ara_str, 'a', 'IgnoreCase', 1)};
    
    % Using the string matching each molecule + the delimiter, construct a
    % regexp
    if iptg_concind==1
        iptg_matchexp = strjoin({iptg_splitchar, '\d*\.*\d*', iptg_splitchar, iptg_str}, '');
        ara_matchexp = strjoin({ara_splitchar, '\d*\.*\d*', ara_splitchar, ara_str}, '');
    else
        % the iptg/ara come before the concentration
        iptg_matchexp = strjoin({iptg_splitchar, iptg_str, iptg_splitchar, '\d*\.*\d*'}, '');
        ara_matchexp = strjoin({ara_splitchar, ara_str, ara_splitchar, '\d*\.*\d*'}, '');
    end

    % Having obtained these indices & split characters, use the getConcNum
    % function on each image name in the list to get the IPTG & ara
    IPTGs = zeros(1, length(filenames));
    aras = zeros(1, length(filenames));
    
    for i = 1:length(filenames)
        % Check if filename has 'pt' instead of '.' and replace if so
        tempfile = filenames{i};
       if(contains(erase(tempfile, iptg_str), 'pt'))
           tempexp = '\d*pt\d*';
           tempfile = erase(filenames{i}, lower(iptg_str));
           tempstr = regexpi(tempfile, tempexp, 'match'); % eg '0pt2'
           tempfile = replace(filenames{i}, tempstr, replace(tempstr, 'pt', '.'));
       end
        
        % Extract representative strings
        tempiptg_str = regexpi(tempfile, iptg_matchexp, 'match');
        tempara_str = regexpi(tempfile, ara_matchexp, 'match');
        
        % If we get the delimiter at the beginning of this expression,
        % remove it
        if strcmpi(tempiptg_str{1}(1), iptg_splitchar)
            tempiptg_str = tempiptg_str{1}(2:end);
        end
        if strcmpi(tempara_str{1}(1), ara_splitchar)
            tempara_str = tempara_str{1}(2:end);
        end
        
        % Extract & store the concentrations
        IPTGs(i) = getConcNum(tempiptg_str, iptg_splitchar, iptg_concind);
        aras(i) = getConcNum(tempara_str, ara_splitchar, ara_concind);
    end
    


    

end

%%


function splitchar = getSplitChar(test_str)
    % Find the delimiter used in filenames
    if contains(test_str, '_')
        splitchar = '_';
    else
        splitchar = input("Type delimiter character (eg . or space)", 's');
        if isempty(splitchar)
            splitchar = " ";
        end
        
    end
end

function concind = getConcInd(test_str, splitchar)
    concind = [];
    % Separate based on delimiter
    temp_split = split(test_str, splitchar);
    
    % Using the split, check which is the number
    for i = 1:length(temp_split)
        if any(isstrprop(temp_split{i}, 'digit'))
            % This is the location of the concentration number
            concind = i;
        end
    end

end


function matchexp = makeRegExp(test_str, splitchar, concind)
    % Given the test string & knowing where to split it, generate a regular
    % expression to use for matching
    splitstrs = split(test_str, splitchar);
    % The concentration should be represented by a digit or set of digits
    splitstrs{concind} = regexptranslate('flexible', splitstrs{concind}, '\d*');
    % Remake the expression now that that one has been replaced
    matchexp = strjoin(splitstrs, splitchar);
    
end


function concnum = getConcNum(test_str, splitchar, concind)
    % Get the actual concentration out of the string
    temp_split = split(test_str, splitchar);
    concstr = temp_split{concind};
    concstr = strrep(concstr, 'pt', '.');
    concnum = str2double(concstr);

end