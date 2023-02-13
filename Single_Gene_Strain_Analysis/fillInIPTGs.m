function allExptData = fillInIPTGs(allExptData);

    inds = find(~cellfun(@any, {allExptData.iptgs}));

    for i = 1:length(inds)   
        ind = inds(i);
        disp(allExptData(ind).exptfile);
        % Get the image files
        imfiles = lower(allExptData(ind).imfiles);
        tempiptgs = NaN(1, length(imfiles));
        % Erase .tif and other parts that may be present
    %     imfiles = erase(imfiles, '.tif');
        imfiles = erase(imfiles, 'plac');
        % iterate over the images
        for j = 1:length(imfiles)
            imname = imfiles{j};
            iptg_str_pat = '_\d*i_';
            repstr = regexpi(imname, iptg_str_pat, 'match');
            if ~isempty(repstr) & contains(imname, '-21_')
                % it's the new format, use this to get the iptg
                conc = repstr;
                conc = erase(conc, '_');
                conc = erase(conc, 'i');
                tempiptgs(j) = str2double(conc);
            else
                % Try other approaches
                
                % split on delimiter
                if contains(imname, '_')
                    % erase numbers representing the replicate
                    repstr = regexpi(imname, '_?[a-z]*\d*.tif', 'match');
                    if ~isempty(repstr)
        %                 if contains(repstr, 'iptg')
        %                     repstr = erase(repstr, 'iptg');
        %                 end
                        imname = erase(imname, repstr);
                    else
                        erase(imname, '.tif');
                    end
                    imname_parts = split(imname, '_');
                else
                    % Maybe has spaces instead?
                    repstr = regexpi(tempstr, '\s*\w*\d.tif', 'match');
                    if ~isempty(repstr)
                        imname = erase(imname, repstr);
                    else
                        erase(imname, '.tif');
                    end
                    imname_parts = split(imname);
                end
            % Find 'mm'
                if any(contains(imname_parts, 'mm'))
                    conc = imname_parts{contains(imname_parts, 'mm')};
                    if contains(conc, 'pt')
                        % need to replace with decimal point
                        conc = strrep(conc, 'pt', '.');
                    end
                    conc = erase(conc, 'mm');            
                    tempiptgs(j) = str2double(conc);
                else
                    % contains something else or no number; look to see if we have
                    % any numbers that are just on their own now?
                    conc = cellfun(@str2num, ...
                            imname_parts, 'UniformOutput', false);
                    if length(find(~cellfun(@isempty, conc)))==1
                        tempiptgs(j) = conc{~cellfun(@isempty, conc)};
                    end

                end
            end
            if isnan(tempiptgs(j))
                disp(imfiles{j});
                tempiptgs(j) = input('Type in iptg number: ');
            end
        % store iptg
        end
        allExptData(ind).iptgs = tempiptgs;
    end

end
