function dateFormatted = formatDateStr(datestr)
    % take in a string representing a timelapse date in mm-dd-yy format and
    % then turn it into a datevec
    if ~isdatetime(datestr)
        datenums = str2double(split(datestr, ["_", "-"]));
        if length(datenums)==3
            datearray = [datenums(3) datenums(1) datenums(2)];
            if datearray(1)<2000
                datearray(1) = datearray(1)+2000;
            end
                
            dateFormatted = datetime(datearray);
        else
            disp(datestr);
            datestr_temp = input('Problem with date. Type it manually in format mm-dd-yy: ', 's');
            datenums = str2double(split(datestr_temp, ["_", "-"]));
            datearray = [datenums(3) datenums(1) datenums(2)];
            if datearray(1)<2000
                datearray(1) = datearray(1)+2000;
            end
            dateFormatted = datetime(datearray);
        end
    else
        dateFormatted = datestr;
    end

end