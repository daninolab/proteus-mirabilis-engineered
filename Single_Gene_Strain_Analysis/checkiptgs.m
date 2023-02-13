function test = checkiptgs(allExptData)
% use this to make a cell with which to check all the iptgs and file names
% together easily

test = cell(length(allExptData), 1);
for i = 1:length(allExptData)
    temp = cell(length(allExptData(i).imfiles), 2);
    for j = 1:length(allExptData(i).imfiles)
        temp{j, 1} = allExptData(i).imfiles{j};
        temp{j, 2} = allExptData(i).iptgs(j);
    end
    test{i} = temp;
end
    

end