function allExptData = fillInGene(allExptData)
    % Given allExptData for single gene experiments, obtain the gene for 
    % each image and store in allExptData.
    
    gene_list = {'wt', 'gfp', 'chew', 'lrp', 'flia', 'umod', 'flgm'};
    
    % Iterate over experiments
    for i = 1:length(allExptData)

        % For each experiment, get list of file names
        imfiles = lower(allExptData(i).imfiles);
        % Set up cells to store gene names, IPTGs
        genes = cell(1, length(imfiles));
        iptgs = zeros(1, length(imfiles));
        
        % Erase extension from all of them
        imfiles = erase(imfiles, '.tif');
        imfiles = erase(imfiles, '.jpg');
        
        % Find the gene in the file name
        for j = 1:length(gene_list)
            if any(contains(imfiles, gene_list{j}))
                genes(contains(imfiles, gene_list{j})) = gene_list(j);
            end
        end
        
        % Store found gene
        allExptData(i).genes = genes;
        

    end
end