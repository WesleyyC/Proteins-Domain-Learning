function [ proteinARG, protein] = GenerateProteinARG( start_sequence, end_sequence, fileName)
    
    cutoff = 8;

    protein = csvread(fileName);
    
    start_idx = NaN;
    end_idx = NaN;
    
    for i = 1:size(protein,1)
        if start_sequence == protein(i,1) && isnan(start_idx)
            start_idx = i;
        end
        
        if end_sequence==protein(i,1) && isnan(end_idx)
            end_idx= i;
        end
    end
        
    protein = protein(start_idx:end_idx,:);
    
    number_of_AA = size(protein,1);
    proteinEdge = cell(number_of_AA);
    proteinNode = cell([1,number_of_AA]);

    for i = 1:number_of_AA
        proteinNode{i} = protein(i,2);
    end
    
    for i = 1:number_of_AA
        for j = i:number_of_AA
            dist = pdist2(protein(i,3:5),protein(j,3:5));
            dist = dist*(dist<cutoff);
            proteinEdge{i,j} = dist;
            proteinEdge{j,i} = dist;
        end
    end
    
    proteinARG = ARG(proteinEdge,proteinNode);    
end