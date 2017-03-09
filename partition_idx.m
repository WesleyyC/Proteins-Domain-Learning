function idx = partition_idx (portion, num_nodes)
    idx = [];

    idx_endpoint = zeros([1,portion]);

    idx_endpoint(1)=1;
    for i = 2:length(idx_endpoint)
        idx_endpoint(i) = round((i-1)*num_nodes/portion);
    end
    
    for i = length(idx_endpoint):-1:1
        start_idx = idx_endpoint(i);
        if i == length(idx_endpoint)
            end_idx = num_nodes;
        else
            end_idx = idx_endpoint(i+1)-1;
        end
        
        idx = [idx,start_idx:end_idx];
    end
end