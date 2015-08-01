function [c] = node_compatibility(node1, node2, BLOSUM)
    % node_compatibility function is used to calculate the similarity
    % between node1 and node2
    
    % the score is between [0,1]
    
    % the higher the score, the more similiarity are there between node1
    % and node2
    
    % this function can be define by the user, but in our case is
    % c(N,n)=1-3|N-n|;
    
    % assume node1 and node2 are node object
        
    c=0;
    
    if ~node1.hasAtrs()||~node2.hasAtrs()
        return;  % if either of the nodes has NaN attribute, set similarity to 0
    else
        c=BLOSUM(min(node1.atrs,node2.atrs),max(node1.atrs,node2.atrs));
    end
end

