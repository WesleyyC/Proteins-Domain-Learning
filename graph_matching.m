function [ match_matrix, match_score ] = graph_matching( ARG1,ARG2,BLOSUM )
%   GRADUATED_ASSIGN_ALGORITHM is a function that compute the best match
%   matrix with two ARGs

    % set up condition and variable
    % beta is the converging for getting the maximize number
    beta_0 = 0.5;
    beta_f = 10;
    beta_r = 1.075;
    % I control the iteration number for each round
    I_0 = 4;
    I_1 = 30;
    % e control a range
    e_B = 0.5;
    e_C=0.05;    
    % node attriubute compatability weight
    alpha = 1;
    
    
    
    % make sure ARG1 is always the smaller graph
    flip = 0;
    if ARG1.num_nodes>ARG2.num_nodes
        tmp = ARG1;
        ARG1 = ARG2;
        ARG2 = tmp;
        flip=1;
    end
    
    
    
    
    % the size of the real matchin matrix
    A=ARG1.num_nodes;
    I=ARG2.num_nodes;
    real_size = [A,I];
    % the size of the matrix with slacks
    augment_size = real_size+1;
    
    
    
    
    % set up the matrix
    % init a guest m_Head with 1+e
    e=1.5;
    m_Init = rand(augment_size)*e;
    m_Head = m_Init;
    % initial beta to beta_0
    beta = beta_0;
    
    
    
    
    
    
    % pre-calculate the node compatability
    
    % create an function handle for calculating compatibility
    node_compat_handle=@(node1,node2)inner_node_compatibility(node1,node2);
    % calculate the compatibility
    C_n=cellfun(node_compat_handle,repmat(ARG1.nodes',1,I),repmat(ARG2.nodes,A,1));
    % times the alpha weight
    C_n=alpha*C_n;
    
    
    
    
    % pre-calculate the edge compatability
        
    % for low connected rate
    weight_handle = @(edge)edge.weight;
    [i_1,j_1,~] = find(sparse(cellfun(weight_handle,ARG1.edges)));
    [i_2,j_2,~] = find(sparse(cellfun(weight_handle,ARG2.edges)));
    arg1_edges_num = length(i_1);
    arg2_edges_num = length(i_2);
    ARG1_edge_index = mat2cell([i_1,j_1],ones([1,arg1_edges_num]));
    ARG2_edge_index = mat2cell([i_2,j_2],ones([1,arg2_edges_num]));

    C_e = mat2cell(zeros([A*A,I*I]),ones([1,A])*A,ones([1,I])*I);    

    for z = 1:arg1_edges_num
        for y = 1:arg2_edges_num
            index1 = ARG1_edge_index{z};
            index2 = ARG2_edge_index{y};
            C_e{index1(1),index2(1)}(index1(2),index2(2)) = inner_edge_compatibility(ARG1.edges{index1(1),index1(2)},ARG2.edges{index2(1),index2(2)});
        end
    end 

    
    
    
    
    
    % start matching  
    while beta<beta_f   % do A until beta is less than beta_f
        converge_B = 0; % a flag for terminating process B
        I_B = 0;    % counting the iteration of B
        while ~converge_B && I_B <= I_0 % do B until B is converge or iteration exceeds
            old_B=m_Head;   % get the old matrix
            I_B = I_B+1;    % increment the iteration counting
            
            % Build the partial derivative matrix Q
            m_Head_realsize = m_Head(1:A,1:I);
            % sum up the terms for partial differentiation
            sum_fun=@(mat)sum(sum(mat.*m_Head_realsize));
            Q=cellfun(sum_fun,C_e);
            
            %add node attribute
            Q=Q+C_n;
            
            % Normalize Q to avoid NaN/0 produce from exp()
            Q=normr(Q);
            % Now update m_Head!
            m_Head(1:A,1:I)=exp(beta*Q);
            
            converge_C = 0; % a flag for terminating process B
            I_C = 0;    % counting the iteration of C
            %m_One = zeros(size(m_Head));    % a middleware for doing normalization
            while ~converge_C && I_C <= I_1    % Begin C
                I_C=I_C+1;  % increment C
                old_C=m_Head;   % get the m_Head before processing to determine convergence
                
                %normalize the row
                s=sum(m_Head,2);
                n=repmat(s,1,I+1);
                %n(A+1,:)=ones(size(n(A+1,:)));
                m_One=m_Head./n;
                
                % normalize the column
                s=sum(m_One,1);
                n=repmat(s,A+1,1);
                %n(:,I+1)=ones(size(n(:,I+1)));
                m_Head=m_One./n;
                
                % check convergence
                convergeC();
            end
            % check convergence
            convergeB();
        end
        % increment beta
        beta=beta_r*beta;
    end
    
    
    
    
    
    % get the match_matrix in real size
    match_matrix = heuristic(m_Head,A,I);
    match_score = m_Head;

    if(flip)
        match_matrix=match_matrix';
    end
    
    
    
    
    
    
    % inner function
    function [] = convergeC()
        converge_C = abs(sum(sum(m_Head-old_C)))<e_C;
    end

    function [] = convergeB()
        converge_B = abs(sum(sum(m_Head(1:A,1:I)-old_B(1:A,1:I))))<e_B;
    end

    function [c] = inner_node_compatibility(node1, node2)
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

    function [c] = inner_edge_compatibility(edge1, edge2)
        % edge_compatibility function is used to calculate the similarity
        % between edges

        % the score is between [0,1]

        % the higher the score, the more similiarity are there between edge1
        % and edge2

        % this function can be define by the user, but in our case is
        % c(E,e)=1-3|E-e|;

        % assume node1 and node2 are node object

        weight_range = 9;  %update with GenerateProteinARGs

        c = 0;

        %The range of the edge rate
        if ~edge1.trueEdge()||~edge2.trueEdge()
            return;
        else     
            %normalize the score
            c = 1-3*abs(edge1.weight-edge2.weight)/weight_range;
        end

    end
     
end

