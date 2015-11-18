function [ match_matrix, match_score ] = graph_matching( ARG1,ARG2,BLOSUM )
%   GRADUATED_ASSIGN_ALGORITHM is a function that compute the best match
%   matrix with two ARGs

    % parallel computing flag
    pFlag = 1;

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
    
    weight_range = 9;
    
    tic()
    
    % pre-calculate the edge compatability
    C_e = sparse(A*A,I*I);  
    
    if pFlag
        edges_ARG1 = sparse(flattern_matrix(ARG1.edges));
        edges_ARG2 = sparse(flattern_matrix(ARG2.edges));
        parfor p = 1:A*A
            C_e(p,:)=sparse((1-3*abs(edges_ARG2-edges_ARG1(p))/weight_range).*(edges_ARG2>0)*(edges_ARG1(p)>0));        
        end
    else
        for p = 1:A*I
            fill_Ce(floor((p-1)/I)+1,p-(floor((p-1)/I))*I);
        end
    end
    toc()
    
    
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
            sum_fun=@(a,i)sum(sum(full(C_e(((a-1)*A+1):((a-1)*A+A),((i-1)*I+1):((i-1)*I+I))).*m_Head_realsize));
            Q=cellfun(sum_fun,num2cell(repmat((1:A)',1,I)),num2cell(repmat((1:I),A,1)));
            
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
    
    
    
    
    % INNER FUNCTION SECTION
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

    function [] = fill_Ce(a,i)
        edge_compat_handle=@(edge1,edge2)(1-3*abs(edge1.weight-edge2.weight)/weight_range).*(edge1.weight>0)*(edge2.weight>0);
        C_e(((a-1)*A+1):((a-1)*A+A),((i-1)*I+1):((i-1)*I+I))=cellfun(edge_compat_handle,...
            repmat(ARG1.edges(a,:)',1,I),...
            repmat(ARG2.edges(i,:),A,1));
    end

    function [flat] = flattern_matrix(edges)
        len = length(edges);
        flat = zeros (1, len*len);
        for flat_p = 1:len*len
            flat(flat_p)=edges{floor((flat_p-1)/len)+1,flat_p-(floor((flat_p-1)/len))*len}.weight;
        end
    end
end

