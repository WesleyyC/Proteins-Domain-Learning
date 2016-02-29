 function [ match_matrix, C_n, C_e] = graph_matching( ARG1,ARG2,BLOSUM )
%   GRADUATED_ASSIGN_ALGORITHM is a function that compute the best match
%   matrix with two ARGs

    % parallel computing flag
    pFlag=1;

    % set up condition and variable
    % beta is the converging for getting the maximize number
    beta_0 = 0.5;
    beta_f = 10;
%     beta_f = 50;
    beta_r = 1.075;
    % I control the iteration number for each round
    I_0 = 4;
%     I_0 = 50;
    I_1 = 30;
    % e control a range
    e_B = 0.5;
    e_C=0.05;    
    % node attriubute compatability weight
    alpha = 0.1;
    
    flip = 0;
    if ~isa(ARG2,'mdl_ARG')
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
    augment_size = real_size+[1,0];
    
    % set up the matrix
    % init a guest m_Head with 1+e
    e=1.5;
    m_Init = rand(augment_size)*e;
    m_Head = m_Init;
    % initial beta to beta_0
    beta = beta_0;
%     m_Head(end,:)=e/2;
%     m_Head(:,end)=e/2;
    
    % pre-calculate the node compatability
    % create an function handle for calculating compatibility
    node_compat_handle=@(node1,node2)inner_node_compatibility(node1,node2);
    % calculate the compatibility
    C_n=cellfun(node_compat_handle,repmat(ARG1.nodes',1,I),repmat(ARG2.nodes,A,1));
    % times the alpha weight
    C_n=alpha*C_n;
    
    % pre-calculate the edge compatability
    C_e = zeros(A*A,I*I);  
    
    if pFlag
        edges_atrs = flattern_matrix_weight(ARG1.edges);
        mdl_edges_atrs = flattern_matrix_weight(ARG2.edges);
        mdl_edges_cov = flattern_matrix_cov(ARG2.edges);
        mdl_edges_cov_inv = flattern_matrix_cov_inv(ARG2.edges);
        edge_num_atrs = 1;
        parfor p = 1:A*A
            % because edge atr is a single number, we can do some
            % modification to our orignal formula
            C_e(p,:)=(exp(-0.5*(edges_atrs(p)-mdl_edges_atrs).*mdl_edges_cov_inv.*(edges_atrs(p)-mdl_edges_atrs))./...
                ((2*pi)^(edge_num_atrs/2)*sqrt(mdl_edges_cov)))...
                .*(edges_atrs(p)>0).*(mdl_edges_atrs>0);      
        end
        
    else
        for p = 1:A*I
            fill_Ce(floor((p-1)/I)+1,p-(floor((p-1)/I))*I);
        end
    end
    
%     figure()
%     imshow(normr(C_n),'InitialMagnification',2000)
%     figure()

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
%             m_Head(1:A,1:I)=exp(beta*Q);
            m_Head(1:A,1:I)=normc(exp(beta*Q));
            
            converge_C = 0; % a flag for terminating process B
            I_C = 0;    % counting the iteration of C
            %m_One = zeros(size(m_Head));    % a middleware for doing normalization
            
            while ~converge_C && I_C <= I_1    % Begin C
                
                I_C=I_C+1;  % increment C
                old_C=m_Head;   % get the m_Head before processing to determine convergence
                
                %normalize the row
                s=sum(m_Head,2);
                n=repmat(s,1,I);
                n(end,:)=ones(size(n(end,:)));
                m_One=m_Head./n;
                
                % normalize the column
                s=sum(m_One,1);
                n=repmat(s,A+1,1);
                n(:,end)=ones(size(n(:,end)));
                m_Head=m_One./n;
                
                % check convergence
               
                convergeC();
            end
            
%             imshow(m_Head,'InitialMagnification',2000);
            % check convergence
            convergeB();
            
        end
        
        % increment beta
        
        beta=beta_r*beta;
    end
%     imshow(m_Head,'InitialMagnification',2000);

    % get the match_matrix in real size
    match_matrix = heuristic(m_Head,A,I);
    
    if(flip)
        match_matrix=match_matrix';
        C_n = C_n';
        C_e = C_e';
        tmp=A;
        A=I;
        I=tmp;
    end
    C_n = C_n/alpha;
    C_e = mat2cell(full(C_e),ones([1,A])*A,ones([1,I])*I);
    
    
    
    
    % INNER FUNCTION SECTION
    function [] = convergeC()
        converge_C = abs(sum(sum(m_Head-old_C)))<e_C;
    end

    function [] = convergeB()
        converge_B = abs(sum(sum(m_Head(1:A,1:I)-old_B(1:A,1:I))))<e_B;
    end

    function [c] = inner_node_compatibility(node1, node2)

        c=0;

        if ~node1.hasAtrs()||~node2.hasAtrs()
            return;  % if either of the nodes has NaN attribute, set similarity to 0
        else
            c=node2.getAtrs()*BLOSUM*node1.getAtrs()';
            % same thing
            %c=sum(sum(BLOSUM.*(node2.atrs'*node1.atrs)));
        end
    end

    function [c] = inner_edge_compatibility(edge, mdl_edge)
    % node_compatibility function is used to calculate the similarity
    % between node1 and node2
    
        c=0;

        if ~edge.trueEdge()||~mdl_edge.trueEdge()
            return;  % if either of the edge does not exist or has NaN attribute
        else

            % get number of attributes
            num_atrs = mdl_edge.numberOfAtrs();

            % get the mean of attributes
            edge_atrs = edge.weight;
            mdl_edge_atrs = mdl_edge.weight;
            % get the covariance matrix of model node
            mdl_edge_cov = mdl_edge.cov;
            mdl_edge_cov_inv = mdl_edge.cov_inv;

            % calculate the score
            c=exp(-0.5*(edge_atrs-mdl_edge_atrs)*mdl_edge_cov_inv*(edge_atrs-mdl_edge_atrs)')/...
                ((2*pi)^(num_atrs/2)*sqrt(det(mdl_edge_cov)));
        end
    end

    function [] = fill_Ce(a,i)
        edge_compat_handle=@(edge1,edge2)inner_edge_compatibility(edge1,edge2);
        C_e(((a-1)*A+1):((a-1)*A+A),((i-1)*I+1):((i-1)*I+I))=cellfun(edge_compat_handle,...
            repmat(ARG1.edges(a,:)',1,I),...
            repmat(ARG2.edges(i,:),A,1));
    end

    function [flat] = flattern_matrix_weight(edges)
        len = length(edges);
        flat = zeros (1, len*len);
        for flat_p = 1:len*len
            flat(flat_p)=edges{floor((flat_p-1)/len)+1,flat_p-(floor((flat_p-1)/len))*len}.getAtrs();
        end
    end

    function [flat] = flattern_matrix_cov(edges)
        len = length(edges);
        flat = zeros (1, len*len);
        for flat_p = 1:len*len
            flat(flat_p)=edges{floor((flat_p-1)/len)+1,flat_p-(floor((flat_p-1)/len))*len}.getCov();
        end
    end

    function [flat] = flattern_matrix_cov_inv(edges)
        len = length(edges);
        flat = zeros (1, len*len);
        for flat_p = 1:len*len
            flat(flat_p)=edges{floor((flat_p-1)/len)+1,flat_p-(floor((flat_p-1)/len))*len}.getCovInv();
        end
    end
end

