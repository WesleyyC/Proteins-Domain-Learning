classdef sprMDL < handle & matlab.mixin.Copyable
    %   SpatialPatternMDL is a generalization model of a kind of spatail
    %   pattern providing by a set of sample ARGs   
    
    %   assume null value is 0
    
    properties
        % Important number in model
        number_of_components = NaN;
        number_of_sample = NaN;
        
        % alpha, weight for each model
        weight = NaN;
        
        % the components
        mdl_ARGs={};
        
        % the sample
        sampleARGs = NaN;
        
        % graph matching return
        node_match_scores={};
        node_compatibilities={};
        edge_compatibilities={};
        
        % sample-component matching score
        component_scores = NaN;
        sample_component_matching_probs=NaN;
        
        % threshold score for confirming pattern
        thredshold_score = NaN;
        z_test_mean = 0;
        z_test_sigma = 0;
        
        % debug var
        c_idx = NaN;
        
        % matching prob debugger
        mp = NaN;
        
        
    end
    
    properties (Constant)
        % debug mode
        debug = 0;
        
        % Maximum EM rounds
        iteration_EM = 60;
        
        % Converging epsilon
        e_mdl_converge = 1e-4;
        
        % Node deleting threshold
        % We don't want to delete nodes in early stage
        % so choose such number carefully
        % e_delete_base - e_delete_iter^iter
        e_delete_base = 1;
        e_delete_iter = 0.85;
        e_delete_thredshold=0.9;

        % z_test properties
        z_test_alpha = 0.05;
        z_test_sample_number = 50;
    end
    
    methods
        % Constructor for the class
        function  obj = sprMDL(sampleARGs, number_of_components)
            
            % Throw error if not enough argument
            if nargin < 1
                error "NotEnoughArgument";
            elseif length(sampleARGs)<number_of_components
                error "NotEnoughSample"
            end
            
            % Pass sample ARGs
            obj.sampleARGs=sampleARGs;
            
            % Get the number of components
            obj.number_of_components = number_of_components;
            obj.number_of_sample=length(sampleARGs);
            obj.mdl_ARGs = cell(1,number_of_components);
            
            % Assigning Weight to 1
            obj.weight = ones([1,number_of_components])/number_of_components;
            
            % Randoming pick component from sampleARGs
            idx = randperm(length(sampleARGs)); % we first permutate the index for randomness
            idx = idx(1:number_of_components);  % take what we need
            comp_ARG = sampleARGs(idx);
            obj.c_idx = idx;
            
            % Now convert it to model ARG
            generate_mdl_ARG=@(A)mdl_ARG(A);
            obj.mdl_ARGs=cellfun(generate_mdl_ARG,comp_ARG,'UniformOutput',false); 
            
            % mp
            obj.mp = cell([1,number_of_components]);
            
            % Train the model with the sample
            obj.trainModel();
        end
            
        % Train the model with the sample
        function trainModel(obj)
            % Set up variable
            converge = false;
            iter = 0;
            % EM iteration
            while ~converge && iter<obj.iteration_EM
                % increment the iter
                iter = iter+1;
                % get the old obj before iteration for testing converging
                old_obj = obj.copy();
                % go through one EM iteration
                obj.EM();
                % check converging condition
                converge = sprMDL.mdl_converge(old_obj,obj,obj.e_mdl_converge);
            end
            
            % update the component structure depends on the node frequency
            obj.updateComponentStructure(true);
            obj.EM();
            
            % get the thredshold minimum score
            obj.getThredsholdScore();
        end
        
        % The EM-alogirthem procedure
        function EM(obj)    
            trainRound=tic();
            
            % get the node matching score
            obj.graphMatching(true);
            
            % get the sample-component matching score and probability
            obj.getMatchingProbs();
            
            % update the weight for each component
            obj.updateComponentWeight();
            
            % update the frequency for each component node
            obj.updateComponentNodeFrequency();
            
            % update the atrs for each component node
            obj.updateComponentNodeAAIndex();
            obj.updateComponentNodeVector();
            
            % update the atrs for each component edge
            obj.updateComponentEdgeAtrs();
            
            % update the covariance matrix for each component edge
            obj.updateComponentEdgeCov();
            
            % monitor component strucutre
            obj.updateComponentStructure(false);
            
            toc(trainRound)
        end
        
        % get the graph matching score for each sample-componennt pair
        function graphMatching(obj,train)
            obj.node_match_scores = cell([obj.number_of_sample,obj.number_of_components]);
            obj.node_compatibilities = cell(size(obj.node_match_scores));
            obj.edge_compatibilities = cell(size(obj.node_match_scores));
            
            for i=1:obj.number_of_sample
                for j = 1:obj.number_of_components
                    [node_match_score,node_compatibility,edge_compatibility] = graph_matching(obj.sampleARGs{i},obj.mdl_ARGs{j},train);
                    obj.node_match_scores{i,j}=node_match_score;
                    obj.node_compatibilities{i,j}=node_compatibility;
                    obj.edge_compatibilities{i,j}=edge_compatibility;
                end
            end
            
            if obj.debug
                obj.debugMatch()
            end
        end
        
        % get the sample-component matching score and probability
        function getMatchingProbs(obj)
            % Use the graphMatching() data calculate the sample-componennt
            % matching score & probability
            handle=@(node_match_score,node_compatibility,edge_compatibility)sprMDL.component_score(node_match_score,node_compatibility,edge_compatibility);
            obj.component_scores=cellfun(handle,obj.node_match_scores,obj.node_compatibilities,obj.edge_compatibilities);
            s=sum(obj.component_scores,2);
            n=repmat(s,1,obj.number_of_components);
            obj. sample_component_matching_probs=obj.component_scores./n;
        end
        
        % update the weight for each component
        function updateComponentWeight(obj)
            % sum up the samples-component probabilities and divide by
            % sample number
            obj.weight = sum(obj.sample_component_matching_probs)/obj.number_of_sample;
        end
        
         
        % update the frequency for each component node
        function updateComponentNodeFrequency(obj)
            % for each component
            for i = 1:obj.number_of_components
                frequency=0;
                sample_node_sum = 0;
                % for each sample
                for j = 1:obj.number_of_sample
                    % we calculate the component node frequency in this
                    % specific sample
                    current_freq =sum(obj.node_match_scores{j,i});
                    frequency=frequency+current_freq*obj.sample_component_matching_probs(j,i);
                    sample_node_sum=sample_node_sum+obj.sampleARGs{j}.num_nodes*obj.sample_component_matching_probs(j,i);
                end
                % update the frequency for all the nodes in the model in
                % the same time
                obj.mdl_ARGs{i}.nodes_freq = frequency/sample_node_sum;
            end    
        end
        
        % update the atrs for each component node
        function updateComponentNodeAAIndex(obj)
            % for each component
            for h = 1:obj.number_of_components
                % for each node
                for n = 1:obj.mdl_ARGs{h}.num_nodes
                    atrs = 0;
                    denominator=0;
                    % we go over the sample
                    for i = 1:obj.number_of_sample
                        [v, a] = max(obj.node_match_scores{i,h}(:,n));
                        if v > denominator
                            atrs = obj.sampleARGs{i}.nodes_vector(a,:);
                        end
                    end
                    % udpate the value
                    obj.mdl_ARGs{h}.nodes_aa_index(n,:)=atrs;
                end
            end       
        end
        
        % update the covariance matrix for each component node
        function updateComponentNodeVector(obj)
            % compute local distribution
            % for each component
            for h = 1:obj.number_of_components
                % for each node
                for n = 1:obj.mdl_ARGs{h}.num_nodes
                    if any(obj.mdl_ARGs{h}.nodes_vector(n,:))
                        vec = 0;
                        % we go over the sample
                        for i = 1:obj.number_of_sample
                            current_vec = zeros(1,20);
                            for j = 1:length(obj.node_match_scores{i,h}(:,n))
                                current_vec(obj.sampleARGs{i}.nodes_vector(j,:))=current_vec(obj.sampleARGs{i}.nodes_vector(j,:))+obj.node_match_scores{i,h}(j,n);
                            end
                            vec  = vec + current_vec*obj.sample_component_matching_probs(i,h);
                        end
                        % udpate the value
                        new_vec = normr(vec).*normr(vec);
                        obj.mdl_ARGs{h}.nodes_vector(n,:) = new_vec;
                    end
                end
            end    
            
            % compute global distribution
            global_vec = zeros(20,20);
            for h = 1:obj.number_of_components
                for n = 1:obj.mdl_ARGs{h}.num_nodes
                    global_vec(obj.mdl_ARGs{h}.nodes_aa_index(n,:),:)=global_vec(obj.mdl_ARGs{h}.nodes_aa_index(n,:),:)+obj.mdl_ARGs{h}.nodes_vector(n,:)*obj.weight(h);
                end
            end
            global_vec = normr(global_vec).*normr(global_vec);
            
            % update local distribution
            for h = 1:obj.number_of_components
                for n = 1:obj.mdl_ARGs{h}.num_nodes
                    new_vec = obj.mdl_ARGs{h}.nodes_vector(n,:)./global_vec(obj.mdl_ARGs{h}.nodes_aa_index(n,:));
                    new_vec = log(new_vec);
                    new_vec = normr(new_vec).*normr(new_vec);
                    obj.mdl_ARGs{h}.nodes_vector(n,:) = new_vec;
                end
            end
                
        end
        
        % update the atrs for each component edge
        function updateComponentEdgeAtrs(obj)
            %for each component
            for h = 1:obj.number_of_components
                %for each edge 
                for o = 1:obj.mdl_ARGs{h}.num_nodes
                    for t = o+1:obj.mdl_ARGs{h}.num_nodes
                        if any(obj.mdl_ARGs{h}.edges_matrix(o,t,:))
                            atrs = 0;
                            denominator=0;
                            %for each sample
                            for i = 1:obj.number_of_sample
                                base = repmat(obj.node_match_scores{i,h}(:,o)*obj.node_match_scores{i,h}(:,t)',1,1,size(obj.sampleARGs{i}.edges_matrix,3));
                                current_sample_atrs = obj.sampleARGs{i}.edges_matrix.*base;
                                current_sample_denominator = any(obj.sampleARGs{i}.edges_matrix,3).*(obj.node_match_scores{i,h}(:,o)*obj.node_match_scores{i,h}(:,t)');
                                
                                current_sample_atrs = sum(sum(current_sample_atrs));
                                current_sample_denominator = sum(sum(current_sample_denominator));
                                
                                atrs=atrs+current_sample_atrs*obj.sample_component_matching_probs(i,h);
                                denominator = denominator + current_sample_denominator*obj.sample_component_matching_probs(i,h);
                            end
                            % update the value
                            new_atr = atrs/denominator;
                            obj.mdl_ARGs{h}.edges_matrix(o,t,:) = new_atr;
                            obj.mdl_ARGs{h}.edges_matrix(t,o,:) = new_atr;
                        end
                    end
                end
            end                     
        end
        
        % update the covariance matrix for each component edge
        function updateComponentEdgeCov(obj)
            get_cov_time=@(r) reshape(reshape(r,1,[])'*reshape(r,1,[]),1,length(r)^2);
            %for each component
            for h = 1:obj.number_of_components
                %for each edge 
                for o = 1:obj.mdl_ARGs{h}.num_nodes
                    for t = o+1:obj.mdl_ARGs{h}.num_nodes
                        if any(obj.mdl_ARGs{h}.edges_matrix(o,t,:))
                            cov = 0;
                            denominator=0;
                            %for each sample
                            for i = 1:obj.number_of_sample
                                z_atrs = obj.sampleARGs{i}.edges_matrix-repmat(obj.mdl_ARGs{h}.edges_matrix(o,t,:),size(obj.sampleARGs{i}.edges_matrix,1),size(obj.sampleARGs{i}.edges_matrix,2),1);
                                z_atrs = repmat(any(obj.sampleARGs{i}.edges_matrix,3),1,1,size(z_atrs,3)).*z_atrs;
                                base = obj.node_match_scores{i,h}(:,o)*obj.node_match_scores{i,h}(:,t)';
                                zz_atrs = NaN(size(z_atrs,1),size(z_atrs,2),size(z_atrs,3)^2);
                                for x = 1:size(z_atrs,1)
                                    for y = 1:size(z_atrs,2)
                                        zz_atrs(x,y,:) = get_cov_time(z_atrs(x,y,:));
                                    end
                                end
                                current_sample_cov = zz_atrs.*repmat(base,1,1,size(zz_atrs,3));
                                current_sample_denominator = any(obj.sampleARGs{i}.edges_matrix,3).*base;
                                
                                current_sample_cov = sum(sum(current_sample_cov));
                                current_sample_denominator = sum(sum(current_sample_denominator));
                                
                                cov=cov+current_sample_cov*obj.sample_component_matching_probs(i,h);
                                denominator = denominator + current_sample_denominator*obj.sample_component_matching_probs(i,h);
                            end
                            % update the value
                            new_cov = cov/denominator;
                            new_cov = sprMDL.normalize_cov(new_cov);
                            obj.mdl_ARGs{h}.edges_cov(o,t,:) = new_cov;
                            obj.mdl_ARGs{h}.edges_cov(t,o,:) = new_cov;
                        end
                    end
                end
            end     
        end
        
        % update the component structure depends on the node frequency
        function updateComponentStructure(obj,modify)
            % for each component
            for w = 1:obj.number_of_components
                av_frequency=0;
                prob_sum = 0;
                % for each sample
                for j = 1:obj.number_of_sample
                    % we calculate the component node frequency in this
                    % specific sample
                    current_freq =sum(obj.node_match_scores{j,w});
                    av_frequency=av_frequency+current_freq*obj.sample_component_matching_probs(j,w);
                    prob_sum=prob_sum+obj.sample_component_matching_probs(j,w);
                end
                % get the average matching probability
                av_matching_prob = av_frequency/prob_sum;
                obj.mp{w}(end+1,1:length(av_matching_prob)) = av_matching_prob;
                if modify
                    % delet the node that is less tha the threshold 1-e^iter
                    deleteIdx = av_matching_prob < obj.e_delete_thredshold;
                    deleteIdx(end)=0; % the null node will always be remained
                    obj.mdl_ARGs{w}.modifyStructure(deleteIdx);
                end
            end
        end
                
        % Get the thredshold score for confimrming pattern
        function getThredsholdScore(obj)                
            random_sample_number = obj.z_test_sample_number;
            random_score = zeros([1,length(random_sample_number)]);

            count = 1;

            for i = 1:length(obj.sampleARGs)
                for j = 1:random_sample_number/length(obj.sampleARGs)
                    tmpARG = obj.sampleARGs{i}.copy;
                    n_idx = randperm(tmpARG.num_nodes);
                    e_idx = randperm(tmpARG.num_nodes);
                    tmpARG.nodes = tmpARG.nodes(n_idx);
                    tmpARG.nodes_vector = tmpARG.nodes_vector(n_idx, :);
                    tmpARG.edges = tmpARG.edges(e_idx, e_idx);
                    tmpARG.edges_matrix = tmpARG.edges_matrix(e_idx, e_idx, :);
                    score = obj.scorePattern(tmpARG);
                    random_score(count) = score;
                    count = count + 1;
                end
            end
            
            obj.z_test_mean = mean(random_score);
            obj.z_test_sigma = std(random_score);
            obj.thredshold_score = obj.z_test_mean + 3 * obj.z_test_sigma;
            
        end
        
        % Detect if a ARG has the same pattern
        function [tf, score] = checkPattern(obj, ARG)
            score = obj.scorePattern(ARG);
            tf = ztest(score, obj.z_test_mean,obj.z_test_sigma,'Alpha',obj.z_test_alpha, 'Tail','Right');
        end
        
        % Detect if a ARG has the same pattern
        function score = scorePattern(obj, ARG)
            score = 0;
            match = cell([1,obj.number_of_components]);
            compatibility = cell([1,obj.number_of_components]);
            for i = 1:obj.number_of_components
               [node_match_score,node_compatibility,edge_compatibility]=graph_matching(ARG,obj.mdl_ARGs{i},false);
               match{i} = node_match_score;
               compatibility{i} = node_compatibility;
               score = score + ...
                    sprMDL.component_score(node_match_score,node_compatibility,edge_compatibility) * obj.weight(i);
            end
            if obj.debug
                figure;
                f_count = 1;
                fig_x = 2;
                fig_y = obj.number_of_components;

                % match score
                for i = 1:obj.number_of_components
                    subplot(fig_x, fig_y, f_count);
                    imshow(match{i})
                    f_count = f_count + 1;
                end
                
                for i = 1:obj.number_of_components
                    subplot(fig_x, fig_y, f_count);
                    imshow(compatibility{i})
                    f_count = f_count + 1;
                end
            end
        end
        
        function [] = debugMatch(obj)
            figure;
            f_count = 1;
            fig_x = obj.number_of_components*2;
            fig_y = obj.number_of_sample;
            
            % match score
            for i = 1:obj.number_of_sample
                for j = 1:obj.number_of_components
                    subplot(fig_x, fig_y, f_count);
                    imshow(obj.node_match_scores{i,j})
                    f_count = f_count + 1;
                end
            end
            
            for i = 1:obj.number_of_sample
                for j = 1:obj.number_of_components
                    subplot(fig_x, fig_y, f_count);
                    imshow(obj.node_compatibilities{i,j})
                    f_count = f_count + 1;
                end
            end
            
            figure;
            fig_x = 1;
            fig_y = obj.number_of_components;
            for i = 1:obj.number_of_components
                subplot(fig_x,fig_y,i);
                imshow(obj.mp{i})
            end
            
        end
    end
    
    methods(Static)
        % scoring for each sample-component pair
        function score = component_score(node_match_score,node_compatibility,edge_compatibility)
            
            % get the size
            [A, I] = size(node_match_score);
           
            % no null node
            node_match_score(:,end)=0;
            
            % edge score
            Q = zeros(A, I);
            for a = 1:A
                for i = 1:I
                    Q(a,i)=sum(sum(edge_compatibility(((a-1)*A+1):((a-1)*A+A),((i-1)*I+1):((i-1)*I+I)).*node_match_score));
                end
            end
            
            % node score
            Q = Q + node_compatibility;
            
            % consider probability
            Q = Q .* node_match_score;
            
            % sum up the score
            score = sum(sum(Q));
        end
        
        % mdl_converge judges if the model is converged
        function converge = mdl_converge( old_mdl, new_mdl, e )
            if length(old_mdl.weight)==length(new_mdl.weight)
                diff = sum(abs(old_mdl.weight-new_mdl.weight))/length(new_mdl.weight);
                converge = diff<e;
            else
                converge = false;
            end
        end
        
        % normalize_cov
        function M = normalize_cov(M)
            M = reshape(M,sqrt(length(M)),sqrt(length(M)));
            M = normr(M).*normr(M);
            M = normc(M).*normc(M);
            M = normr(M).*normr(M);
            M = reshape(M,1,[]);
        end
    end
end

