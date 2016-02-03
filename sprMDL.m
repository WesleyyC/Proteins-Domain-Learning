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
        
        % BLOSUM Matrix
        BLOSUM = NaN;
        
        
        %exp
        counter = 1
        base_counter = 0;
        base = zeros(1,2)
        distance = zeros(1,2)
        node_ab_match = zeros(1,2)
        node_cd_match = zeros(1,2)
        sample_match = zeros(1,2)
        edge_match = zeros(1,2)
        node_deleted = 0;
        
        
    end
    
    properties (Constant)
        % Maximum EM rounds
        iteration_EM = 30;
        % Converging epsilon
        e_mdl_converge = 1e-4;
        % Node deleting threshold
        % We don't want to delete nodes in early stage
        % so choose such number carefully
        % e_delete_base - e_delete_iter^iter
        e_delete_base = 1;
        e_delete_iter = 0.8;
    end
    
    methods
        % Constructor for the class
        function  obj = sprMDL(sampleARGs,number_of_components)
            
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
            % Now convert it to model ARG
            obj.mdl_ARGs=cellfun(@mdl_ARG,comp_ARG,'UniformOutput',false); 
            
            % Build BLOSUM matrix
            BLOSUM_Sigma = 2;
            C=[9,-1,-1,-3,0,-3,-3,-3,-4,-3,-3,-3,-3,-1,-1,-1,-1,-2,-2,-2];
            S=[-1,4,1,-1,1,0,1,0,0,0,-1,-1,0,-1,-2,-2,-2,-2,-2,-3];
            T=[-1,1,4,1,-1,1,0,1,0,0,0,-1,0,-1,-2,-2,-2,-2,-2,-3];
            P=[-3,-1,1,7,-1,-2,-1,-1,-1,-1,-2,-2,-1,-2,-3,-3,-2,-4,-3,-4];
            A=[0,1,-1,-1,4,0,-1,-2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-2,-3];
            G=[-3,0,1,-2,0,6,-2,-1,-2,-2,-2,-2,-2,-3,-4,-4,0,-3,-3,-2];
            N=[-3,1,0,-2,-2,0,6,1,0,0,-1,0,0,-2,-3,-3,-3,-3,-2,-4];
            D=[-3,0,1,-1,-2,-1,1,6,2,0,-1,-2,-1,-3,-3,-4,-3,-3,-3,-4];
            E=[-4,0,0,-1,-1,-2,0,2,5,2,0,0,1,-2,-3,-3,-3,-3,-2,-3];
            Q=[-3,0,0,-1,-1,-2,0,0,2,5,0,1,1,0,-3,-2,-2,-3,-1,-2];
            H=[-3,-1,0,-2,-2,-2,1,1,0,0,8,0,-1,-2,-3,-3,-2,-1,2,-2];
            R=[-3,-1,-1,-2,-1,-2,0,-2,0,1,0,5,2,-1,-3,-2,-3,-3,-2,-3];
            K=[-3,0,0,-1,-1,-2,0,-1,1,1,-1,2,5,-1,-3,-2,-3,-3,-2,-3];
            M=[-1,-1,-1,-2,-1,-3,-2,-3,-2,0,-2,-1,-1,5,1,2,-2,0,-1,-1];
            I=[-1,-2,-2,-3,-1,-4,-3,-3,-3,-3,-3,-3,-3,1,4,2,1,0,-1,-3];
            L=[-1,-2,-2,-3,-1,-4,-3,-4,-3,-2,-3,-2,-2,2,2,4,3,0,-1,-2];
            V=[-1,-2,-2,-2,0,-3,-3,-3,-2,-2,-3,-3,-2,1,3,1,4,-1,-1,-3];
            F=[-2,-2,-2,-4,-2,-3,-3,-3,-3,-3,-1,-3,-3,0,0,0,-1,6,3,1];
            Y=[-2,-2,-2,-3,-2,-3,-2,-3,-2,-1,2,-2,-2,-1,-1,-1,-1,3,7,2];
            W=[-2,-3,-3,-4,-3,-2,-4,-4,-3,-2,-2,-3,-3,-1,-3,-2,-3,1,2,11];
            BLOSUM=[C',S',T',P',A',G',N',D',E',Q',H',R',K',M',I',L',V',F',Y',W'];
            BLOSUM=exp(BLOSUM/BLOSUM_Sigma);
            obj.BLOSUM = BLOSUM;
            
            % Train the model with the sample
            obj.trainModel();
        end
        
        % if has the same pattern, return a new ARG which is the partial
        % ARG from the original ARG that has the similar pattern
        function pattern = getSamePattern(obj,testARG)
            % if the original ARG has the similar pattern
            if obj.checkSamePattern(testARG)
                % get the sample-components weight
                sample_component_weight = zeros([1,obj.number_of_components]);
                for i = 1:obj.number_of_components
                    [node_match_score,node_compatibility,edge_compatibility]=graph_matching(testARG,obj.mdl_ARGs{i},obj.BLOSUM);
                    sample_component_weight(i)=...
                        sprMDL.component_score(node_match_score,node_compatibility,edge_compatibility,false) * obj.weight(i);
                end
                % normalize to one
                sample_component_weight=sample_component_weight/sum(sample_component_weight);
                
                % ge each nodes score
                nodes_score = 0;
                for i = 1:obj.number_of_components
                    node_match_score=graph_matching(testARG,obj.mdl_ARGs{i},obj.BLOSUM);
                    node_match_score = node_match_score.*repmat(obj.mdl_ARGs{i}.getNodeFrequency(),testARG.num_nodes,1);
                    node_match_score = sum(node_match_score,2)';
                    nodes_score = nodes_score+node_match_score*obj.weight(i)*sample_component_weight(i);
                end
                
                % setting up selection thredshol
                getSampleNodeHandle = @(sample)sample.num_nodes;
                totalSampleNode = sum(cellfun(getSampleNodeHandle,obj.sampleARGs));
                e_node_selection = 0.95*obj.number_of_sample/(obj.number_of_components*totalSampleNode);
                
                % select node
                idx = find(nodes_score >= e_node_selection);
                % get the building structure
                patternM = testARG.matrix(idx,idx);
                patternAtrs = testARG.atrs_vector(idx);
                % build a new ARG for represented-pattern
                pattern = ARG(patternM,patternAtrs);
                
            else
                pattern = NaN;                
            end
                
        end
        
        % Detect if a ARG has the same pattern
        function tf = checkSamePattern(obj, ARG)
            scores = zeros(1,obj.number_of_components);
            parfor i = 1:obj.number_of_components
               [node_match_score,node_compatibility,edge_compatibility]=graph_matching(ARG,obj.mdl_ARGs{i},obj.BLOSUM);
               scores(i) = sprMDL.component_score(node_match_score,node_compatibility,edge_compatibility,false) * obj.weight(i);
            end
            tf = sum(scores)>=obj.thredshold_score;
        end
        
        % showing the pattern that this model summarized
        function pattern = summarizedPattern(obj)
            % get the id of the most weighted model
%             component_prob = sum(obj.sample_component_matching_probs).*obj.weight;
%             component_prob = component_prob/sum(component_prob);
            [~,idx]=max(obj.weight);
            % get the model, and in case there are multiple ARG, we only
            % takes the first elemtn
            representMDL = obj.mdl_ARGs{idx(1)};
            % get the representation
            struct = representMDL.showARG();
            pattern = ARG(struct.M(1:end-1,1:end-1),struct.Na(1:end-1));
            view(struct.bg);
        end
        
        % Get the thredshold score for confimrming pattern
        function getThredsholdScore(obj)
            handle=@(node_match_score,node_compatibility,edge_compatibility)sprMDL.component_score(node_match_score,node_compatibility,edge_compatibility,false);
            obj.component_scores=cellfun(handle,obj.node_match_scores,obj.node_compatibilities,obj.edge_compatibilities);
            scores = obj.component_scores.*repmat(obj.weight,obj.number_of_sample,1);
            scores = sum(scores,2);
            obj.thredshold_score = min(scores);
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
                obj.EM(iter);
                % check converging condition
                converge = mdl_converge(old_obj,obj,obj.e_mdl_converge);
            end
            % get the thredshold minimum score
            obj.getThredsholdScore();
        end
        
        % The EM-alogirthem procedure
        function EM(obj,iter)   
            EMRoundStart=tic();
            % get the node matching score
            obj.graphMatching();
            % get the sample-component matching score and probability
            obj.getMatchingProbs();
            % update the weight for each component
            obj.updateComponentWeight();
            % update the frequency for each component node
            obj.updateComponentNodeFrequency();
            % update the atrs for each component node
            obj.updateComponentNodeAtrs();
            % update the atrs for each component edge
            obj.updateComponentEdgeAtrs();
            % update the covariance matrix for each component edge
            obj.updateComponentEdgeCov();
            % update the component structure depends on the node frequency
            obj.updateComponentStructure(iter);
            obj.base_counter = obj.base_counter+1;
            toc(EMRoundStart);
        end
        
        % update the component structure depends on the node frequency
        function updateComponentStructure(obj,iter)
            % for each component
            for w = 1:obj.number_of_components
                av_frequency=zeros(obj.number_of_sample,obj.mdl_ARGs{w}.num_nodes);
                prob_sum = zeros(1,obj.number_of_sample);
                % for each sample
                parfor j = 1:obj.number_of_sample
                    % we calculate the component node frequency in this
                    % specific sample
                    current_freq =sum(obj.node_match_scores{j,w});
                    av_frequency(j,:)=current_freq*obj.sample_component_matching_probs(j,w);
                    prob_sum(j)=obj.sample_component_matching_probs(j,w);
                end
                av_frequency=sum(av_frequency);
                prob_sum=sum(prob_sum);
                % get the average matching probability
                av_matching_prob = av_frequency/prob_sum;
                % delet the node that is less tha the threshold 1-e^iter
                deleteIdx = av_matching_prob < obj.e_delete_base-obj.e_delete_iter^iter;
                deleteIdx(end)=0; % the null node will always be remained
                obj.mdl_ARGs{w}.modifyStructure(deleteIdx);
                
                %exp
                obj.node_deleted=obj.node_deleted+length(find(deleteIdx));
            end
        end
        
        % Edge updating can be done using a list of component instead of traversing
        % all the nodes with SIX nested for loop. But this version is easier to code so we will
        % start it from here.
        
        % update the atrs for each component edge
        function updateComponentEdgeAtrs(obj)
            %for each component
            for h = 1:obj.number_of_components
                %for each edge 
                for o = 1:obj.mdl_ARGs{h}.num_nodes
                    for t = o+1:obj.mdl_ARGs{h}.num_nodes
                        if any(obj.mdl_ARGs{h}.edges{o,t}.getAtrs())
                            atrs = 0;
                            denominator=0;
                            %for each sample
                            for i = 1:obj.number_of_sample
                                current_sample_atrs = 0;
                                current_sample_denominator = 0;
                                %for each edge in sample
                                for c =  1:obj.sampleARGs{i}.num_nodes
                                    for d =  1:obj.sampleARGs{i}.num_nodes
                                        if any(obj.sampleARGs{i}.edges{c,d}.getAtrs())
                                            current_sample_atrs=current_sample_atrs+...
                                                obj.sampleARGs{i}.edges{c,d}.getAtrs()*obj.node_match_scores{i,h}(c,o)*obj.node_match_scores{i,h}(d,t);
                                            current_sample_denominator=current_sample_denominator+obj.node_match_scores{i,h}(c,o)*obj.node_match_scores{i,h}(d,t);
                                        end
                                    end
                                end
                                atrs=atrs+current_sample_atrs*obj.sample_component_matching_probs(i,h);
                                denominator = denominator + current_sample_denominator*obj.sample_component_matching_probs(i,h);
                            end
                            % update the value
                            obj.mdl_ARGs{h}.edges{o,t}.updateAtrs(atrs/denominator);
                            obj.mdl_ARGs{h}.edges{t,o}.updateAtrs(atrs/denominator);
                        end
                    end
                end
            end                     
        end
        
%         function updateComponentEdgeAtrs(obj)
%             update_result = cell(1,obj.number_of_components);
%             
%             %for each component
%             parfor h = 1:obj.number_of_components
%                 sum_atrs = 0;
%                 sum_denominator = 0;
%                 %for each sample
%                 for i = 1:obj.number_of_sample              
%                     sum_atrs=sum_atrs+obj.node_match_scores{i,h}'*obj.sampleARGs{i}.edges_matrix*obj.node_match_scores{i,h}*obj.sample_component_matching_probs(i,h);
%                     sum_denominator=sum_denominator+obj.node_match_scores{i,h}'*obj.node_match_scores{i,h}*obj.sample_component_matching_probs(i,h);
%                 end
%                 
%                 % update the value
%                 update_result{h}=sum_atrs./sum_denominator;
%             end
%             
%             for h = 1:obj.number_of_components
%                 obj.mdl_ARGs{h}.edges_matrix=update_result{h};
%             end
%         end
        
        % update the covariance matrix for each component edge
        function updateComponentEdgeCov(obj)
            %for each component
            for h = 1:obj.number_of_components
                %for each edge 
                for o = 1:obj.mdl_ARGs{h}.num_nodes
                    for t = o+1:obj.mdl_ARGs{h}.num_nodes
                        if any(obj.mdl_ARGs{h}.edges{o,t}.getAtrs())
                            cov = 0;
                            denominator=0;
                            %for each sample
                            for i = 1:obj.number_of_sample
                                current_sample_cov = 0;
                                current_sample_denominator = 0;
                                %for each edge in sample
                                for c =  1:obj.sampleARGs{i}.num_nodes
                                    for d =  1:obj.sampleARGs{i}.num_nodes
                                        if any(obj.sampleARGs{i}.edges{c,d}.getAtrs())
                                            z_atrs=obj.sampleARGs{i}.edges{c,d}.getAtrs()-obj.mdl_ARGs{h}.edges{o,t}.getAtrs();
                                            current_sample_cov=current_sample_cov+...
                                                z_atrs'*z_atrs*obj.node_match_scores{i,h}(c,o)*obj.node_match_scores{i,h}(d,t);
                                            current_sample_denominator=current_sample_denominator+obj.node_match_scores{i,h}(c,o)*obj.node_match_scores{i,h}(d,t);
                                            %update experiment
%                                             obj.distance(obj.counter)=sqrt(z_atrs'*z_atrs);
%                                             obj.node_ab_match(obj.counter)=obj.node_match_scores{i,h}(c,o);
%                                             obj.node_cd_match(obj.counter)=obj.node_match_scores{i,h}(d,t);
%                                             obj.sample_match(obj.counter)=obj.sample_component_matching_probs(i,h);
%                                             obj.edge_match(obj.counter)=exp(-0.5*(z_atrs)*obj.mdl_ARGs{h}.edges{o,t}.getCovInv()*(z_atrs)')/((2*pi)^(1/2)*sqrt(det(obj.mdl_ARGs{h}.edges{o,t}.getCov())));
%                                             obj.base(obj.counter)=obj.base_counter;
%                                             obj.counter = obj.counter+1;
                                        end
                                    end
                                end
                                cov=cov+current_sample_cov*obj.sample_component_matching_probs(i,h);
                                denominator = denominator + current_sample_denominator*obj.sample_component_matching_probs(i,h);
                            end
                            % update the value
                            obj.mdl_ARGs{h}.edges{o,t}.updateCov(cov/denominator);
                            obj.mdl_ARGs{h}.edges{t,o}.updateCov(cov/denominator);
                        end
                    end
                end
            end     
            
        end
        
        % update the atrs for each component node
        % # this function can be easier, but I can think of a better way to
        % do this yet since atrs can be vector and cell operation is not
        % fast/easy
        function updateComponentNodeAtrs(obj)
            % for each component
            for h = 1:obj.number_of_components
                % for each node
                for n = 1:obj.mdl_ARGs{h}.num_nodes
                    if any(obj.mdl_ARGs{h}.nodes{n}.getAtrs())
                        atrs = 0;
                        denominator=0;
                        % we go over the sample
                        for i = 1:obj.number_of_sample
                            current_sample_atrs = 0;
                            current_sample_denominator = 0;
                            % and finds its matching node, calculate the
                            % average atrs
                            for m =  1:obj.sampleARGs{i}.num_nodes
                                if any(obj.sampleARGs{i}.nodes{m}.getAtrs())
                                    current_sample_atrs=current_sample_atrs+obj.sampleARGs{i}.nodes{m}.getAtrs()*obj.node_match_scores{i,h}(m,n);
                                    current_sample_denominator = current_sample_denominator + obj.node_match_scores{i,h}(m,n);
                                end
                            end
                            atrs = atrs + current_sample_atrs*obj.sample_component_matching_probs(i,h);
                            denominator = denominator + current_sample_denominator*obj.sample_component_matching_probs(i,h);
                        end
                        % udpate the value
                        newAtr = atrs/denominator;
                        newAtr = newAtr/sum(newAtr);
                        obj.mdl_ARGs{h}.nodes{n}.updateAtrs(newAtr);
                    end
                end
            end       
        end
        
        % update the frequency for each component node
        function updateComponentNodeFrequency(obj)
            % for each component
            for i = 1:obj.number_of_components
                frequency=zeros(obj.number_of_sample,obj.mdl_ARGs{i}.num_nodes);
                sample_node_sum = zeros(1,obj.number_of_sample);
                % for each sample
                parfor j = 1:obj.number_of_sample
                    % we calculate the component node frequency in this
                    % specific sample
                    current_freq =sum(obj.node_match_scores{j,i});
                    frequency(j,:)=current_freq*obj.sample_component_matching_probs(j,i);
                    sample_node_sum(j)=obj.sampleARGs{j}.num_nodes*obj.sample_component_matching_probs(j,i);
                end
                % update the frequency for all the nodes in the model in
                % the same time
                frequency=sum(frequency);
                sample_node_sum=sum(sample_node_sum);
                obj.mdl_ARGs{i}.updateNodeFrequency(frequency/sample_node_sum);
            end    
        end
        
        % update the weight for each component
        function updateComponentWeight(obj)
            % sum up the samples-component probabilities and divide by
            % sample number
            obj.weight = sum(obj.sample_component_matching_probs)/obj.number_of_sample;
        end
        
        % get the sample-component matching score and probability
        function getMatchingProbs(obj)
            % Use the graphMatching() data calculate the sample-componennt
            % matching score & probability
            handle=@(node_match_score,node_compatibility,edge_compatibility)sprMDL.component_score(node_match_score,node_compatibility,edge_compatibility,true);
            obj.component_scores=cellfun(handle,obj.node_match_scores,obj.node_compatibilities,obj.edge_compatibilities);
            s=sum(obj.component_scores,2);
            n=repmat(s,1,obj.number_of_components);
            obj. sample_component_matching_probs=obj.component_scores./n;
        end
        
        % get the graph matching score for each sample-componennt pair
        % **************************************passing the whole model?
        function graphMatching(obj)
            obj.node_match_scores = cell([obj.number_of_sample,obj.number_of_components]);
            obj.node_compatibilities = cell(size(obj.node_match_scores));
            obj.edge_compatibilities = cell(size(obj.node_match_scores));
            
            for i=1:obj.number_of_sample
                i_node_match_scores=cell(1,obj.number_of_components);
                i_node_compatibilities=cell(1,obj.number_of_components);
                i_edge_compatibilities=cell(1,obj.number_of_components);
                parfor j = 1:obj.number_of_components
                    [node_match_score,node_compatibility,edge_compatibility] = graph_matching(obj.sampleARGs{i},obj.mdl_ARGs{j},obj.BLOSUM);
                    i_node_match_scores{j}=node_match_score;
                    i_node_compatibilities{j}=node_compatibility;
                    i_edge_compatibilities{j}=edge_compatibility;
                end
                obj.node_match_scores(i,:)=i_node_match_scores;
                obj.node_compatibilities(i,:)=i_node_compatibilities;
                obj.edge_compatibilities(i,:)=i_edge_compatibilities;
            end
        end
    end
    
    methods(Static)
        % scoring for each sample-component pair
        function score = component_score(node_match_score,node_compatibility,edge_compatibility,train)
           % calculate the prob from the nodes part
           score = sum(sum(node_match_score.*node_compatibility));           
           % calculate the prob from the edges part
           edge_times_handle = @(mat)sum(sum(bsxfun(@times,node_match_score,mat)));
           first_time = cellfun(edge_times_handle,edge_compatibility);
           % add both part together
           score = score + sum(sum(first_time.*node_match_score));
           % nomalize by the number of nodes in the component if it is
           % during the training
           if train
               % experiment more with other option
                score = score/log(size(node_match_score,2));
           end
        end
       
    end
end

