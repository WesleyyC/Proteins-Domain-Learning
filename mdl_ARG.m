classdef mdl_ARG < handle
    %   mdl_ARG represetns a component in our model
    properties (GetAccess=public,SetAccess=public)
        num_nodes = NaN;
        
        nodes = {};
        nodes_vector = NaN;
        nodes_frequency = NaN;
        
        edges = {};
        edges_matrix = NaN;
        edges_cov = NaN;
        edges_cov_inv = NaN;
    end
    
    methods
        % setting up constructor which will take an sample ARG and build a
        % new component for the model.
        function self = mdl_ARG(ARG)
            % Throw error if not enough argument
            if nargin < 1
                error "NotEnoughArgument";
            end
            
            % Get the number of nodes
            self.num_nodes=ARG.num_nodes+1;
            
            % Build the nodes cell
            self.nodes = cell(1,self.num_nodes);
            mdl_node_handle=@(node)mdl_node(node.ID,self);
            self.nodes(1:self.num_nodes-1) = cellfun(mdl_node_handle,ARG.nodes,'UniformOutput',false);
            self.nodes{self.num_nodes} = mdl_node(self.num_nodes, self);
            % Build the nodes_vector
            self.nodes_vector = cell(1,self.num_nodes);
            self.nodes_vector(1:self.num_nodes-1)=ARG.nodes_vector;
            self.nodes_vector{self.num_nodes}=zeros(1,20);
            % Buil the nodes_frequency
            freq = 1/self.num_nodes;
            self.nodes_frequency = ones(1,self.num_nodes)*freq;
            
            % Build the edges cell
            self.edges = cell(self.num_nodes,self.num_nodes);
            mdl_edge_handle=@(edge)mdl_edge(self,self.nodes{edge.node1.ID},self.nodes{edge.node2.ID});
            self.edges(1:self.num_nodes-1,1:self.num_nodes-1) = cellfun(mdl_edge_handle,ARG.edges,'UniformOutput',false);
            for i=1:self.num_nodes
                self.edges{self.num_nodes,i}=mdl_edge(self,self.nodes{self.num_nodes},self.nodes{i});
                self.edges{i,self.num_nodes}=mdl_edge(self,self.nodes{i},self.nodes{self.num_nodes});
            end
            % Build the edges_matrix
            self.edges_matrix = ARG.edges_matrix;
            for i=1:self.num_nodes
                self.edges_matrix(self.num_nodes,i)=0;
                self.edges_matrix(i,self.num_nodes)=0;
            end
            % Build the edges_cov
            mdl_edge_cov_handle = @(edge)eye(length(edge.getAtrs()));
            self.edges_cov = cell2mat(cellfun(mdl_edge_cov_handle,self.edges,'UniformOutput',false));
            % Build the edges_cov_inv
            mdl_edge_cov_inv_handle = @(edge)inv(edge.getCov());
            self.edges_cov_inv = cell2mat(cellfun(mdl_edge_cov_inv_handle,self.edges,'UniformOutput',false));
                
        end
        
        % update the nodes frequency in the model
        function updateNodeFrequency(obj,frequencies)
            obj.nodes_frequency = frequencies;
        end
        
        % delete nodes in the model according to the given indexes
        function modifyStructure(obj,deletingNodes)
            obj.nodes(deletingNodes)=[];
            obj.edges(deletingNodes,:)=[];
            obj.edges(:,deletingNodes)=[];
            obj.num_nodes=obj.num_nodes-length(find(deletingNodes));
        end
        
        % return a vector of nodes frequency
        function frequencies = getNodeFrequency(obj)
           frequencies = obj.nodes_frequency;
        end
        
        % show the model ARG in matrix
        function model_struct = showARG(obj)
            % get the nodes frequencies
            nodes_frequence = obj.getNodeFrequency();
            
            % get the nodes attributes
            getNodeAttributes = @(node)node.atrs;
            nodes_attributes = cellfun(getNodeAttributes,obj.nodes);
            if~iscell(nodes_attributes)
                nodes_attributes=num2cell(nodes_attributes);
            end
            % build the matrix
            getEdgeAttributes = @(edge)edge.atrs;
            M = cellfun(getEdgeAttributes,obj.edges);
            
            % draw it out
            bg = biograph(triu(M(1:end-1,1:end-1)),[],'ShowArrows','off','ShowWeights','on');

            model_struct = struct(  'Nf',nodes_frequence,... % the node frequency
                                    'M', M, ... % the matrix
                                    'bg', bg...% the graph
                                    );
            model_struct.Na = nodes_attributes;
        end
            
    end
    
end

