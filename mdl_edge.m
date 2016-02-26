classdef mdl_edge < edge
    % mdl_edge is a subclass of edge which will be used in the mdl_ARG
    
    properties(Constant)
        e_inv = 0.2;
        conv_eye = 0.1;
    end
    
    methods
        % Constructor for the class
        function  obj = mdl_edge(ARG,node1,node2)
            % Throw error if not enough argument
            if nargin < 3
                error "NotEnoughArgument";
            end
            
            % Passing original value
            obj=obj@edge(ARG,node1,node2);

        end
        
        function [val] = getAtrs(obj)
            val=obj.ARG.edges_matrix(obj.node1.ID,obj.node2.ID);
        end
        
        % Update Mean
        function updateAtrs(obj,weight)
            obj.ARG.edges_matrix(obj.node1.ID,obj.node2.ID) = weight;
        end
        
        function [val] = getCov(obj)
            val=obj.ARG.edges_cov(obj.node1.ID,obj.node2.ID);
        end
        
        function [val] = getCovInv(obj)
            val=obj.ARG.edges_cov_inv(obj.node1.ID,obj.node2.ID);
        end
        
        % Update Covariance Matrix
        function updateCov(obj,cov)
            obj.ARG.edges_cov(obj.node1.ID,obj.node2.ID) = cov;
            obj.ARG.edges_cov_inv(obj.node1.ID,obj.node2.ID) = obj.inverse(cov);
        end
    end
    
    methods(Static)
        % In case there is a singularity problem
        function cov_inv = inverse(mat)
            if rcond(mat) < mdl_edge.e_inv
                mat=mat+eye(size(mat))*mean2(mat)*mdl_edge.conv_eye;
                cov_inv = mdl_edge.inverse(mat);
            else
                cov_inv = inv(mat);
            end 
%             cov_inv = 1/mat;
        end
    end
    
end

