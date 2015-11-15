classdef mdl_edge < edge
    % mdl_edge is a subclass of edge which will be used in the mdl_ARG
    
    properties (GetAccess=public,SetAccess=protected)
        cov=NaN;
        cov_inv = NaN;
    end
    
    properties(Constant)
        e_inv = 0.2;
        conv_eye = 0.1;
    end
    
    methods
        % Constructor for the class
        function  obj = mdl_edge(weight,node1ID,node2ID,sortedNodes)
            % Throw error if not enough argument
            if nargin < 4
                error "NotEnoughArgument";
            end
            
            % Passing original value
            obj=obj@edge(weight,node1ID,node2ID,sortedNodes);
            
            % Initial covariance matrix as an identtiy matrix
            obj.cov = eye(length(weight));
            obj.cov_inv=inv(obj.cov);
        end
        
        % Update Mean
        function updateAtrs(obj,weight)
            obj.weight = weight;
        end
        
        % Update Covariance Matrix
        function updateCov(obj,cov)
            obj.cov = cov;
            obj.cov_inv = mdl_edge.inverse(obj.cov);
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
        end
    end
    
end

