classdef mdl_node < node
    % mdl_node is a subclass of node which will be used in the mdl_ARG
    
    methods
        % Constructor for the class
        function  obj = mdl_node(ID,ARG)
            % Throw error if not enough argument
            if nargin < 2
                error "NotEnoughArgument";
            end
            
            % Passing original value
            obj=obj@node(ID,ARG);
            
        end
        
        function updateFrequency(obj,freq)
            obj.ARG.nodes_frequency(obj.ID) = freq;
        end
        
        function [val] = getFrequency(obj)
            val = obj.ARG.nodes_frequency(obj.ID);
        end
    
        % Update Mean
        function updateAtrs(obj,atrs)
            obj.ARG.nodes_vector{obj.ID} = atrs;
        end
    end
    
end

