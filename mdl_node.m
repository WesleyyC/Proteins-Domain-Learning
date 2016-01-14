classdef mdl_node < node
    % mdl_node is a subclass of node which will be used in the mdl_ARG
    
    properties (GetAccess=public,SetAccess=protected)
        frequency = NaN;
    end
    
    methods
        % Constructor for the class
        function  obj = mdl_node(ID,atrs,frequency)
            % Throw error if not enough argument
            if nargin < 2
                error "NotEnoughArgument";
            end
            
            % Passing original value
            obj=obj@node(ID,atrs);
            
            % Initial frequency
            obj.frequency = frequency;
        end
        
        function updateFrequency(obj,freq)
            obj.frequency = freq;
        end
    
        % Update Mean
        function updateAtrs(obj,atrs)
            obj.atrs = atrs;
        end
    end
    
end

