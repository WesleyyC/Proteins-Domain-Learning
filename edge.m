classdef edge < handle
    %   edge is the connection between node and 
    %   it will have some assigned weight and two end points (nodes)
    
    properties (GetAccess=public,SetAccess=protected)
        % The attributes
        weight = NaN;   % the weight for the link, and be NaN if there is no link between node1 and node2
        node1 = NaN;
        node2 = NaN;
        node1ID = NaN;
        node2ID = NaN;
    end
    
    methods
        % Constructor for the class
        function  self = edge(weight,node1ID,node2ID,sortedNodes)
            % Throw error if not enough argument
            if nargin < 4
                error "NotEnoughArgument";
            end
            
            % Otherwise, we process the argument
            self.weight = weight;
            
            self.node1 = sortedNodes{node1ID};
            self.node2 = sortedNodes{node2ID};
            
            self.node1ID = node1ID;
            self.node2ID = node2ID;
        end
        
        % Get the similarity between two edges
        function [c] = compatibility(obj,obj2)
            c = edge_compatibility(obj,obj2);
        end
        
        function [tf] = trueEdge(obj)
            tf=obj.weight~=0;
        end
        
        % Get number of attributes
        function [no] = numberOfAtrs(obj)
            no=length(obj.weight);
        end
    end
    
end

