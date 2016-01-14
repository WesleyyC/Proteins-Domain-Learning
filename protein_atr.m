function [ output_args ] = protein_atr( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    output_args = cell(size(input_args));
    input_args = round(input_args);
    
    for i = 1 : length(input_args)
        
        k = input_args(i);
        
        if k<1
            k=1;
        elseif k>20
            k = 20;
        end
            
        
        atr = zeros(1,20);
        atr(k) = 1;
        output_args{i} = atr;
    end
 
end

