function [modelStruct1, modelStruct2] = fixfield(modelStruct1, modelStruct2 )

    f1 = fieldnames(modelStruct1);
    f2 = fieldnames(modelStruct2);
    
    if ~isequal(f1, f2)
        missingFrom1 = ~isfield(modelStruct1, f2);
        missingFrom2 = ~isfield(modelStruct2, f1);
        if any(missingFrom1)
            fields = f2(missingFrom1);
            for i = 1:numel(fields)
                modelStruct1.(fields{i}) = '';
            end
        end
        if any(missingFrom2)
           fields = f1(missingFrom2); 
           for i = 1:numel(fields)
               modelStruct2.(fields{i}) = '';
           end
        end
    end
end

