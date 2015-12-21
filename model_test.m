%% Script for matching the protein
tic();
draw_flag = 0;

%% Set Up Edge and Node Attribute

BLOSUM_Sigma = 2;   % node attribute

distance_cutoff = 13;   % edge attribute cut off

save = 1;

%% Set Up Protein

proteinOneFile = 'new_4D1E_CH.csv';
start_sequence_one = 60;
end_sequence_one = 150;
[proteinOneARG,~] = GenerateProteinARGs(start_sequence_one,end_sequence_one, proteinOneFile,distance_cutoff);

proteinTwoFile = 'new_4Q59_CH.csv';
start_sequence_two = 10;
end_sequence_two = 90;
[proteinTwoARG,~] = GenerateProteinARGs(start_sequence_two,end_sequence_two, proteinTwoFile,distance_cutoff);


%%

trainingSample = cell([1,2]);
trainingSample{1} = proteinOneARG;
trainingSample{2} = proteinTwoARG;
%% Build the model

MDL = sprMDL(trainingSample,2);

toc()

%% Save the result
if (save)
    datetime=datestr(now);
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    save (['Data/' datetime '.mat']) 
end