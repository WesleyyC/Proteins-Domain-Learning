%% Script for matching the protein
tic();
draw_flag = 0;

%% Set Up Edge and Node Attribute

BLOSUM_Sigma = 2;   % node attribute

distance_cutoff = 13;   % edge attribute cut off

save = 0;

%% Set Up Protein

proteinOneFile = 'new_4D1E_CH.csv';
start_sequence_one = 100;
end_sequence_one = 120;

proteinTwoFile = 'new_4Q59_CH.csv';
start_sequence_two = 100;
end_sequence_two = 120;

[proteinThreeARG,p1] = GenerateProteinARGs(start_sequence_one,end_sequence_one, proteinOneFile,distance_cutoff);
[proteinFourARG,p2] = GenerateProteinARGs(start_sequence_two,end_sequence_two, proteinTwoFile,distance_cutoff);


trainingSample = cell([1,2]);
trainingSample{1} = proteinOneARG;
trainingSample{2} = proteinTwoARG;


%% Build the model

MDL = sprMDL(trainingSample,2);

toc()