%% Script for matching the protein
tic();
draw_flag = 0;

%% Set Up Edge and Node Attribute

BLOSUM_Sigma = 2;   % node attribute

distance_cutoff = 13;   % edge attribute cut off

save = 0;

%% Set Up Protein

proteinOneFile = 'new_4D1E_CH.csv';
start_sequence_one = 90;
end_sequence_one = 140;
[proteinOneARG,~] = GenerateProteinARGs(start_sequence_one,end_sequence_one, proteinOneFile,distance_cutoff);

proteinTwoFile = 'new_4Q59_CH.csv';
start_sequence_two = 1;
end_sequence_two = 30;
[proteinTwoARG,~] = GenerateProteinARGs(start_sequence_two,end_sequence_two, proteinTwoFile,distance_cutoff);

% proteinOneFile = 'new_4D1E_CH.csv';
% start_sequence_one = 110;
% end_sequence_one = 135;
% [proteinThreeARG,~] = GenerateProteinARGs(start_sequence_one,end_sequence_one, proteinOneFile,distance_cutoff);
% 
% proteinTwoFile = 'new_4Q59_CH.csv';
% start_sequence_two = 120;
% end_sequence_two = 145;
% [proteinFourARG,~] = GenerateProteinARGs(start_sequence_two,end_sequence_two, proteinTwoFile,distance_cutoff);
% 
% proteinOneFile = 'new_4D1E_CH.csv';
% start_sequence_one = 140;
% end_sequence_one = 165;
% [proteinFiveARG,~] = GenerateProteinARGs(start_sequence_one,end_sequence_one, proteinOneFile,distance_cutoff);
% 
% proteinTwoFile = 'new_4Q59_CH.csv';
% start_sequence_two = 150;
% end_sequence_two = 175;
% [proteinSixARG,~] = GenerateProteinARGs(start_sequence_two,end_sequence_two, proteinTwoFile,distance_cutoff);

%%

trainingSample = cell([1,2]);
trainingSample{1} = proteinOneARG;
trainingSample{2} = proteinTwoARG;
% trainingSample{3} = proteinThreeARG;
% trainingSample{4} = proteinFourARG;
% trainingSample{5} = proteinFiveARG;
% trainingSample{6} = proteinSixARG;
%% Build the model

MDL = sprMDL(trainingSample,2);

toc()