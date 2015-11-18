%% Script for matching the protein
tic();
draw_flag = 0;

%% Set Up Edge and Node Attribute

BLOSUM_Sigma = 2;   % node attribute

distance_cutoff = 13;   % edge attribute cut off

save = 1;

%% Set Up Protein

proteinOneFile = 'new_1G83_65_147.csv';
start_sequence_one = 60;
end_sequence_one = 150;
[proteinOneARG,~] = GenerateProteinARGs(start_sequence_one,end_sequence_one, proteinOneFile,distance_cutoff);

proteinTwoFile = 'new_1WQU_18_88.csv';
start_sequence_two = 10;
end_sequence_two = 90;
[proteinTwoARG,~] = GenerateProteinARGs(start_sequence_two,end_sequence_two, proteinTwoFile,distance_cutoff);

proteinThreeFile = 'new_1XA6_34_107.csv';
start_sequence_three = 30;
end_sequence_three = 110;
[proteinThreeARG,~] = GenerateProteinARGs(start_sequence_three,end_sequence_three, proteinThreeFile,distance_cutoff);

proteinFourFile = 'new_2ABL_72_147.csv';
start_sequence_four = 70;
end_sequence_four = 150;
[proteinFourARG,~] = GenerateProteinARGs(start_sequence_four,end_sequence_four, proteinFourFile,distance_cutoff);

proteinFiveFile = 'new_2CR4_20_100.csv';
start_sequence_five = 20;
end_sequence_five = 100;
[proteinFiveARG,~] = GenerateProteinARGs(start_sequence_five,end_sequence_five, proteinFiveFile,distance_cutoff);

proteinSixFile = 'new_2EKX_13_94.csv';
start_sequence_six = 10;
end_sequence_six = 100;
[proteinSixARG,~] = GenerateProteinARGs(start_sequence_six,end_sequence_six, proteinSixFile,distance_cutoff);

proteinSevenFile = 'new_2LQW_14_88.csv';
start_sequence_seven = 10;
end_sequence_seven = 90;
[proteinSevenARG,~] = GenerateProteinARGs(start_sequence_seven,end_sequence_seven, proteinSevenFile,distance_cutoff);

proteinEightFile = 'new_3OV1_6_81.csv';
start_sequence_eight = 1;
end_sequence_eight = 90;
[proteinEightARG,~] = GenerateProteinARGs(start_sequence_eight,end_sequence_eight, proteinEightFile,distance_cutoff);

proteinNineFile = 'new_3VS0_65_147.csv';
start_sequence_nine = 60;
end_sequence_nine = 147;
[proteinNineARG,~] = GenerateProteinARGs(start_sequence_nine,end_sequence_nine, proteinNineFile,distance_cutoff);

proteinTenFile = 'new_4EIH_7_82.csv';
start_sequence_ten = 10;
end_sequence_ten = 90;
[proteinTenARG,~] = GenerateProteinARGs(start_sequence_ten,end_sequence_ten, proteinTenFile,distance_cutoff);


%%

trainingSample = cell([1,10]);
trainingSample{1} = proteinOneARG;
trainingSample{2} = proteinTwoARG;
trainingSample{3} = proteinThreeARG;
trainingSample{4} = proteinFourARG;
trainingSample{5} = proteinFiveARG;
trainingSample{6} = proteinSixARG;
trainingSample{7} = proteinSevenARG;
trainingSample{8} = proteinEightARG;
trainingSample{9} = proteinNineARG;
trainingSample{10} = proteinTenARG;

%% Build the model

MDL = sprMDL(trainingSample,3);

toc()

%% Save the result
if (save)
    datetime=datestr(now);
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    save (['Data/' datetime '.mat']) 
end