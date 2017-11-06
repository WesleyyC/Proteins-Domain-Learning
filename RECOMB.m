% %% same 95%
% 
% proteinARGs = cell(0);
% proteinARGs{end+1} = GenerateProteinARG('test-protein/test1/2g7iA01.csv');
% proteinARGs{end+1} = GenerateProteinARG('test-protein/test1/4ontF01.csv');
% proteinARGs{end+1} = GenerateProteinARG('test-protein/test1/3r62B01.csv');
% proteinARGs{end+1} = GenerateProteinARG('test-protein/test1/3kzjA01.csv');
% proteinARGs{end+1} = GenerateProteinARG('test-protein/test1/3oxuD01.csv');
% 
% mdl_1 = sprMDL(proteinARGs, 1);

%% same sfam

proteinARGs = cell(0);
proteinARGs{end+1} = GenerateProteinARG('test-protein/test3/1cklC01.csv');
proteinARGs{end+1} = GenerateProteinARG('test-protein/test3/1ojwB01.csv');
proteinARGs{end+1} = GenerateProteinARG('test-protein/test3/2w81A02.csv');
proteinARGs{end+1} = GenerateProteinARG('test-protein/test3/2wiiC01.csv');
proteinARGs{end+1} = GenerateProteinARG('test-protein/test3/3o8eD02.csv');

mdl_2 = sprMDL(proteinARGs, 2);

%%

mdl_test = mdl_2;

%% 95%
proteinARGs_test1 = cell(0);
proteinARGs_test1{end+1} = GenerateProteinARG('test-protein/test1/3oxuF01.csv');
proteinARGs_test1{end+1} = GenerateProteinARG('test-protein/test1/2bzmA01.csv');
proteinARGs_test1{end+1} = GenerateProteinARG('test-protein/test1/3r62A01.csv');
proteinARGs_test1{end+1} = GenerateProteinARG('test-protein/test1/3rj3E01.csv');
proteinARGs_test1{end+1} = GenerateProteinARG('test-protein/test1/4ontE01.csv');

%%

test1_result = zeros([1,length(proteinARGs_test1)]);
test1_score = zeros([1,length(proteinARGs_test1)]);

for i = 1:length(proteinARGs_test1)
    [result,score] = mdl_test.checkPattern(proteinARGs_test1{i});
    test1_result(i) = result;
    test1_score(i) = score;
end

%% 60%

proteinARGs_test2 = cell(0);
proteinARGs_test2{end+1} = GenerateProteinARG('test-protein/test2/2xqwC01.csv');
proteinARGs_test2{end+1} = GenerateProteinARG('test-protein/test2/3kxvA01.csv');
proteinARGs_test2{end+1} = GenerateProteinARG('test-protein/test2/3zd1A01.csv');
proteinARGs_test2{end+1} = GenerateProteinARG('test-protein/test2/3zd1B01.csv');
proteinARGs_test2{end+1} = GenerateProteinARG('test-protein/test2/4j38B01.csv');

%%

test2_result = zeros([1,length(proteinARGs_test2)]);
test2_score = zeros([1,length(proteinARGs_test2)]);

for i = 1:length(proteinARGs_test2)
    [result,score] = mdl_test.checkPattern(proteinARGs_test2{i});
    test2_result(i) = result;
    test2_score(i) = score;
end

%% sfam

proteinARGs_test3 = cell(0);
proteinARGs_test3{end+1} = GenerateProteinARG('test-protein/test3/1cklC01.csv');
proteinARGs_test3{end+1} = GenerateProteinARG('test-protein/test3/1ojwB01.csv');
proteinARGs_test3{end+1} = GenerateProteinARG('test-protein/test3/2w81A02.csv');
proteinARGs_test3{end+1} = GenerateProteinARG('test-protein/test3/2wiiC01.csv');
proteinARGs_test3{end+1} = GenerateProteinARG('test-protein/test3/3o8eD02.csv');

%%

test3_result = zeros([1,length(proteinARGs_test3)]);
test3_score = zeros([1,length(proteinARGs_test3)]);

for i = 1:length(proteinARGs_test3)
    [result,score] = mdl_test.checkPattern(proteinARGs_test3{i});
    test3_result(i) = result;
    test3_score(i) = score;
end

%% topol

proteinARGs_test4 = cell(0);
proteinARGs_test4{end+1} = GenerateProteinARG('test-protein/test4/3k3tA02.csv');
proteinARGs_test4{end+1} = GenerateProteinARG('test-protein/test4/1rlhA01.csv');
proteinARGs_test4{end+1} = GenerateProteinARG('test-protein/test4/3ci0J02.csv');
proteinARGs_test4{end+1} = GenerateProteinARG('test-protein/test4/3vwoA02.csv');
proteinARGs_test4{end+1} = GenerateProteinARG('test-protein/test4/2zycA02.csv');

%%

test4_result = zeros([1,length(proteinARGs_test4)]);
test4_score = zeros([1,length(proteinARGs_test4)]);

for i = 1:length(proteinARGs_test4)
    [result,score] = mdl_test.checkPattern(proteinARGs_test4{i});
    test4_result(i) = result;
    test4_score(i) = score;
end

%% arch

proteinARGs_test5 = cell(0);
proteinARGs_test5{end+1} = GenerateProteinARG('test-protein/test5/2y5gL00.csv');
proteinARGs_test5{end+1} = GenerateProteinARG('test-protein/test5/2a06V00.csv');
proteinARGs_test5{end+1} = GenerateProteinARG('test-protein/test5/2r33B00.csv');
proteinARGs_test5{end+1} = GenerateProteinARG('test-protein/test5/2y5fL00.csv');
proteinARGs_test5{end+1} = GenerateProteinARG('test-protein/test5/3gbmA01.csv');

%%

test5_result = zeros([1,length(proteinARGs_test5)]);
test5_score = zeros([1,length(proteinARGs_test5)]);

for i = 1:length(proteinARGs_test5)
    [result,score] = mdl_test.checkPattern(proteinARGs_test5{i});
    test5_result(i) = result;
    test5_score(i) = score;
end

%% class

proteinARGs_test6 = cell(0);
proteinARGs_test6{end+1} = GenerateProteinARG('test-protein/test6/1ei5A02.csv');
proteinARGs_test6{end+1} = GenerateProteinARG('test-protein/test6/1yvyB02.csv');
proteinARGs_test6{end+1} = GenerateProteinARG('test-protein/test6/2sh1A00.csv');
proteinARGs_test6{end+1} = GenerateProteinARG('test-protein/test6/3ku3A01.csv');
proteinARGs_test6{end+1} = GenerateProteinARG('test-protein/test6/3rioA01.csv');

%%

test6_result = zeros([1,length(proteinARGs_test6)]);
test6_score = zeros([1,length(proteinARGs_test6)]);

for i = 1:length(proteinARGs_test6)
    [result,score] = mdl_test.checkPattern(proteinARGs_test6{i});
    test6_result(i) = result;
    test6_score(i) = score;
end