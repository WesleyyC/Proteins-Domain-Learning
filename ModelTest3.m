
%{
Note for domains:
Test1: 95% similarity; 
Test2: 60%  similarity;
Test3: sfam similarity;
Test4: same Topol
Test5: same Arch
Test6: same Class
%}
%% sprMDL auto testing
clear

%% Test1: 95% similarity;
proteinARGs_1 = cell(0);
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/1fjsL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/1z6eL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2bokL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2jkhL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2p3tA00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2p3uA00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2pr3B00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2vwlL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2vwnL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2vwoL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2w3iB00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2xbvL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2xbwL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2xbxL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2xbyL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2xc0L00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2xc5L00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/2y5fL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/3cenL00.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein/test1/3m36L00.csv');

%% Test2: 60%  similarity;
proteinARGs_2 = cell(0);
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/1lpgA00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/1nfuB00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/1nfyB00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/1xkaL02.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/1xkbA02.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/1xkbB02.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/2j2uB00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/2j4iB00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/2uwlB00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/2uwoB00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/2vh0B00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/2wygB00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/2y7zB00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/2y81B00.csv');
proteinARGs_2{end+1} = GenerateProteinARG('protein/test2/3ensA02.csv');

%% Test3: sfam similarity;
proteinARGs_3 = cell(0);
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/1adxA00.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/1cqeB01.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/1dqbA02.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/1dx5L02.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/1xkaL01.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/2bo2A02.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/2bouA02.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/2boxA02.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/3ensA01.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/3hptA01.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/3nt1A01.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/3ntbB01.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/3qcwA06.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/3r05A06.csv');
proteinARGs_3{end+1} = GenerateProteinARG('protein/test3/4ishL00.csv');
%% Test4: same Topol
proteinARGs_4 = cell(0);
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1kkeA01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1kkeB01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1kkeC01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1lk9A01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1lk9B01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1qiuC01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1qiuE01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1v1hA01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1v1hF01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/1v1iB01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/2horA01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/2hoxA01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/2hoxB01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/2hoxC01.csv');
proteinARGs_4{end+1} = GenerateProteinARG('protein/test4/2hoxD01.csv');

%% Test5: same Arch
proteinARGs_5 = cell(0);
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1bx8A00.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1c9tH00.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1cklE01.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1cklF01.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1eakA03.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1h8pA02.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1l6jA03.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1pdcA00.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1skzA01.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/1skzA02.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/2o39C01.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/2o39D01.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/2uwpB00.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/3bg4D00.csv');
proteinARGs_5{end+1} = GenerateProteinARG('protein/test5/4ayeB01.csv');

%% Test6: same Class
proteinARGs_6 = cell(0);
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/1b8wA00.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/1nh2C00.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/1nvpC00.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/1rm1B02.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/1rm1C02.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/1ytfC00.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/2de6C03.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/2sh1A00.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/3gtqI01.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/3rzdI01.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/4a3bI01.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/4gv5B00.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/4nbbB03.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/4nbcC03.csv');
proteinARGs_6{end+1} = GenerateProteinARG('protein/test6/4nbfB03.csv');
%% Training
mdl = sprMDL(proteinARGs_1, 2);
save model.mat mdl;

%% Scores
original_score_1 = scores(proteinARGs_1, mdl);
original_score_2 = scores(proteinARGs_2, mdl);
original_score_3 = scores(proteinARGs_3, mdl);
original_score_4 = scores(proteinARGs_4, mdl);
original_score_5 = scores(proteinARGs_5, mdl);
original_score_6 = scores(proteinARGs_6, mdl);

%% Score function

function original_score = scores(proteinARGs, mdl)
original_score = zeros([1,length(proteinARGs)]);
for i = 1:length(proteinARGs)
    [result, score] = mdl.checkPattern(proteinARGs{i});
    original_score(i) = score;
end
end



%% Detect Rate Calculation
%{
original_result = zeros([1,length(proteinARGs)]);
original_score = zeros([1,length(proteinARGs)]);
for i = 1:length(proteinARGs)
    [result, score] = mdl.checkPattern(proteinARGs{i});
    original_result(i) = result;
    original_score(i) = score;
end
original_detect_rate = sum(original_result)/length(original_score);
%}

%{
%% compare reuslt
fig = figure;
hax = axes;
nbin = 50;

histogram(original_score_1,nbin)
hold on
%line([mdl.thredshold_score mdl.thredshold_score],get(hax,'YLim'),'Color','g','LineWidth', 2)
%hold on
%legend('original')

histogram(original_score_2,nbin)
hold on
histogram(original_score_3,nbin)
hold on
histogram(original_score_4,nbin)
hold on
histogram(original_score_5,nbin)
hold on
%histogram(original_score_6,nbin)
%hold on
%line([mdl.thredshold_score mdl.thredshold_score],get(hax,'YLim'),'Color','g','LineWidth', 2)
hold on
%legend('100%','60%','sfam','same_topol', 'same_arch', 'same_class')
legend('100%','60%','sfam','same_topol', 'same_arch')
%}


%{ 
ORIGINAL
%% sprMDL auto testing
clear

%% build protiens

proteinARGs = cell(0);

% entire CH1
% proteinARGs = [proteinARGs, GenerateProteinARG(32,137,'protein/CH1/1sjj.csv')];
% proteinARGs = [proteinARGs, GenerateProteinARG(45,150,'protein/CH1/1tjt.csv')];
% proteinARGs = [proteinARGs, GenerateProteinARG(45,150,'protein/CH1/1wku.csv')];
% proteinARGs = [proteinARGs, GenerateProteinARG(31,136,'protein/CH1/2eyi.csv')];
% proteinARGs = [proteinARGs, GenerateProteinARG(38,143,'protein/CH1/4d1e.csv')];
% proteinARGs{end+1} = GenerateProteinARG(26,140,'protein/CH1/1sjj.csv');
% proteinARGs{end+1} = GenerateProteinARG(42,150,'protein/CH1/1tjt.csv');
% proteinARGs{end+1} = GenerateProteinARG(42,150,'protein/CH1/1wku.csv');
% proteinARGs{end+1} = GenerateProteinARG(26,140,'protein/CH1/2eyi.csv');
% proteinARGs{end+1} = GenerateProteinARG(34,150,'protein/CH1/4d1e.csv');

% partial SH3
proteinARGs{end+1} = GenerateProteinARG('protein/test1/1ojyA03.csv');
proteinARGs{end+1} = GenerateProteinARG('protein/test1/1ojyB03.csv');
proteinARGs{end+1} = GenerateProteinARG('protein/test1/1ojyC03.csv');
proteinARGs{end+1} = GenerateProteinARG('protein/test1/1ojyD03.csv');
%proteinARGs{end+1} = GenerateProteinARG('protein/newDomain/1vveA01.csv');

% mix sequence
% proteinARGs{end+1} = GenerateProteinARG(1,40,'protein/Mix/2nch.csv');
% proteinARGs{end+1} = GenerateProteinARG(11,51 ,'protein/Mix/3jcu.csv');
% proteinARGs{end+1} = GenerateProteinARG(123,149 ,'protein/Mix/3wtq.csv');
% proteinARGs{end+1} = GenerateProteinARG(61,101,'protein/Mix/4pc7.csv');
% proteinARGs{end+1} = GenerateProteinARG(24,65 ,'protein/Mix/5fzi.csv');
% proteinARGs{end+1} = GenerateProteinARG(3,53,'protein/Mix/5h9l.csv');
% proteinARGs{end+1} = GenerateProteinARG(3,53,'protein/Mix/5tlq.csv');

%% load the helper function

mdl = sprMDL(proteinARGs, 2);

%% check original smaple

original_result = zeros([1,length(proteinARGs)]);
original_score = zeros([1,length(proteinARGs)]);

for i = 1:length(proteinARGs)
    [result, score] = mdl.checkPattern(proteinARGs{i});
    original_result(i) = result;
    original_score(i) = score;
end

original_detect_rate = sum(original_result)/length(original_result)

%%file printout
fileID = fopen('test1.txt','w');
fprintf(fileID, 'original_score = \n');
fprintf(fileID,'%f\n',original_score);
fprintf(fileID, 'original_sresult = \n');
fprintf(fileID,'%f\n',original_result);
fprintf(fileID, 'original_detect_rate = \n');
fprintf(fileID,'%f\n',original_detect_rate);

%% reverse order

reverse_result = zeros([1,length(proteinARGs)]);
reverse_score = zeros([1,length(proteinARGs)]);

for i = 1:length(proteinARGs)
    tmpARG = proteinARGs{i}.copy;
    r_idx = tmpARG.num_nodes:-1:1;
    tmpARG.nodes = tmpARG.nodes(r_idx);
    tmpARG.nodes_vector = tmpARG.nodes_vector(r_idx, :);
    tmpARG.edges = tmpARG.edges(r_idx, r_idx);
    tmpARG.edges_matrix = tmpARG.edges_matrix(r_idx, r_idx, :);
    [result,score] = mdl.checkPattern(tmpARG);
    reverse_result(i) = result;
    reverse_score(i) = score;
end

reverse_detect_rate = sum(reverse_result)/length(reverse_result)

%% partition

partition_num = 3;

partition_result = zeros([1,length(proteinARGs)]);
partition_score = zeros([1,length(proteinARGs)]);

for i = 1:length(proteinARGs)
    tmpARG = proteinARGs{i}.copy;
    p_idx = partition_idx(partition_num,tmpARG.num_nodes);
    tmpARG.nodes = tmpARG.nodes(p_idx);
    tmpARG.nodes_vector = tmpARG.nodes_vector(p_idx, :);
    tmpARG.edges = tmpARG.edges(p_idx, p_idx);
    tmpARG.edges_matrix = tmpARG.edges_matrix(p_idx, p_idx, :);
    [result, score] = mdl.checkPattern(tmpARG);
    partition_result(i) = result;
    partition_score(i) = score;
end

partition_detect_rate = sum(partition_result)/length(partition_result)

%% random test

random_sample_number = 10;

random_result = zeros([1,length(random_sample_number)]);
random_score = zeros([1,length(random_sample_number)]);

count = 1;

for i = 1:length(proteinARGs)
    for j = 1:random_sample_number/length(proteinARGs)
        tmpARG = proteinARGs{i}.copy;
        n_idx = randperm(tmpARG.num_nodes);
        e_idx = randperm(tmpARG.num_nodes);
        tmpARG.nodes = tmpARG.nodes(n_idx);
        tmpARG.nodes_vector = tmpARG.nodes_vector(n_idx, :);
        tmpARG.edges = tmpARG.edges(e_idx, e_idx);
        tmpARG.edges_matrix = tmpARG.edges_matrix(e_idx, e_idx, :);
        [result, score] = mdl.checkPattern(tmpARG);
        random_result(count) = result;
        random_score(count) = score;
        count = count + 1;
    end
end

random_detect_rate = sum(random_result)/length(random_result)

%% compare reuslt
fig = figure;
hax = axes;
nbin = 50;

histogram(original_score,nbin)
hold on
line([mdl.thredshold_score mdl.thredshold_score],get(hax,'YLim'),'Color','g','LineWidth', 2)
hold on
legend('original')

histogram(reverse_score,nbin)
hold on
histogram(partition_score,nbin)
hold on
histogram(random_score,nbin)
hold on
line([mdl.thredshold_score mdl.thredshold_score],get(hax,'YLim'),'Color','g','LineWidth', 2)
hold on
legend('original','reverse','partition','random')
%% visual   
% visual(32,57,32,40,'protein/CH1/1sjj.csv',45,70,45,53,'protein/CH1/1wku.csv',mdl.node_match_scores{1,1})
% visual(32,57,NaN,NaN,'protein/CH1/1sjj.csv',45,70,NaN,NaN,'protein/CH1/1wku.csv',mdl.node_match_scores{1,1})
%}
