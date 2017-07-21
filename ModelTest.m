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
proteinARGs{end+1} = GenerateProteinARG(140,230,'protein/SH2/1ad5.csv');
proteinARGs{end+1} = GenerateProteinARG(1,90,'protein/SH2/1aou.csv');
proteinARGs{end+1} = GenerateProteinARG(1,90,'protein/SH2/1fbz.csv');
proteinARGs{end+1} = GenerateProteinARG(140,230,'protein/SH2/1qcf.csv');
proteinARGs{end+1} = GenerateProteinARG(120,210,'protein/SH2/4d8k.csv');

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
