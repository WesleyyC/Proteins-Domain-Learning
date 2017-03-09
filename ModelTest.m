%% sprMDL auto testing
clear

%% build protiens

proteinARGs = cell(0);
% proteinARGs = [proteinARGs, GenerateProteinARG(32,137,'protein/CH1/1sjj.csv')];
% proteinARGs = [proteinARGs, GenerateProteinARG(45,150,'protein/CH1/1tjt.csv')];
% proteinARGs = [proteinARGs, GenerateProteinARG(45,150,'protein/CH1/1wku.csv')];
% proteinARGs = [proteinARGs, GenerateProteinARG(31,136,'protein/CH1/2eyi.csv')];
% proteinARGs = [proteinARGs, GenerateProteinARG(38,143,'protein/CH1/4d1e.csv')];
proteinARGs{end+1} = GenerateProteinARG(26,61,'protein/CH1/1sjj.csv');
proteinARGs{end+1} = GenerateProteinARG(42,67,'protein/CH1/1tjt.csv');
proteinARGs{end+1} = GenerateProteinARG(42,67,'protein/CH1/1wku.csv');
proteinARGs{end+1} = GenerateProteinARG(26,51,'protein/CH1/2eyi.csv');
proteinARGs{end+1} = GenerateProteinARG(34,59,'protein/CH1/4d1e.csv');

%% load the helper function

mdl = sprMDL(proteinARGs, 2);

%% check original smaple

original_result = zeros([1,length(proteinARGs)]);

for i = 1:length(proteinARGs)
    original_result(i)=mdl.checkPattern(proteinARGs{i});
end

original_detect_rate = sum(original_result)/length(original_result)

%% reverse order

reverse_result = zeros([1,length(proteinARGs)]);

for i = 1:length(proteinARGs)
    tmpARG = proteinARGs{i}.copy;
    r_idx = tmpARG.num_nodes:-1:1;
    tmpARG.nodes = tmpARG.nodes(r_idx);
    tmpARG.nodes_vector = tmpARG.nodes_vector(r_idx, :);
    tmpARG.edges = tmpARG.edges(r_idx, r_idx);
    tmpARG.edges_matrix = tmpARG.edges_matrix(r_idx, r_idx, :);
    reverse_result(i)=mdl.checkPattern(tmpARG);
end

reverse_detect_rate = sum(reverse_result)/length(reverse_result)

%% partition

partition_num = 4;

partition_result = zeros([1,length(proteinARGs)]);

for i = 1:length(proteinARGs)
    tmpARG = proteinARGs{i}.copy;
    p_idx = partition_idx(partition_num,tmpARG.num_nodes);
    tmpARG.nodes = tmpARG.nodes(p_idx);
    tmpARG.nodes_vector = tmpARG.nodes_vector(p_idx, :);
    tmpARG.edges = tmpARG.edges(p_idx, p_idx);
    tmpARG.edges_matrix = tmpARG.edges_matrix(p_idx, p_idx, :);
    partition_result(i)=mdl.checkPattern(tmpARG);
end

partition_detect_rate = sum(partition_result)/length(partition_result)

%% random test

random_sample_number = 10;

random_result = zeros([1,length(random_sample_number)]);

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
        random_result(count) = mdl.checkPattern(tmpARG);
        count = count + 1;
    end
end

random_detect_rate = sum(random_result)/length(random_result)
