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

% partial SH2
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

%%
proteinARGs_1 = cell(0);
proteinARGs_1{end+1} = GenerateProteinARG('protein-new/test1/1ojyA03.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein-new/test1/1ojyB03.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein-new/test1/1ojyC03.csv');
proteinARGs_1{end+1} = GenerateProteinARG('protein-new/test1/1ojyD03.csv');

%% load the helper function

mdl = sprMDL(proteinARGs_1, 2);

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

%%
testProteinARGs = cell(0);
% testProteinARGs{end+1} = GenerateProteinARG(145,235,'protein/SH2/test/1g83.csv');
% testProteinARGs{end+1} = GenerateProteinARG(123,210,'protein/SH2/test/1lcj.csv');
% testProteinARGs{end+1} = GenerateProteinARG(140,230,'protein/SH2/test/2hck.csv');
% testProteinARGs{end+1} = GenerateProteinARG(140,230,'protein/SH2/test/2ptk.csv');
% testProteinARGs{end+1} = GenerateProteinARG(21,110,'protein/SH2/test/4u1p.csv');

% testProteinARGs{end+1} = GenerateProteinARG(50,160,'protein/CH1/test/2r0o.csv');
% testProteinARGs{end+1} = GenerateProteinARG(35,145,'protein/CH1/test/5a36.csv');
% testProteinARGs{end+1} = GenerateProteinARG(35,145,'protein/CH1/test/5a4b.csv');
% testProteinARGs{end+1} = GenerateProteinARG(31,136,'protein/CH1/test/2eyn.csv');
% testProteinARGs{end+1} = GenerateProteinARG(45,150,'protein/CH1/test/3lue.csv');

test_result = zeros([1,length(testProteinARGs)]);
test_score = zeros([1,length(testProteinARGs)]);

count = 1;

for i = 1:length(testProteinARGs)
    tmpARG = testProteinARGs{i}.copy;
    [result, score] = mdl.checkPattern(tmpARG);
    test_result(count) = result;
    test_score(count) = score;
    count = count + 1;
end

test_detect_rate = sum(test_result)/length(test_result)

%% compare reuslt
fig = figure;
hax = axes;
binWidth = 5;

histogram(original_score,'binWidth',binWidth)
hold on
histogram(partition_score,'binWidth',binWidth)
hold on
histogram(random_score,'binWidth',binWidth)
hold on
histogram(test_score,'binWidth',binWidth)
hold on
line([mdl.thredshold_score mdl.thredshold_score],get(hax,'YLim'),'Color','k','LineWidth', 2)
hold on
legend('Original Training Protein','Partition Reversed Original Training Protein','Random Protein','Untrained Protein with the Same Domain')
title('SH2 Domain Model Similarity Score Test')
xlabel('Model Similarity Score h(G|S)');
ylabel('Occurence')
set(findall(fig,'-property','FontSize'),'FontSize',22)

%%
% Y = [original_score;
%     partition_score;
%     test_score;
%     random_score(1:length(original_score))];
% 
% figure;
% bar(Y,'FaceColor',[0/255,173/255,76/255])
% 
% % legend('1AD5','1AOU','1FBZ','1QCF','4D8K','Location','northeast')
% set(findall(gca,'-property','FontSize'),'FontSize',22)
% set(gca,'xticklabel',{'Original','Distorted','Test','Random'})
% ylabel('Motif Detection Score')
% title('CH1 Motif Detection Score')
% ylim([0 370])

%%
X = 1:length(original_score);
len = length(X);
bar(X,original_score,'grouped','FaceColor',[255/255,102/255,178/255])
hold on

X = 1:length(partition_score);
X = X+len;
len = len+length(X);
bar(X,partition_score,'grouped','FaceColor',[0/255,204/255,204/255])
hold on

X = 1:length(test_score);
X = X+len;
len = len+length(X);
bar(X,test_score,'grouped','FaceColor',[255/255,255/255,0/255])
hold on

X = 1:length(random_score);
X = [1];
X = X+len;
len = len+length(X);
bar(X,mean(random_score),'grouped','FaceColor',[0/255,102/255,204/255])

xlim([0,len+1])
% legend('Original','Distorted','Test','Random')
ylim([200 400])
% ylim([200 750])
ylabel('Motif Detection Score')
title('SH2 Motif Detection Score')
set(findall(gca,'-property','FontSize'),'FontSize',32)
set(gca,'xticklabel',{[]})
set(gca,'xtick',[])

hold on
y_line = mdl.thredshold_score;
plot(xlim,[y_line y_line],'Color','k','LineWidth', 2,'LineStyle',':')


%% graph matching
G1 = GenerateProteinARG(26,140,'protein/CH1/1sjj.csv');
G2 = GenerateProteinARG(42,150,'protein/CH1/1tjt.csv');

[ match_matrix, C_n, C_e ] = graph_matching(G1, mdl_ARG(G2), false);


%% visual   
visual(26, 140, 32, 137, 'protein/CH1/1sjj.csv', 42, 150, 45, 150, 'protein/CH1/1tjt.csv', match_matrix)
% visual(140, 230, 148, 230, 'protein/SH2/1ad5.csv', 1, 90, 7, 89, 'protein/SH2/1aou.csv', match_matrix)

%%

P1 = getpdb('1sjj');
P2 = getpdb('2eyi');

%%

pdbsuperpose(P1, P2);
h3 = findobj('Tag', 'BioinfoMolviewer'); % retrieve handle for molviewer
evalrasmolscript(h3, ['select all; zoom 200; center selected']);
evalrasmolscript(h3, ['select all; cartoons off; ' ...
                      'select model = 1; strands on; color red; ' ...% ubiquitin
                      'select model = 2; strands on; color blue;']); % SUMO

%%
start = [0.5,0.25,0.25];
trans = [0.5,0.25,0.25;
         0.3,0.7,0;
         0.3,0,0.7;];
     
     
start*trans^8