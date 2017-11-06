%% same sfam
% 
proteinARGs = cell(0);
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1aipD03.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1g3iU02.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1g4aF02.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1hjpA03.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1ixrB03.csv');

proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1kyiU02.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1mn3A00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1nv9A01.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1oaiA00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1oksA00.csv');

proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1oqyA02.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1otrA00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1pgyA00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1t6oA00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test9/1tr8B02.csv');

mdl_test = sprMDL(proteinARGs, 4);

%% sfam

proteinARGs_sfam = cell(0);

proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/1tteA02.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/1v92A00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/1wj7A01.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/1xb2B01.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2dnaA00.csv');

proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2dzlA00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2ekkA00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2g3qA00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2l2dA00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2l4fA00.csv');

proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2lbcA02.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2lvaA01.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/2mj5B00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/3e21A00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test9/3ihpB05.csv');
%%

sfam_score = zeros([1,length(proteinARGs_sfam)]);

for i = 1:length(proteinARGs_sfam)
    [result,score] = mdl_test.checkPattern(proteinARGs_sfam{i});
    sfam_score(i) = score;
end

mean(sfam_score)


%% other

proteinARGs_other = cell(0);

proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1e9rA02.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1gvnC00.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1j09A04.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1jb0F00.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1mu5A02.csv');

proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1no1A00.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1oaoC01.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1oizA01.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1pjqA03.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/1v33A02.csv');

proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/2cruA01.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/2gnoA02.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/2j5yA00.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/3futA02.csv');
proteinARGs_other{end+1} = GenerateProteinARG('test2-protein/test10/4qosA02.csv');

%%

other_score = zeros([1,length(proteinARGs_other)]);

for i = 1:length(proteinARGs_other)
    [result,score] = mdl_test.checkPattern(proteinARGs_other{i});
    other_score(i) = score;
end

mean(other_score)


%%
figure;
X = 1:length(sfam_score);
len = length(X);
bar(X,sfam_score,'grouped','FaceColor',[255/255,102/255,178/255])
hold on

X = 1:length(other_score);
X = X+len;
len = len+length(X);
bar(X,other_score,'grouped','FaceColor',[0/255,204/255,204/255])
hold on

legend('Same Superfamily','Other Superfamily')
ylabel('Motif Detection Score')
title('2.10.25.10 Motif Detection Score')
set(findall(gca,'-property','FontSize'),'FontSize',32)
set(gca,'xticklabel',{[]})
set(gca,'xtick',[])
