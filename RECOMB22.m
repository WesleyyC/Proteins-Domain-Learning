% %% same sfam
% 
proteinARGs = cell(0);
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/1hx2A00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/1k37A00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/1n7dA08.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/1npeB02.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/1nt0G02.csv');

proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/1tpgA02.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/1uzqA03.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/1yo8A03.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/2nprA02.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/2rhpA02.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/2rnlA00.csv');

proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/3c9aD00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/3fbyB01.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/3gcxE00.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/3poyA06.csv');
proteinARGs{end+1} = GenerateProteinARG('test2-protein/test8/3qcwB03.csv');

mdl_2 = sprMDL(proteinARGs, 1);

%%
mdl_test = mdl_2;

%% sfam

proteinARGs_sfam = cell(0);
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1ataA00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1autL01.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1bf9A00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1g1rC02.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1hj7A02.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1hx2A00.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1k37A00.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1n7dA08.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1npeB02.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1nt0G02.csv');
% 
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1tpgA02.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1uzqA03.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1yo8A03.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/1z6cA00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/2bz6L00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/2k2tA00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/2mgrA01.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/2nprA02.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/2rhpA02.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/2rnlA00.csv');
% 
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/3c9aD00.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/3fbyB01.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/3gcxE00.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/3poyA06.csv');
% proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/3qcwB03.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/3sovA02.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/3u7uL00.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/4bxwB01.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/4cudA03.csv');
proteinARGs_sfam{end+1} = GenerateProteinARG('test2-protein/test8/4o02B05.csv');

%%

sfam_score = zeros([1,length(proteinARGs_sfam)]);

for i = 1:length(proteinARGs_sfam)
    [result,score] = mdl_test.checkPattern(proteinARGs_sfam{i});
    sfam_score(i) = score;
end

mean(sfam_score)


%% topol

proteinARGs_topol = cell(0);
proteinARGs_topol{end+1} = GenerateProteinARG('test2-protein/test4/1kkeA01.csv');
proteinARGs_topol{end+1} = GenerateProteinARG('test2-protein/test4/1lk9A01.csv');
proteinARGs_topol{end+1} = GenerateProteinARG('test2-protein/test4/1qiuC01.csv');
proteinARGs_topol{end+1} = GenerateProteinARG('test2-protein/test4/1v1hA01.csv');
proteinARGs_topol{end+1} = GenerateProteinARG('test2-protein/test4/1v1iB01.csv');
proteinARGs_topol{end+1} = GenerateProteinARG('test2-protein/test4/2horA01.csv');
proteinARGs_topol{end+1} = GenerateProteinARG('test2-protein/test4/2hoxA01.csv');

%%

topol_score = zeros([1,length(proteinARGs_topol)]);

for i = 1:length(proteinARGs_topol)
    [result,score] = mdl_test.checkPattern(proteinARGs_topol{i});
    topol_score(i) = score;
end

mean(topol_score)


%% arch

proteinARGs_arch = cell(0);
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/1bx8A00.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/1c9tH00.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/1cklE01.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/1eakA03.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/1h8pA02.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/1l6jA03.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/1pdcA00.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/1skzA01.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/2o39C01.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/3bg4D00.csv');
proteinARGs_arch{end+1} = GenerateProteinARG('test2-protein/test5/4ayeB01.csv');

%%

arch_score = zeros([1,length(proteinARGs_arch)]);

for i = 1:length(proteinARGs_arch)
    [result,score] = mdl_test.checkPattern(proteinARGs_arch{i});
    arch_score(i) = score;
end

mean(arch_score)

%% class

proteinARGs_class = cell(0);
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/1b8wA00.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/1nh2C00.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/1nvpC00.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/1rm1B02.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/1ytfC00.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/2de6C03.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/2sh1A00.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/3gtqI01.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/3rzdI01.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/4a3bI01.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/4gv5B00.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/4nbbB03.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/4nbcC03.csv');
proteinARGs_class{end+1} = GenerateProteinARG('test2-protein/test6/4nbfB03.csv');

%%

class_score = zeros([1,length(proteinARGs_class)]);

for i = 1:length(proteinARGs_class)
    [result,score] = mdl_test.checkPattern(proteinARGs_class{i});
    class_score(i) = score;
end

mean(class_score)

%%
figure;
X = 1:length(sfam_score);
len = length(X);
bar(X,sfam_score,'grouped','FaceColor',[255/255,102/255,178/255])
hold on

X = 1:length(topol_score);
X = X+len;
len = len+length(X);
bar(X,topol_score,'grouped','FaceColor',[0/255,204/255,204/255])
hold on

X = 1:length(arch_score);
X = X+len;
len = len+length(X);
bar(X,arch_score,'grouped','FaceColor',[255/255,255/255,0/255])
hold on

X = 1:length(class_score);
X = X+len;
len = len+length(X);
bar(X,class_score,'grouped','FaceColor',[0/255,102/255,204/255])
hold on

% 
% y_line = mdl_test.thredshold_score;
% plot(xlim,[y_line y_line],'Color','k','LineWidth', 2,'LineStyle',':')
% y_line = mean(sfam_score);
% plot(xlim,[y_line y_line],'Color',[255/255,102/255,178/255],'LineWidth', 2,'LineStyle',':')
% y_line = mean(topol_score);
% plot(xlim,[y_line y_line],'Color',[0/255,204/255,204/255],'LineWidth', 2,'LineStyle',':')
% y_line = mean(arch_score);
% plot(xlim,[y_line y_line],'Color',[255/255,255/255,0/255],'LineWidth', 2,'LineStyle',':')
% y_line = mean(class_score);
% plot(xlim,[y_line y_line],'Color',[0/255,102/255,204/255],'LineWidth', 2,'LineStyle',':')


legend('Same Superfamily','Same Fold but Different Superfamily','Same Architecture but Different Fold','Same Class but Different Architecture')
ylabel('Motif Detection Score')
title('2.10.25.10 Motif Detection Score')
set(findall(gca,'-property','FontSize'),'FontSize',32)
set(gca,'xticklabel',{[]})
set(gca,'xtick',[])
