%% Script for matching the protein
tic();
draw_flag = 0;

%% Set Up Edge and Node Attribute

BLOSUM_Sigma = 5;   % node attribute

distance_cutoff = 13;   % edge attribute cut off

save = 0;

%% Set Up Protein

proteinOneFile = 'new_1EAZ_2_99.csv';
start_sequence_one = 1;
end_sequence_one = 100;

proteinTwoFile = 'new_1U2B_5_121.csv';
start_sequence_two = 5;
end_sequence_two = 121;

% proteinOneFile = 'new_4D1E_CH.csv';
% start_sequence_one = 50;
% end_sequence_one = 150;
% 
% proteinTwoFile = 'new_4Q59_CH.csv';
% start_sequence_two = 50;
% end_sequence_two = 150;

[proteinOneARG,p1] = GenerateProteinARGs(start_sequence_one,end_sequence_one, proteinOneFile,distance_cutoff);
[proteinTwoARG,p2] = GenerateProteinARGs(start_sequence_two,end_sequence_two, proteinTwoFile,distance_cutoff);


%% Generate BLOSSUM

C=[9,-1,-1,-3,0,-3,-3,-3,-4,-3,-3,-3,-3,-1,-1,-1,-1,-2,-2,-2];
S=[-1,4,1,-1,1,0,1,0,0,0,-1,-1,0,-1,-2,-2,-2,-2,-2,-3];
T=[-1,1,4,1,-1,1,0,1,0,0,0,-1,0,-1,-2,-2,-2,-2,-2,-3];
P=[-3,-1,1,7,-1,-2,-1,-1,-1,-1,-2,-2,-1,-2,-3,-3,-2,-4,-3,-4];
A=[0,1,-1,-1,4,0,-1,-2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-2,-3];
G=[-3,0,1,-2,0,6,-2,-1,-2,-2,-2,-2,-2,-3,-4,-4,0,-3,-3,-2];
N=[-3,1,0,-2,-2,0,6,1,0,0,-1,0,0,-2,-3,-3,-3,-3,-2,-4];
D=[-3,0,1,-1,-2,-1,1,6,2,0,-1,-2,-1,-3,-3,-4,-3,-3,-3,-4];
E=[-4,0,0,-1,-1,-2,0,2,5,2,0,0,1,-2,-3,-3,-3,-3,-2,-3];
Q=[-3,0,0,-1,-1,-2,0,0,2,5,0,1,1,0,-3,-2,-2,-3,-1,-2];
H=[-3,-1,0,-2,-2,-2,1,1,0,0,8,0,-1,-2,-3,-3,-2,-1,2,-2];
R=[-3,-1,-1,-2,-1,-2,0,-2,0,1,0,5,2,-1,-3,-2,-3,-3,-2,-3];
K=[-3,0,0,-1,-1,-2,0,-1,1,1,-1,2,5,-1,-3,-2,-3,-3,-2,-3];
M=[-1,-1,-1,-2,-1,-3,-2,-3,-2,0,-2,-1,-1,5,1,2,-2,0,-1,-1];
I=[-1,-2,-2,-3,-1,-4,-3,-3,-3,-3,-3,-3,-3,1,4,2,1,0,-1,-3];
L=[-1,-2,-2,-3,-1,-4,-3,-4,-3,-2,-3,-2,-2,2,2,4,3,0,-1,-2];
V=[-1,-2,-2,-2,0,-3,-3,-3,-2,-2,-3,-3,-2,1,3,1,4,-1,-1,-3];
F=[-2,-2,-2,-4,-2,-3,-3,-3,-3,-3,-1,-3,-3,0,0,0,-1,6,3,1];
Y=[-2,-2,-2,-3,-2,-3,-2,-3,-2,-1,2,-2,-2,-1,-1,-1,-1,3,7,2];
W=[-2,-3,-3,-4,-3,-2,-4,-4,-3,-2,-2,-3,-3,-1,-3,-2,-3,1,2,11];

BLOSUM=[C',S',T',P',A',G',N',D',E',Q',H',R',K',M',I',L',V',F',Y',W'];

BLOSUM=exp(BLOSUM/BLOSUM_Sigma);

% normalize the symmetric matrix
% Not perfect
% s=sum(BLOSUM,2);
% n=repmat(s,1,20);
% BLOSUM=BLOSUM./n;

%% Match the prote

[match,score] = graph_matching(proteinOneARG, proteinTwoARG,BLOSUM);
figure; imshow(match);
if(draw_flag)
    draw(p1,p2,match,score);
end

%% Save the result
if (save)
    datetime=datestr(now);
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore

    save (datetime) 
end

toc()