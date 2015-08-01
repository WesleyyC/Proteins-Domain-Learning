%% A script for generating proteins ARG


%% Set Up Variable

BLOSUM_Sigma = 2;

distance_cutoff = 9;

start_sequence = 1;

end_sequence = 115;


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
s=sum(BLOSUM,2);
n=repmat(s,1,20);
BLOSUM=BLOSUM./n;

%% Input the File

protein = csvread('new_4D1E.csv');
protein = protein(start_sequence:end_sequence,:);
number_of_AA = size(protein,1);

%% Ready for ARGs

proteinStructure = zeros(number_of_AA);

proteinAtrs = zeros([number_of_AA,1]);

for i = 1:number_of_AA
    
    proteinAtrs(i) = protein(i,2);
    
    for j = i+1:number_of_AA
        x1=protein(i,3);
        y1=protein(i,4);
        z1=protein(i,5);
        x2=protein(j,3);
        y2=protein(j,4);
        z2=protein(j,5);
        dist = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
        dist = dist*(dist<=distance_cutoff);
        proteinStructure(i,j)=dist;
        proteinStructure(j,i)=dist;
    end
end

