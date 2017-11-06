% loda the data

load('match_visual_1sjj_1tjt.mat')

%%
P1 = getpdb('1sjj'); % we match 26-140 and domain is 32-137
P2 = getpdb('1tjt'); % we match 42-150 and domain is 45-150
P3 = getpdb('4d1e');
P4 = getpdb('1wku');
P5 = getpdb('2eyi');


%% Generate Matching

G1 = GenerateProteinARG(26,140,'protein/CH1/1sjj.csv');
G2 = GenerateProteinARG(42,150,'protein/CH1/1tjt.csv');
[ match_matrix, C_n, C_e ] = graph_matching(G1, mdl_ARG(G2), false);

%% Setup Molviwer for pdbsuperpose of special area

[Dist, RMSD, Transf, PBD2TX] = pdbsuperpose(P1, P2, nan, 'Segment',{'26-140:A','42-150:A'});


h = findobj('Tag', 'BioinfoMolviewer');
evalrasmolscript(h(1), [
    'select all; cartoons off; trace on; '... 
    'select model = 1 ; color red;'... 
    'select model = 2 ; color green;'...
    'center selected; '
    ])


%% Setup Molviwer for pdbsuperpose of domain 

[Dist, RMSD, Transf, PBD2TX] = pdbsuperpose(P1, P2, nan, 'Segment',{'26-140:A','42-150:A'});

h = findobj('Tag', 'BioinfoMolviewer');
evalrasmolscript(h(1), [
    'select all; cartoons off; trace on; '... 
    'select model = 1 ; color red;'... 
    'select model = 2 ; color green;'...
    'center selected; '
    ])

evalrasmolscript(h(1), 'display (model = 2 AND 44-150:A) OR (model = 1 AND 31-137:A);');
evalrasmolscript(h(1), 'select model = 1 AND 31:A; color yellow;');
evalrasmolscript(h(1), 'select model = 2 AND 44:A; color yellow;');

%%

[Dist, RMSD, Transf, PBD2TX] = pdbsuperpose(P1, P2, nan,'display','off', 'Segment',{'31-137:A','44-150:A'});
[Dist, RMSD, Transf, PBD3TX] = pdbsuperpose(P1, P3, nan,'display','off', 'Segment',{'31-137:A','37-143:A'});
[Dist, RMSD, Transf, PBD4TX] = pdbsuperpose(P1, P4, nan,'display','off', 'Segment',{'31-137:A','44-150:A'});
[Dist, RMSD, Transf, PBD5TX] = pdbsuperpose(P1, P5, nan,'display','off', 'Segment',{'31-137:A','30-136:A'});

%% Multiple Alignment
modelStruct1 = P1.Model(1);
modelStruct2 = PBD2TX.Model(1);
modelStruct3 = PBD3TX.Model(1);
modelStruct4 = PBD4TX.Model(1);
modelStruct5 = PBD5TX.Model(1);
pdbStruct1 = P1;
pdbStruct2 = PBD2TX;
pdbStruct3 = PBD3TX;
pdbStruct4 = PBD4TX;
pdbStruct5 = PBD5TX;

%=== make sure that the structures are similar
[modelStruct1,modelStruct2]=fixfield(modelStruct1, modelStruct2);
[modelStruct1,modelStruct3]=fixfield(modelStruct1, modelStruct3);
[modelStruct2,modelStruct3]=fixfield(modelStruct2, modelStruct3);
[modelStruct1,modelStruct4]=fixfield(modelStruct1, modelStruct4);
[modelStruct2,modelStruct4]=fixfield(modelStruct2, modelStruct4);
[modelStruct3,modelStruct4]=fixfield(modelStruct3, modelStruct4);
[modelStruct1,modelStruct5]=fixfield(modelStruct1, modelStruct5);
[modelStruct2,modelStruct5]=fixfield(modelStruct2, modelStruct5);
[modelStruct3,modelStruct5]=fixfield(modelStruct3, modelStruct5);
[modelStruct4,modelStruct5]=fixfield(modelStruct4, modelStruct5);

%=== write PDB struct for combo
pdbCombo = rmfield(pdbStruct1, 'Model'); % in case there are many models
pdbCombo.Model(1) = modelStruct1;

f1 = fieldnames(modelStruct1);
pdbCombo.Model(2) = orderfields(modelStruct2, f1);
pdbCombo.Model(3) = orderfields(modelStruct3, f1);
pdbCombo.Model(4) = orderfields(modelStruct4, f1);
pdbCombo.Model(5) = orderfields(modelStruct5, f1);

pdbCombo = orderfields(pdbCombo, pdbStruct1);

if isfield(pdbCombo.Model(1), 'MDLSerNo')
    pdbCombo.Model(1).MDLSerNo = 1;
    pdbCombo.Model(2).MDLSerNo = 2;
    pdbCombo.Model(3).MDLSerNo = 3;
    pdbCombo.Model(4).MDLSerNo = 4;
    pdbCombo.Model(5).MDLSerNo = 5;
end

h = molviewer(pdbCombo);

evalrasmolscript(h, 'select all; wireframe off; spacefill off;');
evalrasmolscript(h, 'color chain; cartoon on;');

evalrasmolscript(h, [
    'select all; cartoons off; trace on; '... 
    'select model = 1 ; color red;'... 
    'select model = 2 ; color green;'...
    'select model = 3 ; color pink;'...
    'select model = 4 ; color cyan;'...
    'select model = 5 ; color grey;'...
    'display (model = 5 AND 30-136:A) OR(model = 4 AND 44-150:A) OR (model = 3 AND 37-143:A) OR (model = 2 AND 44-150:A) OR (model = 1 AND 31-137:A);'...
    'center selected; '
    ])

evalrasmolscript(h(1), 'select model = 1 AND 31:A; color yellow;');
evalrasmolscript(h(1), 'select model = 2 AND 44:A; color yellow;');
evalrasmolscript(h(1), 'select model = 3 AND 37:A; color yellow;');
evalrasmolscript(h(1), 'select model = 4 AND 44:A; color yellow;');
evalrasmolscript(h(1), 'select model = 5 AND 30:A; color yellow;');

%% Setup our own protein for our own visual

P1TX = P1;
P2TX = PBD2TX;
P1_ATOM = P1TX.Model.Atom;
P2_ATOM = P2TX.Model.Atom;

%% Generate Protein

P = P1_ATOM;
PT = [];
res_num = P(1).resSeq;
x = [];
y = [];
z = [];
for i=1:length(P)
    rs = P(i).resSeq;
    if res_num ~= rs
        PT(end+1,1)=res_num;
        PT(end,2)=0;
        PT(end,3)=mean(x);
        PT(end,4)=mean(y);
        PT(end,5)=mean(z);
        res_num = rs;
        x = [P(i).X];
        y = [P(i).Y];
        z = [P(i).Z];
    else
        x = [x,P(i).X];
        y = [y,P(i).Y];
        z = [z,P(i).Z];
    end
end

P1_VISUAL = PT;

P = P2_ATOM;
PT = [];
res_num = P(1).resSeq;
x = [];
y = [];
z = [];
for i=1:length(P)
    rs = P(i).resSeq;
    if res_num ~= rs
        PT(end+1,1)=res_num;
        PT(end,2)=0;
        PT(end,3)=mean(x);
        PT(end,4)=mean(y);
        PT(end,5)=mean(z);
        res_num = rs;
        x = [P(i).X];
        y = [P(i).Y];
        z = [P(i).Z];
    else
        x = [x,P(i).X];
        y = [y,P(i).Y];
        z = [z,P(i).Z];
    end
end

P2_VISUAL = PT;

%% visual

visual(26, 140, 32, 137, P1_VISUAL, 42, 150, 45, 150, P2_VISUAL, match)


