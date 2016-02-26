a=zeros(mdl.mdl_ARGs{1,1}.num_nodes,20);

for i = 1:mdl.mdl_ARGs{1,1}.num_nodes
    a(i,:)=mdl.mdl_ARGs{1,2}.nodes_vector{1,i};
end

imshow(a)

%%
test=find(mdl.mdl_ARGs{1,1}.edges_matrix)
plot(mdl.mdl_ARGs{1,1}.edges_matrix(test), sqrt(mdl.mdl_ARGs{1,1}.edges_cov(test)), '.');

%%
inspect_index = 1:100:length(mdl.distance);

%%
plot(mdl.distance, mdl.node_ab_match.*mdl.node_cd_match.*mdl.sample_match.*mdl.edge_match, '.');

%%
plot(mdl.distance(inspect_index), mdl.edge_match(inspect_index).*mdl.distance(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), mdl.node_ab_match(inspect_index).*mdl.node_cd_match(inspect_index).*mdl.sample_match(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), min(mdl.node_ab_match(inspect_index).*mdl.node_cd_match(inspect_index).*mdl.sample_match(inspect_index), mdl.edge_match(inspect_index)*2), '.');

%%
plot(mdl.distance(inspect_index), min(mdl.node_ab_match(inspect_index).*mdl.node_cd_match(inspect_index)*1.7,1)+k(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), mdl.node_ab_match(inspect_index).*mdl.node_cd_match(inspect_index).*mdl.sample_match(inspect_index).*mdl.distance(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), mdl.edge_match(inspect_index).*mdl.distance(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), min(mdl.node_ab_match(inspect_index).*mdl.node_cd_match(inspect_index).*mdl.sample_match(inspect_index),0.06).*mdl.distance(inspect_index), '.');


%%
plot(mdl.distance(inspect_index),mdl.node_ab_match(inspect_index),'.')



%%
x = 0:0.1:12;
scale = 5;
shift = 1;
% scale = 1;
% shift = 3;

a_sigmf=@(v)p_sigmf(v,scale,shift);
y = (1-arrayfun(a_sigmf,x));

plot(x,y)
%%
scale = 1;
shift = 3;
a_sigmf=@(v)p_sigmf(v,scale,shift);
plot(mdl.distance(inspect_index),1-arrayfun(a_sigmf,mdl.distance(inspect_index)),'.');

%%
n=1:10;

lambda=@(n)nchoosek(2*n,n);
prob=arrayfun(lambda,n) .* 0.05.^(n);
plot(n,prob)

%%
