a=zeros(mdl.mdl_ARGs{1,1}.num_nodes,20);

for i = 1:mdl.mdl_ARGs{1,1}.num_nodes
    a(i,:)=mdl.mdl_ARGs{1,2}.nodes_vector{1,i};
end

imshow(a)


%%
test=find(mdl.mdl_ARGs{1,1}.edges_matrix)
plot(mdl.mdl_ARGs{1,1}.edges_matrix(test), sqrt(mdl.mdl_ARGs{1,1}.edges_cov(test)), '.');


%%
plot(mdl.distance, mdl.node_ab_match.*mdl.node_cd_match.*mdl.sample_match.*mdl.edge_match, '.');

%%
plot(mdl.distance(inspect_index), mdl.edge_match(inspect_index)+k(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), min(mdl.node_ab_match(inspect_index).*mdl.node_cd_match(inspect_index).*mdl.sample_match(inspect_index)*3,1)+k(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), min(mdl.node_ab_match(inspect_index).*mdl.node_cd_match(inspect_index)*1.7,1)+k(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), mdl.node_ab_match(inspect_index).*mdl.node_cd_match(inspect_index).*mdl.sample_match(inspect_index).*mdl.distance(inspect_index), '.');

%%
plot(mdl.distance(inspect_index), mdl.edge_match(inspect_index).*mdl.distance(inspect_index), '.');
