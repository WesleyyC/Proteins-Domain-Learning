gamma = 0.23;
epsilon = 0.15;
tau = 0.01;

trans = [   0, 1-2*gamma-tau, gamma, gamma, tau;
            0, 1-2*gamma-tau, gamma, gamma, tau;
            0, 1-epsilon-tau, epsilon, 0, tau;
            0, 1-epsilon-tau, 0, epsilon, tau;
            0, 0, 0, 0, 1];
        
        
sample_no = 100000;
chain_length = 1000;
starting_value =1;
total_distance = zeros(1, sample_no);
insert_distance = zeros(1, sample_no);

for j=1:sample_no
    chain = zeros(1,chain_length);
    chain(1)=starting_value;
    for i=2:chain_length
        this_step_distribution = trans(chain(i-1),:);
        cumulative_distribution = cumsum(this_step_distribution);
        r = rand();
        chain(i) = find(cumulative_distribution>r,1);
    end
    total_distance(j)=sum(chain<5);
    insert_distance(j)=sum(chain==3);
end

%%

p = insert_distance./total_distance;

%%
figure
h1 = histfit(total_distance,100,'exponential');
hold on
h2 = histfit(insert_distance./0.1769,100,'exponential');

%%
start = [0, 1-2*gamma-tau, gamma, gamma];
start = start/sum(start);
g_trans = trans(1:4,1:4)/(1-tau);

start*g_trans^15

%%
pn=zeros(1,chain_length);
pi=zeros(1,chain_length);
for x = 1:chain_length
    pi(x)=tau*(1-tau)^x/0.1769;
    pn(x)=tau*(1-tau)^x;
end
plot(1:chain_length, pn)
hold on
plot(1:chain_length, pi)
