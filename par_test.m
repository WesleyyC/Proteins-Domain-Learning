


a = cell(1,1000);
tic()
parfor i = 1:1000
    frequency=0;
    for j = 1:100000
        frequency = frequency + 1;
    end
    a{i} = frequency/i;
end
toc()