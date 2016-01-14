edges=cell(3);

count = 1;
for a=1:3
    for b=1:3
        edges{a,b}=count;
        count=count+1;
    end
end

len = length(edges);
flat = cell (1, len*len);
for flat_p = 1:len*len
    flat{flat_p}=edges{floor((flat_p-1)/len)+1,flat_p-(floor((flat_p-1)/len))*len};
end