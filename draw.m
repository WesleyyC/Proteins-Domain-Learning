%%
hold on
for i=1:115
    x1=protein(i,3);
    y1=protein(i,4);
    z1=protein(i,5);
    x2=new(i,1);
    y2=new(i,2);
    z2=new(i,3);
    X=[x1,x2]';
    Y=[y1,y2]';
    Z=[z1,z2]';
    plot3(X,Y,Z);
    hold on
end

%%