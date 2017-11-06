%% random

random = rand(6,7);
random = normr(random).*normc(random);

partial_correct = random;
partial_correct(3,2) = 2;
partial_correct(4,3) = 2;
partial_correct = normr(partial_correct).*normc(partial_correct);

correct = zeros(6,7);
correct(1,7)=1;
correct(2,1)=1;
correct(3,2)=1;
correct(4,3)=1;
correct(5,4)=1;
correct(6,5)=1;
correct(6,6)=1;


fsize = 500;
RI = imref2d([6,7]);

figure;
subplot(1,3,1)
imshow(random, RI,'InitialMagnification',fsize);
subplot(1,3,2)
imshow(partial_correct, RI,'InitialMagnification',fsize);
subplot(1,3,3)
imshow(correct, RI,'InitialMagnification',fsize);

%%

init = eye(50);

not_clear = init;
not_clear(1,:)=0;
not_clear(1,3)=1;
not_clear(2,:)=0;
not_clear(2,1)=1;
not_clear(3,:)=0;
not_clear(3,5)=1;
not_clear(4,:)=0;
not_clear(4,2)=1;
not_clear(5,:)=0;
not_clear(5,4)=1;
not_clear(50,:)=0;
not_clear(50,48)=1;
not_clear(49,:)=0;
not_clear(49,50)=1;
not_clear(48,:)=0;
not_clear(48,46)=1;
not_clear(47,:)=0;
not_clear(47,49)=1;
not_clear(46,:)=0;
not_clear(46,47)=1;

correct = init;
correct(1:5,1:5)=0;
correct(1:5,end)=1;
correct(46:50,46:50)=0;
correct(end,46:50)=1;

fsize = 500;
RI = imref2d([50,50]);

figure;
subplot(1,2,1)
imshow(not_clear,RI,'InitialMagnification',fsize);
subplot(1,2,2)
imshow(correct,RI,'InitialMagnification',fsize);

%%
range = -10:0.1:10;

z = zeros(length(range),length(range));

xcount = 1;
for x = -10:0.1:10
    ycount = 1;
    for y = -10:0.1:10
        z(xcount,ycount)=x^2-y^2;
        ycount = ycount+1;
    end
    xcount = xcount+1;
end

surf(range, range, z);colorbar

%%

histogram(y, 'BinWidth',1)
xlim([0,50])
xticks(0:4:50)
xlabel('Distance(Å)')
ylabel('Number of Instances')
title('Sample Distance between all Amino Acids')
set(findall(gca,'-property','FontSize'),'FontSize',32)

