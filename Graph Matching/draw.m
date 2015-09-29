function [] = draw(p1,p2,match,score,matchStart,matchEnd)

    % pre-processed
    if(size(p1,1)<size(p2,1))
        tmp=p1;
        p1=p2;
        p2=tmp;
        match=match';
    end
    
    % shift to the center
    p1(:,3)=p1(:,3)-mean(p1(:,3));
    p1(:,4)=p1(:,4)-mean(p1(:,4));
    p1(:,5)=p1(:,5)-mean(p1(:,5));
    p2(:,3)=p2(:,3)-mean(p2(:,3))+max(p1(:,3))*1.5;
    p2(:,4)=p2(:,4)-mean(p2(:,4))+max(p1(:,4))*1.5;
    p2(:,5)=p2(:,5)-mean(p2(:,5))+max(p1(:,5))*1.5;
    
    % draw
    figure()
    for i=1:size(p1,1)-1
        x1=p1(i,3);
        y1=p1(i,4);
        z1=p1(i,5);
        x2=p1(i+1,3);
        y2=p1(i+1,4);
        z2=p1(i+1,5);
        X=[x1,x2]';
        Y=[y1,y2]';
        Z=[z1,z2]';
        plot3(X,Y,Z,'--r','LineWidth',0.1);
        hold on
    end
    
    for i=1:size(p2,1)-1
        x1=p2(i,3);
        y1=p2(i,4);
        z1=p2(i,5);
        x2=p2(i+1,3);
        y2=p2(i+1,4);
        z2=p2(i+1,5);
        X=[x1,x2]';
        Y=[y1,y2]';
        Z=[z1,z2]';
        plot3(X,Y,Z,'--g','LineWidth',0.1);
        hold on
    end
    
    if nargin~=6
        matchStart = 1;
        matchEnd = size(p1,1);
    end
    
    for i=matchStart:matchEnd
        x1=p1(i,3);
        y1=p1(i,4);
        z1=p1(i,5);
        p2i = find(match(i,:));
        color_depth = score(i, p2i);
        % round color_depth to 0.33/0.66/0.99
        color_depth = ceil(3*color_depth)/3;
        if ~isempty(p2i)
            x2=p2(p2i,3);
            y2=p2(p2i,4);
            z2=p2(p2i,5);
            X=[x1,x2]';
            Y=[y1,y2]';
            Z=[z1,z2]';
            plot3([x1,x1],[y1,y1],[z1,z1],'.r','MarkerSize',20);
            plot3([x2,x2],[y2,y2],[z2,z2],'.g','MarkerSize',20);
            plot3(X,Y,Z,':','LineWidth',0.8,'Color',[color_depth color_depth color_depth]);
            hold on
        end
    end
   
    set(gca,'Color',[0 0 0]);

end
