function [] = visual(start_sequence_1, end_sequence_1, start_highlight_sequence_1, end_highlight_sequence_1, fileName_1, start_sequence_2, end_sequence_2, start_highlight_sequence_2, end_highlight_sequence_2, fileName_2, match_result)
%function [] = visual(p1,p2,match,score,matchStart,matchEnd,highlight)
    
    protein_1 = csvread(fileName_1);
    protein_2 = csvread(fileName_2);

    start_idx_1 = NaN;
    end_idx_1 = NaN;
    start_highlight_idx_1 = NaN;
    end_highlight_idx_1 = NaN;
    start_idx_2 = NaN;
    end_idx_2 = NaN;
    start_highlight_idx_2 = NaN;
    end_highlight_idx_2 = NaN;

    % extract index
    for i = 1:size(protein_1,1)
        if start_sequence_1 == protein_1(i,1) && isnan(start_idx_1)
            start_idx_1 = i;
        end
        if end_sequence_1==protein_1(i,1) && isnan(end_idx_1)
            end_idx_1= i;
        end
        if start_highlight_sequence_1 == protein_1(i,1) && isnan(start_highlight_idx_1)
            start_highlight_idx_1 = i;
        end
        if end_highlight_sequence_1==protein_1(i,1) && isnan(end_highlight_idx_1)
            end_highlight_idx_1= i;
        end
    end
    start_highlight_idx_1 = start_highlight_idx_1-start_idx_1+1;
    end_highlight_idx_1 = end_highlight_idx_1-start_idx_1+1;
    for i = 1:size(protein_2,1)
        if start_sequence_2 == protein_2(i,1) && isnan(start_idx_2)
            start_idx_2 = i;
        end
        if end_sequence_2==protein_2(i,1) && isnan(end_idx_2)
            end_idx_2= i;
        end
        if start_highlight_sequence_2 == protein_2(i,1) && isnan(start_highlight_idx_2)
            start_highlight_idx_2 = i;
        end
        if end_highlight_sequence_2==protein_2(i,1) && isnan(end_highlight_idx_2)
            end_highlight_idx_2= i;
        end
        
    end
    start_highlight_idx_2 = start_highlight_idx_2-start_idx_2+1;
    end_highlight_idx_2 = end_highlight_idx_2-start_idx_2+1;
    
    protein_1 = protein_1(start_idx_1:end_idx_1,:);
    protein_2 = protein_2(start_idx_2:end_idx_2,:);
    
    score=match_result(:,1:end-1);
    match=zeros(size(score));
    for i = 1:size(score,2)
        [~,idx]=max(match_result(:,i));
        match(idx,i)=1;
    end
    
    % shift to the center
    protein_1(:,3)=protein_1(:,3)-mean(protein_1(:,3));
    protein_1(:,4)=protein_1(:,4)-mean(protein_1(:,4));
    protein_1(:,5)=protein_1(:,5)-mean(protein_1(:,5));
    protein_2(:,3)=protein_2(:,3)-mean(protein_2(:,3))+max(protein_1(:,3))*1.5;
    protein_2(:,4)=protein_2(:,4)-mean(protein_2(:,4))+max(protein_1(:,4))*1.5;
    protein_2(:,5)=protein_2(:,5)-mean(protein_2(:,5))+max(protein_1(:,5))*1.5;
    
    % draw
    figure()
    
    for i=1:size(protein_1,1)-1
        x1=protein_1(i,3);
        y1=protein_1(i,4);
        z1=protein_1(i,5);
        x2=protein_1(i+1,3);
        y2=protein_1(i+1,4);
        z2=protein_1(i+1,5);
        X=[x1,x2]';
        Y=[y1,y2]';
        Z=[z1,z2]';
        if start_highlight_idx_1<=i && end_highlight_idx_1>i
            plot3(X,Y,Z,'--b','LineWidth',0.1);
        else
            plot3(X,Y,Z,'--r','LineWidth',0.1);
        end
        hold on
    end
    
    for i=1:size(protein_2,1)-1
        x1=protein_2(i,3);
        y1=protein_2(i,4);
        z1=protein_2(i,5);
        x2=protein_2(i+1,3);
        y2=protein_2(i+1,4);
        z2=protein_2(i+1,5);
        X=[x1,x2]';
        Y=[y1,y2]';
        Z=[z1,z2]';
        if start_highlight_idx_2<=i && end_highlight_idx_2>i
            plot3(X,Y,Z,'--b','LineWidth',0.1);
        else
            plot3(X,Y,Z,'--g','LineWidth',0.1);
        end
        hold on
    end
    
    for i=1:size(match,1)
        x1=protein_1(i,3);
        y1=protein_1(i,4);
        z1=protein_1(i,5);
        p2i = find(match(i,:));
        color_depth = score(i, p2i);
        % round color_depth to 0.33/0.66/0.99
        color_depth = ceil(3*color_depth)/3;
        if ~isempty(p2i)
            x2=protein_2(p2i,3);
            y2=protein_2(p2i,4);
            z2=protein_2(p2i,5);
            X=[x1,x2]';
            Y=[y1,y2]';
            Z=[z1,z2]';
            plot3([x1,x1],[y1,y1],[z1,z1],'.r','MarkerSize',10);
            plot3([x2,x2],[y2,y2],[z2,z2],'.g','MarkerSize',10);
            plot3(X,Y,Z,':','LineWidth',0.8,'Color',[color_depth color_depth color_depth]);
            hold on
        end   
    end
    
    set(gca,'Color',[0 0 0]);

end