function [xx, yy, xxx, yyy]=scan_columns_ltor(im_thre, len_cs, len_se)
% This function scans the thresholded images column by column to find the 
% clusters.
% len_cs is the minimun number of elements in a column segment
% len_se is the minimum number of same elements between two column segments


[y,x] = find(im_thre);

min_x = min(x);
max_x = max(x);

y_array = [];
for i = min_x:max_x
    ind = find(x==i);
    num_ind = length(ind);
    
    if i == min_x
        yy_1 = y(ind);
        diff_yy = diff(yy_1);
        ind_2 = find(diff_yy>1);
        if ~isempty(ind_2)
            diff_ind_2 = diff([0 ind_2' num_ind]);
            ind_3 = find(diff_ind_2>=len_cs);
            num_ind_3 = length(ind_3);
            
            xxx = cell(num_ind_3,1);
            yyy = cell(num_ind_3,1);
            xx = cell(num_ind_3,1);
            yy = cell(num_ind_3,1);
            y_array = cell(num_ind_3,1);
            counter = 0;
            for j = 1:length(ind_2)+1
                if (j==1 && ind_2(1)>=len_cs)
                    counter = counter+1;
                    xxx{counter, 1} = [xxx{counter, 1} i];
                    xx{counter, 1} = [xx{counter, 1} i*ones(1,ind_2(1))];
                    yyy{counter, 1}=[yyy{counter, 1} mean(yy_1(1:ind_2(1)))];
                    yy{counter, 1} = [yy{counter, 1} yy_1(1:ind_2(1))'];
                    y_array{counter, 1} = yy_1(1:ind_2(1));
                elseif j>1 && j<=length(ind_2)
                    if ind_2(j)-ind_2(j-1) >= len_cs
                        counter = counter+1;
                        xxx{counter,1} = [xxx{counter,1} i];
                        xx{counter,1} = [xx{counter,1} i*ones(1,ind_2(j)-ind_2(j-1))];
                        yyy{counter,1} = [yyy{counter,1} mean(yy_1(ind_2(j-1)+1:ind_2(j)))];
                        yy{counter,1} = [yy{counter,1} yy_1(ind_2(j-1)+1:ind_2(j))'];
                        y_array{counter,1} = yy_1(ind_2(j-1)+1:ind_2(j));
                    end
                elseif j == length(ind_2)+1
                    if length(yy_1)-ind_2(j-1) >= len_cs
                        counter = counter+1;
                        xxx{counter,1} = [xxx{counter,1} i];
                        xx{counter,1} = [xx{counter,1} i*ones(1,length(yy_1)-ind_2(j-1))];
                        yyy{counter,1} = [yyy{counter,1} mean(yy_1(ind_2(j-1)+1:end))];
                        yy{counter,1} = [yy{counter,1} (yy_1(ind_2(j-1)+1:end))'];
                        y_array{counter,1} = yy_1(ind_2(j-1)+1:end);
                    end
                end
            end
            
        else
            if length(yy_1) >= len_cs
                xxx{1, 1} = i;
                xx{1, 1} = i*ones(1,length(yy_1));
                yyy{1, 1} = mean(yy_1);
                yy{1, 1} = yy_1';
                y_array{1, 1} = yy_1;
            end
        end
        
    else
        if num_ind < len_cs
            for ii = 1:size(y_array,1)
                y_array{ii,1} = [];
            end
        else
            y_array_1 = y_array;
            vec_count = zeros(size(y_array,1),1);
            yy_1 = y(ind);
            diff_yy = diff(yy_1);
            ind_2 = find(diff_yy>1);
            
            if isempty(ind_2) 
                for k = 1:size(y_array_1,1)
                    if length(yy_1)>=len_cs && find_num_same_elements(yy_1,y_array_1{k,1})>=len_se
                        xxx{k,1} = [xxx{k,1} i];
                        yyy{k,1} = [yyy{k,1} mean(yy_1)];
                        xx{k,1} = [xx{k,1} i*ones(1,length(yy_1))];
                        yy{k,1} = [yy{k,1} yy_1'];
                        y_array{k,1} = yy_1;
                        vec_count(k,1) = 1;
                    end
                end
                if sum(vec_count) == 0
                    if length(yy_1) >= len_cs
                        xxx{size(y_array,1)+1, 1} = i;
                        yyy{size(y_array,1)+1, 1} = mean(yy_1);
                        xx{size(y_array,1)+1, 1} = i*ones(1,length(yy_1));
                        yy{size(y_array,1)+1, 1} = yy_1';
                        y_array{size(y_array,1)+1, 1} = yy_1;
                    end
                    for m = 1:length(vec_count)
                        y_array{m, 1} = [];
                    end
                else
                    for m = 1:length(vec_count)
                        if vec_count(m) == 0
                            y_array{m, 1} = [];
                        end
                    end
                end                
            else                
                for j = 1:length(ind_2)+1
                    if j==1 && ind_2(1)>=len_cs
                        counter = 0;
                        for k = 1:size(y_array_1, 1)
                            if find_num_same_elements(yy_1(1:ind_2(1)), y_array_1{k,1})>=len_se
                                xxx{k,1} = [xxx{k,1} i];
                                yyy{k,1} = [yyy{k,1} mean(yy_1(1:ind_2(1)))];
                                xx{k,1} = [xx{k,1} i*ones(1,ind_2(1))];
                                yy{k,1} = [yy{k,1} (yy_1(1:ind_2(1)))'];
                                y_array{k,1} = yy_1(1:ind_2(1));
                                vec_count(k,1) = 1;
                                counter = counter+1;
                            end
                        end
                        if counter == 0
                            xxx{size(y_array,1)+1, 1} = i;
                            yyy{size(y_array,1)+1, 1} = mean(yy_1(1:ind_2(1)));
                            xx{size(y_array,1)+1, 1} = i*ones(1,ind_2(1));
                            yy{size(y_array,1)+1, 1} = (yy_1(1:ind_2(1)))';
                            y_array{size(y_array,1)+1, 1} = yy_1(1:ind_2(1));
                        end
                    elseif j>1 && j<=length(ind_2) && ind_2(j)-ind_2(j-1)>=len_cs
                        counter = 0;
                        for k = 1:size(y_array_1,1)
                            if find_num_same_elements(yy_1(ind_2(j-1)+1:ind_2(j)),y_array_1{k,1}) >= len_se
                                if vec_count(k,1) == 0
                                    xxx{k,1} = [xxx{k,1} i];
                                    yyy{k,1} = [yyy{k,1} mean(yy_1(ind_2(j-1)+1:ind_2(j)))];
                                    xx{k,1} = [xx{k,1} i*ones(1,ind_2(j)-ind_2(j-1))];
                                    yy{k,1} = [yy{k,1} (yy_1(ind_2(j-1)+1:ind_2(j)))'];
                                    y_array{k,1} = yy_1(ind_2(j-1)+1:ind_2(j));
                                    vec_count(k,1) = 1;
                                    counter = counter+1;
                                else
                                    kk = size(y_array,1)+1;
                                    xxx{kk,1} = xxx{k,1};
                                    yyy{kk,1} = yyy{k,1};
                                    xx{kk,1} = xx{k,1};
                                    yy{kk,1} = yy{k,1};
                                    indd = xx{kk,1}==i;
                                    xx{kk,1}(indd) = [];
                                    yy{kk,1}(indd) = [];
                                    yyy{kk,1}(end) = mean(yy_1(ind_2(j-1)+1:ind_2(j)));
                                    xx{kk,1} = [xx{kk,1} i*ones(1,ind_2(j)-ind_2(j-1))];
                                    yy{kk,1} = [yy{kk,1} (yy_1(ind_2(j-1)+1:ind_2(j)))'];
                                    y_array{kk,1} = yy_1(ind_2(j-1)+1:ind_2(j));
                                    counter = counter+1;
                                end
                            end
                        end
                        if counter == 0
                            xxx{size(y_array,1)+1,1} = i;
                            yyy{size(y_array,1)+1,1} = mean(yy_1(ind_2(j-1)+1:ind_2(j)));
                            xx{size(y_array,1)+1,1} = i*ones(1,ind_2(j)-ind_2(j-1));
                            yy{size(y_array,1)+1,1} = (yy_1(ind_2(j-1)+1:ind_2(j)))';
                            y_array{size(y_array,1)+1,1} = yy_1(ind_2(j-1)+1:ind_2(j));
                        end
                    elseif j==length(ind_2)+1 && length(yy_1)-ind_2(j-1)>=len_cs
                        counter = 0;
                        for k = 1:size(y_array_1,1)
                            if find_num_same_elements(yy_1(ind_2(j-1)+1:end),y_array_1{k,1})>=len_se
                                if vec_count(k,1) == 0
                                    xxx{k,1} = [xxx{k,1} i];
                                    yyy{k,1} = [yyy{k,1} mean(yy_1(ind_2(j-1)+1:end))];
                                    xx{k,1} = [xx{k,1} i*ones(1,length(yy_1)-ind_2(j-1))];
                                    yy{k,1} = [yy{k,1} (yy_1(ind_2(j-1)+1:end))'];
                                    y_array{k,1} = yy_1(ind_2(j-1)+1:end);
                                    vec_count(k,1) = 1;
                                    counter = counter+1;
                                else
                                    kk = size(y_array,1)+1;
                                    xxx{kk,1} = xxx{k,1};
                                    yyy{kk,1} = yyy{k,1};
                                    yyy{kk,1}(end) = mean(yy_1(ind_2(j-1)+1:end));
                                    xx{kk,1} = xx{k,1};
                                    yy{kk,1} = yy{k,1};
                                    indd = xx{kk,1}==i;
                                    xx{kk,1}(indd) = [];
                                    yy{kk,1}(indd) = [];
                                    xx{kk,1} = [xx{kk,1} i*ones(1,length(yy_1)-ind_2(j-1))];
                                    yy{kk,1} = [yy{kk,1} (yy_1(ind_2(j-1)+1:end))'];
                                    y_array{kk,1} = yy_1(ind_2(j-1)+1:end);
                                    counter = counter+1;
                                end
                            end
                        end
                        
                        if counter==0
                            xxx{size(y_array,1)+1, 1} = i;
                            yyy{size(y_array,1)+1, 1} = mean(yy_1(ind_2(j-1)+1:end));
                            xx{size(y_array,1)+1, 1} = i*ones(1,length(yy_1)-ind_2(j-1));
                            yy{size(y_array,1)+1, 1} = (yy_1(ind_2(j-1)+1:end))';
                            y_array{size(y_array, 1)+1,1} = yy_1(ind_2(j-1)+1:end);
                        end
                    end
                end
                for m = 1:length(vec_count)
                    if vec_count(m)==0
                        y_array{m, 1} = [];
                    end
                end
            end
        end
    end
end