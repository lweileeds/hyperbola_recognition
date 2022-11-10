function num = find_num_same_points(xx1, yy1, xx2, yy2)
% function num = find_num_same_elements(a,b) is used to find the number of
% same points from two group of points.
% xx1 contains the x coordinates of the points from group one
% yy1 contains the y coordinates of the points from group one
% xx2 contains the x coordinates of the points from group two
% yy2 contains the y coordinates of the points from group two

len_1 = length(xx1);
len_2 = length(xx2);
num = 0;
if len_1<len_2
    for i = 1:len_1
        for j = 1:len_2
            if xx1(i)==xx2(j) && yy1(i)==yy2(j)
                num = num+1;
                continue;
            end
        end
    end
else
    for i = 1:len_2
        for j = 1:len_1
            if xx2(i)==xx1(j) && yy2(i)==yy1(j)
                num = num+1;
                continue;
            end
        end
    end
end
