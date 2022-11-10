function out_v = cluster_cleaning(xx, yy)
% function out_v=cluster_cleaning(xx, yy) is used to clean the detected
% clusters with large overlapping (such as 85 percent)

len = size(xx,1);
out_v = ones(len,1);

for i = 1:len-1
    if out_v(i) == 1
        upper_b = min(i+5, len);
        num_1 = length(xx{i, 1});
        for j = i+1:upper_b
            num_2 = length(xx{j, 1});
            num_3 = find_num_same_points(xx{i,1}, yy{i,1},xx{j,1}, yy{j,1});
            if num_3>4*num_2/5 && num_1>5*num_2/4 
                out_v(j) = 0;
            elseif num_3>4*num_1/5 && num_2>5*num_1/4
                out_v(i) = 0;
                continue;
            elseif num_3>4*num_1/5 && num_3>4*num_2/5 && num_1<=num_2
                out_v(j) = 0;
                elseif num_3>4*num_1/5 && num_3>4*num_2/5 && num_1>num_2
                    out_v(i) = 0;
                    continue;
            end
        end
    end
end        