function num = find_num_same_elements(a,b)
% function num = find_num_same_elements(a,b) is used to find the number of
% same elements in vectors a and b.
% a is a vector with elements different to each ther oand in ascending order.
% b has the same properties with a.

len_a = length(a);
len_b = length(b);
counter = 0;
j0 = 0;
for i = 1:len_a
    j = j0 + 1;
   while (j<=len_b && a(i)~=b(j))
        j = j+1;
    end
    if j <= len_b
        counter = counter+1;
        j0 = j;
    end
end
num = counter;
        
    