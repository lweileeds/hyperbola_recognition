function [ncc_1, ncc_2, o_x, o_y] = ncc_values_v2(x, y)
% This function calculates the normalized cross correlation (NCC) value of
% the parameters of hyperbolae.

% The input x, y are row vectors. They are the x and y coordinates of the
% detected central strings.

% ncc_1 is the normalized cross correlation value of the first derivatives 
% of y with respect to x against the first derivatives from a template.

% ncc_2 is the normalized cross correlation value of the second derivatives 
% of y with respect to x against the second derivatives from a template.

% o_x is the truncated version of x and o_y is the corresponding truncated
% version of y.

% since the elements of y are recorded as the the row number of related 
% pixels, we use the negative values of them to calculate the corresponding 
% normalized cross correlation values
y = -y;

% Smooth y
y = smooth(y,3);

% Truncate the vector y to avoid the subeffect of convolution
len = length(y);
dy = diff(y);
dy = smooth(dy,3);
ind = find(y==max(y),1,'first');
dy_plus = dy(find(dy>0));
dy_minus = dy(find(dy<0));

if ~isempty(dy_plus)
    mean_plus = mean(dy_plus);
else
    mean_plus = inf;
end

if ~isempty(dy_minus)
    mean_minus = mean(dy_minus);
else
    mean_minus = -inf;
end

[dy_max, loc_max] = max(dy);
[dy_min, loc_min] = min(dy);

if dy_max/mean_plus>4 && dy_min/mean_minus<=4
    if loc_max<=ind && loc_max<len/2 && len-loc_max>11
        y = y(loc_max+1:end);
        x = x(loc_max+1:end);
    elseif loc_max>ind && loc_max>=1/2*len && loc_max>11
        y = y(1:loc_max);
        x = x(1:loc_max);
    end
end
if dy_max/mean_plus<=4 && dy_min/mean_minus>4
    if loc_min<=ind && loc_min<len/2 && len-loc_min>11
        y = y(loc_min+1:end);
        x = x(loc_min+1:end);
    elseif loc_min>ind && loc_min>=1/2*len && loc_min>11
        y = y(1:loc_min);
        x = x(1:loc_min);
    end
end

if dy_max/mean_plus>4 && dy_min/mean_minus>4
    if loc_max<=ind && loc_min<=ind
        max_v = max(loc_max,loc_min);
        if max_v<len/2 && len-max_v>11
            y = y(max_v+1:end);
            x = x(max_v+1:end);
        end
    elseif loc_max>ind && loc_min>ind
        min_v = min(loc_max,loc_min);
        if min_v>=1/2*len && min_v>11
            y = y(1:min_v);
            x = x(1:min_v);
        end
    else
        max_v = max(loc_max,loc_min);
        min_v = min(loc_max,loc_min);
        if min_v<1/2*len && max_v>1/2*len && max_v-min_v>11
            y = y(min_v+1:max_v);
            x = x(min_v+1:max_v);
        end
    end
end

len=length(y);

if rem(length(ind),2) == 0
    vertex_loc = ind(length(ind)/2);
else
    vertex_loc = ind((length(ind)+1)/2);
end

% Generating the x and y coordinates of template 
t_x = 1-vertex_loc:len-vertex_loc;
t_y = -5*sqrt(1+t_x.^2/16);
d_v = diff(t_y);
dd_v = diff(d_v);

d_y = diff(y);
d_y = smooth(d_y);
dd_y = diff(d_y);
dd_y = smooth(dd_y);

ncc_1 = d_v*d_y/(norm(d_v)*norm(d_y));
ncc_2 = dd_v*dd_y/(norm(dd_y)*norm(dd_v));
o_x = x;
o_y = y;