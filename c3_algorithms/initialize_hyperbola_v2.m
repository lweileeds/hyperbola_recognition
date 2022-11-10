function [a_i, b_i, x_i, y_i, init_p]=initialize_hyperbola_v2(x,y,win)
% function [a_i, b_i, x_i, y_i, init_p]=initialize_hyperbola(x,y, win) is used to
% compute the parameters of south facing hyperbola based on the given noisy
% points
% x is a vector containing the x coordinates of the given points
% y is a vector containing the y coordinates of the given points
% win is the length of a window along x axis such as 2
% (x_i,y_i) is the center of the hyperbola.
% a_i is the length of semi-major axis
% b_i is the length of semi-imaginary axis
% The final hyperbola should have the following form:
% (y-y_i)^2/a_i^2 - (x-x_i)^2/b_i^2 = 1

[y_max, ind] = max(y);
max_x = max(x);
min_x = min(x);
leng = max_x - min_x + 1;

if win > leng/4
    win = leng/4;
end

y_t = y_max;
x_t = x(ind);

if x_t-min_x < leng/4
    step = (max_x - x_t)/2;
    IND = find(x > x_t+step-win/2 & x < x_t+step+win/2);
    y_ri = mean(y(IND));
    if y_ri > y_t-1
        [my, I] = min(y(IND));
        y_ri = my;
        xs = x(IND);
        x_ri = xs(I);
    else
        x_ri = mean(x(IND));
    end
   
    IND = find(x > x_t+2*step-win/2 & x < x_t+2*step+win/2);
    
    y_le = mean(y(IND));
    if y_le > y_ri-1
        [my, I] = min(y(IND));
        y_le = my;
        xs = x(IND);
        x_le = xs(I);
    else
        x_le = mean(x(IND));
    end
    % Reflect it to the left side
    x_le = 2*x_t-x_le;
elseif x_t-min_x>=leng/4 && max_x-x_t>=leng/4 && x_t-min_x <= max_x-x_t
    IND = find(x > (x_t+2*min_x)/3-win/2 & x < (x_t+2*min_x)/3+win/2);
    
    y_le = mean(y(IND));
    if y_le > y_t-1
        [my, I] = min(y(IND));
        y_le = my;
        xs = x(IND);
        x_le = xs(I);
    else
        x_le = mean(x(IND));
    end
    
    IND = find(x>(2*x_t-x_le+2*max_x)/3-win/2 & x<(2*x_t-x_le+2*max_x)/3+win/2);
    
    y_ri = mean(y(IND));
    if y_ri > y_le-1
        [my, I] = min(y(IND));
        y_ri = my;
        xs = x(IND);
        x_ri = xs(I);
    else
        x_ri = mean(x(IND));
    end
     x_reflect = 2*x_t-x_ri;
    y_reflect = y_ri;
    x_ri = 2*x_t-x_le;
    y_ri = y_le;
    x_le = x_reflect;
    y_le = y_reflect;
elseif x_t-min_x>=leng/4 && max_x-x_t>=leng/4 && x_t-min_x >= max_x-x_t
    IND = find(x>(x_t+2*max_x)/3-win/2 & x<(x_t+2*max_x)/3+win/2);
    
    y_ri = mean(y(IND));
    if y_ri > y_t-1
        [my, I] = min(y(IND));
        y_ri = my;
        xs = x(IND);
        x_ri = xs(I);
    else
        x_ri = mean(x(IND));
    end
    
    IND = find(x>(2*x_t-x_ri+2*min_x)/3-win/2 & x<(2*x_t-x_ri+2*min_x)/3+win/2);
    
    y_le = mean(y(IND));
    if y_le > y_ri-1
        [my, I] = min(y(IND));
        y_le = my;
        xs = x(IND);
        x_le = xs(I);
    else
        x_le = mean(x(IND));
    end
else
    step = (x_t-min_x)/2;
    IND = find(x>x_t-step-win/2 & x<x_t-step+win/2);
    
    y_ri = mean(y(IND));
    if y_ri > y_t-1
        [my,I] = min(y(IND));
        y_ri = my;
        xs = x(IND);
        x_ri = xs(I);
    else
        x_ri = mean(x(IND));
    end
     x_ri = 2*x_t-x_ri;
    IND = find(x>x_t-2*step-win/2 & x<x_t-2*step+win/2);
    
    y_le = mean(y(IND));
    if y_le > y_ri-1
       [my,I] = min(y(IND));
       y_le = my;
       xs = x(IND);
       x_le = xs(I);
    else
        x_le = mean(x(IND));
    end
    
end

x_i = x_t;
r_t = (x_t-x_ri)^2;
l_t = (x_t-x_le)^2;
 
upper_b=(l_t*y_t-r_t*(y_t-y_le))/l_t;   
lower_b=(x_t-x_ri)/(x_t-x_le)*(y_t-y_le)+y_t;
 

if y_ri<=lower_b || y_ri>=upper_b
    y_ri=(lower_b+upper_b)/2;
end

y_i = (l_t*y_ri^2-r_t*y_le^2+(r_t-l_t)*y_t^2)/(2*(y_ri*l_t-y_le*r_t+y_t*(r_t-l_t)));


init_p = [x_le y_le x_t y_t x_ri y_ri];

a_i = abs(y_i-y_t);
b_i = sqrt(l_t*a_i^2/((y_le-y_i)^2-a_i^2));


