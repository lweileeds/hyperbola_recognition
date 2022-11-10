function [a, b, xc,yc, a_ini, b_ini, xc_ini, yc_ini, time] = G_N_hyperbola_fitting_v2(x, y, win)
% function [a, b, xc,yc] = Gaussian_Newton_hyperbola_fitting(x,y) is used to
% fit a hyperbola to the given noisy points

% x is a vector containing the x coordinates of the given points 
% y is a vector containing the y coordinates of the given points
% win is the length of the window when doing initialization
% (xc,yc) is the center of the hyperbola.
% a is the length of semi-major axis
% b is the length of semi-imaginary axis

tic
if size(x,1) ~= 1
    x = x';            %x should be a row vector
end

if size(y,1) ~= 1
    y = y';            %y should be a row vector
end

xy = [x;y];

% Rotation matrix
r_matr = [0 1; -1 0];

% Initialize the hyperbola
[a, b, xc, yc, ~] = initialize_hyperbola_v2(x, y, win);

a_ini = a;
b_ini = b;
xc_ini = xc;
yc_ini = yc;
xyc = [xc; yc];
delta_v = ones(4,1);

% Fitting
count = 0;
count_max = 51;
while norm(delta_v)>1e-3 && count<count_max
    count=count+1;

if count == count_max
    a = -10;
    break
end

XY = r_matr*(xy-repmat(xyc,1,size(xy,2)));


XY0 = zeros(size(XY));

for i=1:size(XY,2)
    if abs(XY(1,i)) <= a
        XY0(:,i) = [sign(XY(1,i))*a;0];
    elseif abs(XY(1,i))>a && a>=b
        XY0(:,i) = [XY(1,i), sign(XY(2,i))*b*sqrt(XY(1,i)^2-a^2)/a];
    else
        XY0(:,i) = [sign(XY(1,i))*a*sqrt(XY(2,i)^2+b^2)/b, XY(2,i)];
    end
    
end   

for i=1:size(XY,2)
    
    delta_x = [1;0];
    record_max = 101;
    record = 0;
    
    while norm(delta_x) > 1e-4
        record = record+1;
        if record > record_max
            break
        end
    q = [-b^2*XY0(1,i) a^2*XY0(2,i); (a^2+b^2)*XY0(2,i)-b^2*XY(2,i) (a^2+b^2)*XY0(1,i)-a^2*XY(1,i)];
    if sum(sum(isnan(q)))>0 || sum(sum(isinf(q)))>0
        break
    end    
    f(1, 1) = 0.5*(a^2*XY0(2,i)^2-b^2*XY0(1,i)^2+a^2*b^2);
    f(2, 1) = b^2*XY0(1,i)*(XY0(2,i)-XY(2,i))+a^2*XY0(2,i)*(XY0(1,i)-XY(1,i));
    delta_x = pinv(q)*(-f);
    XY0(:,i) = XY0(:,i) + 0.5*delta_x;
    end
end

xy1 = r_matr'*XY0 + repmat(xyc, 1, size(XY0,2));
xy2 = xy - xy1;
xy2 = xy2(:);

J = zeros(2*size(xy, 2), 4);
B = zeros(2, 4);
for j = 1:size(xy, 2)
    B(1, 1) = -a^2*XY0(2, j);
    B(2,1) = a^2*(XY(1, j) - XY0(1, j));
    B(1,2) = -b^2*XY0(1, j);
    B(2,2) = -b^2*(XY(2, j) - XY0(2, j));
    B(1,3) = -a*(b^2 + XY0(2, j)^2);
    B(2,3) = 2*a*XY0(2, j)*(XY(1, j)-XY0(1, j));
    B(1,4) = -b*(a^2 - XY0(1, j)^2);
    B(2,4) = 2*b*XY0(1,j)*(XY(2,j)-XY0(2,j));
    q=[-b^2*XY0(1,j) a^2*XY0(2,j); (a^2+b^2)*XY0(2,j)-b^2*XY(2,j) (a^2+b^2)*XY0(1,j)-a^2*XY(1,j)];
    if sum(sum(isnan(q)))>0 || sum(sum(isinf(q)))>0
        break
    end
    J(2*j-1:2*j,:)=r_matr'*pinv(q)*B;
end

delta_v = pinv(J)*xy2;
if count<1000
    xyc = xyc + 0.3*delta_v(1:2);
    a = a + 0.3*delta_v(3);
    b = b + 0.3*delta_v(4);
    if ~isreal(a) || ~isreal(b) || a<0 || b<0 || isnan(a) || isnan(b) || a==inf || b==inf
        break
    end
elseif count>=1000 && count <1250
    xyc = xyc + 0.3*delta_v(1:2);
    a = a + 0.3*delta_v(3);
    b = b + 0.3*delta_v(4);
    if ~isreal(a) || ~isreal(b) || a<0 || b<0 || isnan(a) || isnan(b) || a==inf || b==inf
        break
    end
else
    xyc = xyc + 0.5*delta_v(1:2);
    a = a + 0.5*delta_v(3);
    b = b + 0.5*delta_v(4);
    if ~isreal(a) || ~isreal(b) || a<0 || b<0 || isnan(a) || isnan(b) || a==inf || b==inf
        break
    end
end
end

xc = xyc(1);
yc = xyc(2);

time = toc;

