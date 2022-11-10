function [xx, yy, xxx, yyy, time] = column_connection_clustering_v2( ...
    real_im, edge_thre, len_cs, len_se, ratio)
% This function is used to detect hyperbola regions in an image.  
% real_im is the input image
% edge_thre is the threshold value used in detecting edges from the
% input image with Canny method
% len_cs is the minimum number of elements in a column segment
% len_se is the minimum number of same elements between two column segments

% xx contains the x-coordinates of the clusters
% yy contains the y-coordinates of the clusters
% xxx contains the x-coordinates of the central string of clusters
% yyy contains the y-coordinates of the central string of clusters


tic
% Detect edges
s_real_im = conv2(real_im, 1/25*ones(5,5),'same');
mea = mean(s_real_im,2);
s_real_im = s_real_im - repmat(mea, 1, size(s_real_im, 2));
ed = edge(s_real_im, 'canny', edge_thre);
ed_im = ed.*s_real_im;

% Identify bright areas
max_v = max(max(ed_im));
ind = ed_im>max_v*ratio;
thre = sum(sum(s_real_im.*ind))/sum(sum(ind));
im_thre = s_real_im>thre;

% Scan the thresholded image to find clusters
[xx, yy, xxx, yyy]=scan_columns_ltor(im_thre, len_cs, len_se);

xx_1 = [];
yy_1 = [];
xxx_1 = [];
yyy_1 = [];
for i = 1:size(xxx,1)    
    integ = 0;
    if length(xxx{i,1})>11
        for j = 2:length(xxx{i,1})
            integ = integ + sqrt((yyy{i,1}(j) - yyy{i,1}(j-1))^2 + 1);
        end
        if length(xx{i,1})/integ^2 < 0.3
            num_row = size(xx_1,1);
            xx_1{num_row+1,1} = xx{i,1};
            yy_1{num_row+1,1} = yy{i,1};
            xxx_1{num_row+1,1} = xxx{i,1};
            yyy_1{num_row+1,1} = yyy{i,1};
        end
    end
end

xx = [];
yy = [];
xxx = [];
yyy = [];
for i = 1:size(xxx_1,1)
    num_row = size(xx,1);
    xx{num_row+1,1} = xx_1{i,1};
    yy{num_row+1,1} = yy_1{i,1};
    xxx{num_row+1,1} = xxx_1{i,1};
    yyy{num_row+1,1} = yyy_1{i,1};
    yyys = smooth(yyy_1{i,1});
    grad = diff(yyys);
    angles = zeros(size(grad));
    IND = abs(grad)<=1e-6;
    angles(IND) = pi/2;
    IND = find(abs(grad)>1e-6);
    for j = 1:length(IND)
        if grad(IND(j))<0
            angles(IND(j)) = -atan(1/grad(IND(j)));
        else
            angles(IND(j)) = pi-atan(1/grad(IND(j)));
        end
    end
    d_angles = angles(2:end) - angles(1:end-1);
    IND = find(d_angles<-pi/18);
    if ~isempty(IND)
        for k = 1:length(IND)+1
            if k == 1
                INDX = xx_1{i,1}<xxx_1{i,1}(IND(k))+1;

                num_row = size(xx,1);
                xx{num_row+1,1} = xx_1{i,1}(INDX);
                yy{num_row+1,1} = yy_1{i,1}(INDX);
                xxx{num_row+1,1} = xxx_1{i,1}(1:IND(k));
                yyy{num_row+1,1} = yyy_1{i,1}(1:IND(k));
            elseif k>1 && k<length(IND)+1
                INDX = xx_1{i,1}>xxx_1{i,1}(IND(k-1)) & xx_1{i,1}<xxx_1{i,1}(IND(k))+1;
                
                num_row = size(xx,1);
                xx{num_row+1,1} = xx_1{i,1}(INDX);
                yy{num_row+1,1} = yy_1{i,1}(INDX);
                xxx{num_row+1,1} = xxx_1{i,1}(IND(k-1)+1:IND(k));
                yyy{num_row+1,1} = yyy_1{i,1}(IND(k-1)+1:IND(k));
            else
                INDX = xx_1{i,1}>xxx_1{i,1}(IND(k-1));
                num_row = size(xx,1);
                xx{num_row+1,1} = xx_1{i,1}(INDX);
                yy{num_row+1,1} = yy_1{i,1}(INDX);
                xxx{num_row+1,1} = xxx_1{i,1}(IND(k-1)+1:end);
                yyy{num_row+1,1} = yyy_1{i,1}(IND(k-1)+1:end);
            end
        end
    end
end
xx_1 = [];
yy_1 = [];
xxx_1 = [];
yyy_1 = [];
for i = 1:size(xxx,1)
    if length(xxx{i,1})>11
        integ = sum(sqrt((yyy{i,1}(2:end)-yyy{i,1}(1:end-1)).^2+ones(1,length(yyy{i,1})-1)));
        
        if length(xx{i,1})/integ^2 < 0.3 && min(yyy{i,1})~=yyy{i,1}(1) && min(yyy{i,1})~=yyy{i,1}(end)
            num_row = size(xx_1,1);
            xx_1{num_row+1,1} = xx{i,1};
            yy_1{num_row+1,1} = yy{i,1};
            xxx_1{num_row+1,1} = xxx{i,1};
            yyy_1{num_row+1,1} = yyy{i,1};
        elseif length(xx{i,1})/integ^2 < 0.3 && min(yyy{i,1})==yyy{i,1}(1)
            yyys = smooth(yyy{i,1});
            if yyys(1)+yyys(3)-2*yyys(2)<0
                num_row = size(xx_1,1);
                xx_1{num_row+1,1} = xx{i,1};
                yy_1{num_row+1,1} = yy{i,1};
                xxx_1{num_row+1,1} = xxx{i,1};
                yyy_1{num_row+1,1} = yyy{i,1};
            end
        elseif length(xx{i,1})/integ^2 < 0.3 && min(yyy{i,1})==yyy{i,1}(end)
            yyys = smooth(yyy{i,1});
            if yyys(end)+yyys(end-2)-2*yyys(end-1)<0
                num_row = size(xx_1,1);
                xx_1{num_row+1,1} = xx{i,1};
                yy_1{num_row+1,1} = yy{i,1};
                xxx_1{num_row+1,1} = xxx{i,1};
                yyy_1{num_row+1,1} = yyy{i,1};
            end
        end
    end
end
xx = xx_1;
yy = yy_1;
xxx = xxx_1;
yyy = yyy_1;
cleaned_im = zeros(size(real_im));
for i = 1:size(xx,1)
    for j = 1:length(xx{i})
        cleaned_im(yy{i}(j),xx{i}(j)) = 1;
    end
end
time = toc;

% Display the input image
subplot(2,3,1);
imagesc(real_im)
colormap gray(256);
title('Input Image');

% Display processed images
subplot(2, 3, 2);
imagesc(s_real_im)
colormap gray(256);
title('After pre-processing');

subplot(2, 3, 3);
imagesc(im_thre);
title('After Thresholding');

subplot(2,3,4);
imagesc(cleaned_im);
title('Output of C3 algorithm');