%% This is an example to demonstrate how to use the C3 hyperbola 
% detection and fitting algorithm 

addpath('c3_algorithms/')  

% load a GPR image
real_im = imread('img1.png');

% detect and fit hyperbolae using the proposed C3 algorithm;
% the ouput is a list of the coordinates of fitted hyperbolae
hyperbolae = c3_hyperbola_fitting(real_im); 
