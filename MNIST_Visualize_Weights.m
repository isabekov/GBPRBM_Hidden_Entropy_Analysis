function MNIST_Visualize_Weights(W, varargin)
% This function plots a sequence of MNIST images from matrix W in a grid 
% by reshaping 784x1 column vectors into 28x28 images concatenated side by side.
%                 Inputs      
% ======================================================
%       W  : 784 x N_Plots    A sequence of MNIST images
% varargin : string           Figure and axes title
%                 Output
% ======================================================
% A figure with images plotted using imagesc function

% Number of images to reconstruct and plot
[V, N_Plots] = size(W);
if V ~= 784
   error('Weight matrix does not have 784 rows!'); 
end
if length(varargin) >= 1
   str = varargin{1};
else
   str = 'MNIST: Images';
end
N_Horiz = floor(sqrt(N_Plots));
N_Vert  =  ceil(sqrt(N_Plots));
if N_Horiz*N_Vert < N_Plots
   N_Horiz = N_Horiz + 1; 
end
% Vertical padding
V_Pad = 1; % pixel
% Horizontal padding
H_Pad = 1; % pixel
% Number of pixels in width (height as well, this is a square image)
N = 28; % pixels
% Allocate memory for the compined image
A = ones(N*N_Horiz + H_Pad*(N_Horiz-1), N*N_Vert + V_Pad*(N_Vert-1));
% Counter for the image sample to be plotted
cnt = 1;
for i=1:N_Horiz
    for j=1:N_Vert         
        % Calculate image position on the grid
        idx_H1 = (i-1)*(N+H_Pad)+1;
        idx_H2 =     i*(N+H_Pad)-H_Pad;
        idx_V1 = (j-1)*(N+V_Pad)+1;
        idx_V2 =     j*(N+V_Pad)-V_Pad;
        % Reshape a 784x1 vector into a 28x28 image. Do the same operation
        % for the reconstructed image
        A(idx_H1:idx_H2, idx_V1:idx_V2) = reshape(W(:,cnt),N,N);
        cnt = cnt + 1;
        if cnt > N_Plots
           break 
        end
    end    
end
figure('Name', str, 'NumberTitle', 'Off');
imagesc(A);
set(gca, 'DataAspectRatio', [1 1 1]);
axis off;
colormap gray;    
title(str);