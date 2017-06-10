function MNIST_Plot_Pairs_of_Images(Img_Orig,Img_Recon, varargin)
% This function plots pairs of MNIST images from two matrices in a grid
% by reshaping 784x1 column vectors from supplied matrices into pairs of 
% two 28x28 images concatenated side by side.
%                 Inputs 
% ============================================================
% Img_Orig : 784 x N_Plots    First image
% Img_Recon: 784 x N_Plots    Second image
% varargin : string           Figure and axes title
%                 Output
% ============================================================
% A figure with pairs of images plotted using imagesc function

if sum(size(Img_Orig)==size(Img_Recon)) ~= 2
    error('Dimension mismatch');
end
if size(Img_Orig,1) ~= 784
   error('Weight matrix does not have 784 rows!'); 
end
if length(varargin) >= 1
   str = varargin{1};
else
   str = 'MNIST: Pairs of images';
end
% Number of images to reconstruct and plot
N_Plots = size(Img_Orig,2);

N_Horiz = floor(sqrt(N_Plots));
N_Vert  = ceil(sqrt(N_Plots));
if N_Horiz*N_Vert < N_Plots
   N_Horiz = N_Horiz + 1; 
end
% Vertical padding
V_Pad = 2; % pixels
% Horizontal padding
H_Pad = 2; % pixels
% Number of pixels in width (height as well, this is a square image)
N = 28; % pixels
% Allocate memory for the compined image
A = ones(N*N_Horiz + H_Pad*(N_Horiz-1), 2*N*N_Vert + V_Pad*(N_Vert-1));
% Counter for the image sample to be plotted
cnt = 1;
for i=1:N_Horiz
    for j=1:N_Vert   
        % Calculating image position on the grid
        idx_H1 = (i-1)*(N+H_Pad)+1;
        idx_H2 =     i*(N+H_Pad)-H_Pad;
        idx_V1 = (j-1)*(2*N+V_Pad)+1;
        idx_V2 =     j*(2*N+V_Pad)-V_Pad;
        % Reshape a 784x1 vector into a 28x28 image. Do the same operation
        % for the reconstructed image
        A(idx_H1:idx_H2, idx_V1:idx_V2) = [reshape(Img_Orig(:,cnt),N,N) reshape(Img_Recon(:,cnt),N,N)];
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