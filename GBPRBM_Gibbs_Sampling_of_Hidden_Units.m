function [h, P_h_m1_v, P_h_p1_v] = GBPRBM_Gibbs_Sampling_of_Hidden_Units(opt, v)
% Inputs:
% Name       Size    Explanation  
% opt.H:   integer, Number of hidden units
% opt.W:   (V x H), Weight matrix
% opt.z_v: (V x 1), Standard deviations for the visible units
% opt.b_h: (H x 1), Hidden biases
% v:       (V x 1), Visible vector

% Outputs
% Name        Size    Explanation  
% P_h_0_v:   (H x 1), Probability(h=0|v)
% P_h_1_v:   (H x 1), Probability(h=1|v)
% h:         (H x 1), Sampled hidden units vector

H = opt.H;
W = opt.W;
z_v = opt.z_v;
b_h = opt.b_h;

% Probability(h =-1|v)
P_h_m1_v = exp(-(W'*(v./exp(z_v)) + b_h)) ./ (2*cosh(W'*(v./exp(z_v)) + b_h)); 
% Probability(h = 1|v)
P_h_p1_v = exp( (W'*(v./exp(z_v)) + b_h)) ./ (2*cosh(W'*(v./exp(z_v)) + b_h)); 
%P_h_p1_v = 1 - P_h_m1_v;

% Sample hidden units
r = rand(H,1);
h = 2*(r > P_h_m1_v) - 1 ;