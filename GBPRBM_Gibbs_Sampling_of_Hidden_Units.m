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

% Probability(h_j =-1|v) for all j
P_h_m1_v = sigmf(-2*(W'*(v./exp(z_v)) + b_h), [1 0]); 

% Probability(h_j =+1|v) for all j
P_h_p1_v = 1 - P_h_m1_v;

% Sample hidden units, example for a single unit
%  _____________________________________________
% 0                    0.5                      1
%  _____________________________________________
% |        Pr(h=-1) + Pr(h=+1) = Area = 1       |
%  _____________________________________________
% |       A random number ~U(h) falls here      |
% |============================-----------------| 
% |________ P_h_m1_v _________||___ P_h_p1_v ___|
% |         Event "-1"        ||    Event "1"   |
%  _____________________________________________
% If the random number is greater than P_h_m1_v, i.e. 
% rand > P_h_m1_v = 1, then decide that h=+1
% If the random number is smaller than P_h_m1_v, i.e. 
% rand > P_h_m1_v = 0, then decide that h=-1
%       Mapping: 
% Logical "0" -> "-1"
% Logical "1" -> "+1" 
% i.e. f(x) = 2*x -1
r = rand(H,1);
h = 2*(r > P_h_m1_v) - 1 ;
