function [v_mean, v_samp] = GBPRBM_Gibbs_Sampling_of_Visible_Units(opt, h)
% Inputs:
% Name       Size    Explanation  
% opt.V:   integer, Number of visible units
% opt.W:   (V x H), Weight matrix
% opt.z_v: (V x 1), Standard deviations for the visible units
% opt.b_v: (V x 1), Visible biases 
% h:       (H x 1), Hidden vector

% Outputs
% Name        Size    Explanation  
% v_mean:   (V x 1), Mean of the estimated signal
% v_samp:   (V x 1), Sample of the signal (v_mean + Gaussian noise)

W = opt.W;
z_v = opt.z_v;
b_v = opt.b_v;
V = opt.V;

v_mean = (b_v + W*h);
v_samp = v_mean + sqrt(exp(z_v)).*randn(V,1);
                  