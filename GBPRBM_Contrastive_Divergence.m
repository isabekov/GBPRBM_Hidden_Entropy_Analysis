function [out, v_cd] = GBPRBM_Contrastive_Divergence(opt)
% Inputs:
% Name       Size    Explanation  
% opt.k:   integer, Number of steps "k" in Gibbs sampling 
% opt.V:   integer, Number of visible units
% opt.H:   integer, Number of hidden units
% opt.v_0: (V x 1), Initial visible vector "v_0"
% opt.W:   (V x H), Weight matrix of size 
% opt.z_v: (V x 1), Standard deviations for the visible units
% opt.b_v: (V x 1), Visible biases 
% opt.b_h: (H x 1), Hidden biases

% Outputs
% Name        Size    Explanation  
% out.d_W:   (V x H), Contrastive Divergence update for weight matrix W
% out.d_z_v: (V x 1), Contrastive Divergence update for visible variances vector
% out.d_b_v: (V x 1), Contrastive Divergence update for visible biases vector
% out.d_b_h: (H x 1), Contrastive Divergence update for hidden biases vector
%            Gibbs sampling
%        ======================
%	     h_0   h_1	h_(k-1)
%	      *     *     *
%	     / \   / \   * \
%	    /   \ /   \ *   \
%	   *     *     *     *
%	 v_0    v_1   v_2    v_k
v_0 = opt.v_0;
k = opt.CD_order;
V = opt.V;
H = opt.H;
W = opt.W;
b_h = opt.b_h;
b_v = opt.b_v;
z_v = opt.z_v;

% Visible vector in contrastive divergence
v_cd = zeros(V,k+1);
% Hidden vector in contrastive divergence
h_cd = zeros(H,k);
v_cd(:,0 +1) = v_0;
% Gibbs sampling
for t = 0:k-1             
   % Sample hidden units
   [h_cd(:,t+1), ~, ~] = GBPRBM_Gibbs_Sampling_of_Hidden_Units(opt, v_cd(:,t+1));
   % Sample visible units by taking just the mean of the multivariate
   % Gaussian. Do not add extra noise.
   [v_cd(:,t+2), ~]    = GBPRBM_Gibbs_Sampling_of_Visible_Units(opt, h_cd(:,t+1));   
end
v_k = v_cd(:,k+1);
% Hyperbolic tangent terms which will be used in the CD-updates
tanh_v_0 = tanh(W'*(v_0./exp(z_v)) + b_h);
tanh_v_k = tanh(W'*(v_k./exp(z_v)) + b_h);
% Contrastive divergence updates
out.d_W =  repmat(1./exp(z_v),1,H).*(v_0*tanh_v_0' - v_k*tanh_v_k');
out.d_b_v =  1./exp(z_v).*(v_0 - v_k);
out.d_b_h = tanh_v_0 - tanh_v_k;
out.d_z_v = 1./exp(z_v).*(v_0.^2/2 - v_0.*b_v - v_0.*(W*tanh_v_0) -...
                         (v_k.^2/2 - v_k.*b_v - v_k.*(W*tanh_v_k)));                                        