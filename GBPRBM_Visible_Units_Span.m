function [opt] = GBPRBM_Visible_Units_Span(opt)
% ======= Determination of the range of  X-,Y- and Z-axes ======
% Determine the range for X-axis (v1 dimension), Y-axis (v2 dimension)
% Y-axis (v3 dimension) of the boundary values of the centroids. 

V = opt.V;
H = opt.H;
W = opt.W;
if size(W,1)~=V
    error('Matrix W is not of size VxH. Quitting.')
end
b_v = opt.b_v;       
VRatio = opt.VRatio;

% Each cell Centroid{j} contains a matrix of size (V x 2^j) with the
% centroid coordinates.
Centroid = cell(H,1);
% Each cell Hid_Vec{j} contains a matrix of size (H x 2^j) with 
% the truth table. This is needed to label each centroid with its hidden
% configuration.
Hid_Vec  = cell(H,1);
for j=1:H
    % "Centroid" and "Hid_Vec" are explained above.
    Centroid{j} = zeros(V,2^j);
    Hid_Vec{j}  = zeros(j,2^j);
    if j==1
        % If j==1, then start walking away from the visible bias point "b_v".
        % Take a step of size |W(:,j)| in the direction of W(:,j)
        % away from the visible bias point "b_v", i.e. step "forward"
        Centroid{j}(:,1) = b_v + W(:,j);
        % Label: "+1"
        Hid_Vec{j}(1,1)  = +1;
        % Take a step of size |W(:,j)| in the direction of -W(:,j)
        % away from the visible bias point "b_v", i.e. step "backward"
        Centroid{j}(:,2) = b_v - W(:,j);
        % Label: "-1"
        Hid_Vec{j}(1,2)  = -1;        
    else
        % After the walk to the previous level (j-1) of propagation,
        % there are 2^(j-1) centroids.
        for k=1:2^(j-1)
            % From previous centroid "k" walk by W(:,j) forward to
            % obtain new level (j) centroid coordinates which will be
            % stored at odd index "2*(k-1)+1"
            Centroid{j}(:,2*(k-1)+1) = Centroid{j-1}(:,k) + W(:,j);
            % Append to the truth table of the previous level (j-1) the new
            % hidden neuron configuration "+1"
            Hid_Vec{j}(:,2*(k-1)+1)  = [Hid_Vec{j-1}(:,k); +1];
            % From previous centroid "k" walk by W(:,j) backward to
            % obtain new level (j) centroid coordinates which will be
            % stored at even index "2*k"
            Centroid{j}(:,2*k)       = Centroid{j-1}(:,k) - W(:,j);
            % Append to the truth table of the previous level (j-1) the new
            % hidden neuron configuration "-1"
            Hid_Vec{j}(:,2*k)  = [Hid_Vec{j-1}(:,k); -1];
        end
    end
end

opt.Centroid = Centroid;
opt.Hid_Vec = Hid_Vec;
opt.v = struct;
for i=1:V
    opt.v(i).min = min(Centroid{end}(i,:));
    opt.v(i).max = max(Centroid{end}(i,:));
    opt.v(i).range = opt.v(i).max - opt.v(i).min;
    opt.v(i).AxesLimMax = opt.v(i).max + VRatio(i)*opt.v(i).range;
    opt.v(i).AxesLimMin = opt.v(i).min - VRatio(i)*opt.v(i).range;
end