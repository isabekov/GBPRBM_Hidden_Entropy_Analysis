function [h_fig, NEEoIHU] = GBPRBM_NEEoIHU(H_Test, varargin)
% This function computes Normalized Empirical Entropy of Individual Hidden
% Units based on activations (i.e. realizations [-1,+1]) of the test hidden 
% vectors by taking a histogram of individual hidden units and calculating
% hidden entropy for each unit.

[h_count, ~] = hist(H_Test',2);
[H, N_Samples_Test] = size(H_Test);
if length(varargin) >= 1
   str = varargin{1};
else
   str = sprintf('GBPRBM, H=%d: Normalized Empirical Entropy of Individual Hidden Units', H);
end

p_h_m = h_count(2,:)/N_Samples_Test;
p_h_p = h_count(1,:)/N_Samples_Test;
Entropy = zeros(H,1);
for j=1:H
    if p_h_m(j)==1
       p_h_m(j)=1-eps; 
    elseif p_h_m(j)==0
       p_h_m(j)=1+eps; 
    end

    if p_h_p(j)==0
       p_h_p(j)=1+eps; 
    elseif p_h_p(j)==1
       p_h_p(j)=1-eps; 
    end
    Entropy(j) = -(p_h_p(j).*log2(p_h_p(j)) + p_h_m(j).*log2(p_h_m(j)));
end
% Normalized Empirical Entropy of Individual Hidden Units
NEEoIHU = sum(Entropy)/H;

h_fig = figure('Name', str, 'NumberTitle', 'Off', ...
       'Units', 'normalized', 'Position', [ 0.1432    0.4648    0.6792    0.4324]); 
   
subplot(131);
imagesc(H_Test);
xlabel('Test samples');
ylabel('Hidden unit index');
if length(varargin) == 2
   RMSE_Test = varargin{2};
   title(sprintf('Activations, RMSE (test) = %0.3f', RMSE_Test));
else
   title('Activations');
end

subplot(132);
bar(p_h_m);
set(gca, 'XLim', [0.5, H+0.5]);
set(gca, 'YLim',[0 1]);
set(gca, 'View', [90 90]);
xlabel('Hidden unit index');
ylabel('p(h_j)');
title('Normalized histogram (Probability)');

subplot(133);
bar(Entropy);
set(gca, 'View', [90 90]);
set(gca, 'XLim', [0.5, H+0.5]);
set(gca, 'YLim',[0 1]);
xlabel('Hidden unit index');
ylabel('Bits');
title(sprintf('Sum of EEoIHUs = %0.2f, Max = %d,\nNEEoIHU = %0.2f', ...
               sum(Entropy), H, NEEoIHU));    
drawnow;