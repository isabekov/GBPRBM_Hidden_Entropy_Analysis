function GBLRBM_Visible_Units_Plot_Geometry(opt)
%% Gaussian-Bernoulli Restricted Boltzmann Machine
if isstruct(opt)
    V = opt.V;
    H = opt.H;
    W = opt.W;
    if size(W,1)~=V
        error('Matrix W is not of size VxH. Quitting.')
    end
    b_v = opt.b_v; 
    Centroid = opt.Centroid;   
end

text(0, 0, 'Origin', 'BackgroundColor', 'yellow', 'FontSize', 8, ...
                   'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
line([0 b_v(1)], [0 b_v(2)], 'Color', 'magenta', 'Marker', 'o');
% Truth table for the hidden neurons
% Superimpose these two matrices
%          Values: "+1" and "0"
TT =  double(dec2bin(0:(2^H)-1)-'0');
% For every configuration of the hidden layer  (2^H different configurations)
for g=TT' % Hidden vector "g" of size: (H x 1)
    v = b_v + W*g;    
    if sum(g)==1
        LineColor = 'black';
    else
        LineColor = 'blue';
    end
    
    line([b_v(1) v(1)], [b_v(2) v(2)], 'Color', LineColor, 'Marker', 'sq', 'LineWidth', 2);
    if sum(g)==0
        str = ['b_v,[' sprintf('%i',g) ']'];
    else
        str = ['[' sprintf('%i',g) ']'];
    end
    text(v(1),v(2), str, 'BackgroundColor', 'yellow', 'FontSize', 8);
    
end

%% Grid with dotted lines
% Draw lines to represent the process of obtaining coordinates for
% the Gaussian pdf, i.e. for a so called "the walk" process.
% Mean of each Gaussian will be called "centroid".
for j=1:H
    if j==1
        % Draw a line between points "b_v + W(:,1)" and "b_v - W(:,1)"
        % Deobfuscation of the indices: 
        % Centroid{j}(1,2) = v1 coordinate (X-axis) of the "-1" hidden
        % neuron configuration of the depth level {j}, because configurations
        % "h_j = -1" (backward) are stored at even indices. 
        % In other words, since j=1, that is an X-coordinate of "b_v - W(:,1)"
        
        % Centroid{1}(2,1) = v2 coordinate (Y-axis) of the "+1" hidden
        % neuron configuration of the depth level {1}, because configurations
        % "h_j = +1" (forward) are stored at odd indices.
        % In other words, that is an Y-coordinate of "b_v + W(:,1)"
        
        xx = [Centroid{1}(1,1) Centroid{j}(1,2)];
        yy = [Centroid{1}(2,1) Centroid{j}(2,2)];
        line(xx,yy, 'Color', 'red', 'LineStyle', ':');
        % Mark visible bias "b_v" point with a dot,
        % i.e. draw a line from "b_v" to "b_v" of zero length.
        %line([0 b_v(1)], [0 b_v(2)], 'Color', 'magenta', 'Marker', 'sq', 'LineStyle', '--');
        %text(b_v(1), b_v(2), 'b_v', 'BackgroundColor', 'yellow', 'FontSize', 8);
    else
        % Set random color in the RGB format
        Color = 'black';
        for k=1:2^(j-1)
            % Dimension: v1
            xx = [Centroid{j}(1,2*(k-1)+1) Centroid{j}(1,2*k)];
            % Dimension: v2
            yy = [Centroid{j}(2,2*(k-1)+1) Centroid{j}(2,2*k)];
            line(xx,yy, 'Color', Color, 'LineStyle', ':', 'LineWidth', 2);
        end
    end
end
