function  GBPRBM_Visible_Units_Plot_Geometry(opt)
%% Gaussian-Bipolar Restricted Boltzmann Machine
if isstruct(opt)
    V = opt.V;
    H = opt.H;
    W = opt.W;
    if size(W,1)~=V
        error('Matrix W is not of size VxH. Quitting.')
    end
    b_v = opt.b_v; 
    if isfield(opt, 'put_labels')
        put_labels = opt.put_labels;
    else
        put_labels = false;
    end
    if isfield(opt, 'LineStyle')
        LineStyle = opt.LineStyle;
    else
        LineStyle = '-';
    end
    if isfield(opt,'FontSize')==1
        FontSize = opt.FontSize;
    else
        FontSize = 8;
    end
    Hid_Vec = opt.Hid_Vec;
    Centroid = opt.Centroid;   
end
% Set random color in the RGB format
% H hidden units + visible bias, total (H+1) colors
if isfield(opt,'Color')==1
    Color = opt.Color;
else      
    Color = {[0, 0.6, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};
end
LC = length(Color);
if LC<=H+1
    Color = [Color, mat2cell(rand(1,3*(H+1-LC)),1,3*ones(1,H+1-LC))];
end

if V==1
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
            StepY = 1/(H+1)*diff(get(gca, 'YLim'))/2;
            yyP = ones(1,2)*diff(get(gca, 'YLim'))*2/3*(H + 1)/H;
            % Draw a line from "0" to "b_v".
            line([0 b_v],yyP, 'Color', Color{j}, ...
                'Marker', 'sq', 'LineStyle', LineStyle);
            if put_labels
                text(b_v,yyP(1)+0.2*StepY, 'b_v', 'BackgroundColor', 'yellow',...
                     'HorizontalAlignment','center','FontSize', FontSize,...
                     'VerticalAlignment','bottom');        
            end
            yyL = yyP - 1.3*StepY;
            % Vertical line
            line(ones(1,2)*b_v, [yyP(1) yyL(1)], 'Color', Color{j+1}, ... 
                 'Marker', 'sq','LineStyle', LineStyle);
            % Horizontal line
            xx = [Centroid{1}(1,1) Centroid{j}(1,2)];
            yyP = yyL;
            line(xx,yyP, 'Color', 'red', 'Marker', 'sq', 'LineStyle', LineStyle);
            % Right
            text(Centroid{1}(1,1)-1/2*W(1), yyP(1)+0.35*StepY, '+W(1)', ...
                 'BackgroundColor', 'white', 'FontSize', FontSize);
            text(Centroid{1}(1,2)+1/2*W(1), yyP(1)+0.35*StepY, '-W(1)', ...
                 'BackgroundColor', 'white', 'FontSize', FontSize);
        else
            yyL = yyP - 1.5*StepY;
            % Set random color in the RGB format        
            for k=1:2^(j-1)                
                xx = ones(1,2)*Centroid{j-1}(1,k);
                %yy = yyL - 0.5*StepY*rand;
                yc = [yyP(1) yyL(1)];
                % Vertical line
                line(xx,yc, 'Color', Color{j+1}, 'LineStyle', LineStyle);

                % Dimension: v1
                xx = [Centroid{j}(1,2*(k-1)+1) Centroid{j}(1,2*k)];
                line(xx,yyL, 'Color', Color{j+1}, 'Marker', 'sq', 'LineStyle', LineStyle);
                % Right
                str = sprintf('+W(%d)', j);
                text(Centroid{j}(1,2*(k-1) +1)-1/2*W(j), yyL(1)+0.35*StepY, str, ...
                                 'BackgroundColor', 'white', 'FontSize', FontSize, ...
                                 'HorizontalAlignment', 'center');
                str = sprintf('-W(%d)', j);
                text(Centroid{j}(1,2*k)+1/2*W(j), yyL(1)+0.35*StepY, str, ...
                                 'BackgroundColor', 'white', 'FontSize', FontSize, ...
                                 'HorizontalAlignment','center');
                if (j==H) && put_labels
                    x = Centroid{j}(1,2*(k-1)+1);
                    y = yyL(1)-0.5*StepY;
                    str = num2str(Hid_Vec{j}(:,2*(k-1)+1),'%+i');
                    text(x,y, str, 'BackgroundColor', 'yellow', ...
                         'HorizontalAlignment','center', 'FontSize', FontSize, ...
                         'VerticalAlignment','top');

                    x = Centroid{j}(1,2*k);
                    str = num2str(Hid_Vec{j}(:,2*k),'%+i');
                    text(x,y, str, 'BackgroundColor', 'yellow',...
                         'HorizontalAlignment','center', 'FontSize', FontSize,...
                         'VerticalAlignment','top');
               end
            end
            yyP = yyL;
        end
    end
else
    %%                V=2 or V =3
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
            if V == 3
                % Mark visible bias "b_v" point with a dot,
                % i.e. draw a line from "b_v" to "b_v" of zero length.
                line([0 b_v(1)], [0 b_v(2)], [0 b_v(3)], 'LineWidth', 2, ...
                     'Color', Color{j}, 'Marker', 'o', 'LineStyle', LineStyle);                
                zz = [Centroid{1}(3,1) Centroid{j}(3,2)];
                line(xx,yy,zz, 'LineWidth', 2, 'Color', Color{j+1}, ...
                    'Marker', 'o', 'LineStyle', LineStyle);
                if put_labels
                    text(0, 0, 0, 'Origin', 'BackgroundColor', 'yellow', ...
                        'FontSize', FontSize, 'VerticalAlignment', 'top',...
                        'HorizontalAlignment', 'right');
                    text(b_v(1), b_v(2), b_v(3),  'b_v', 'BackgroundColor', 'yellow',...
                        'FontSize', FontSize);
                end
            elseif V == 2
                % Mark visible bias "b_v" point with a dot,
                % i.e. draw a line from "b_v" to "b_v" of zero length.
                line([0 b_v(1)], [0 b_v(2)], 'LineWidth', 2, 'LineStyle', LineStyle, ...
                             'Color', Color{j}, 'Marker', 'o');
                line(xx,yy, 'LineWidth', 2, 'Color', Color{j+1}, 'Marker', 'o', 'LineStyle', LineStyle);
                if put_labels
                    text(0, 0, 'Origin', 'BackgroundColor', 'yellow', ...
                        'FontSize', FontSize, 'VerticalAlignment', 'top', ...
                        'HorizontalAlignment', 'right');
                    text(b_v(1), b_v(2), 'b_v', 'BackgroundColor', 'yellow',...
                        'FontSize', FontSize);
                end
            end        
        else
            for k=1:2^(j-1)
                % Dimension: v1
                xx = [Centroid{j}(1,2*(k-1)+1) Centroid{j}(1,2*k)];
                % Dimension: v2
                yy = [Centroid{j}(2,2*(k-1)+1) Centroid{j}(2,2*k)]; 
                if V == 3
                    zz = [Centroid{j}(3,2*(k-1)+1) Centroid{j}(3,2*k)];
                    line(xx,yy,zz, 'LineWidth', 2, 'Color', Color{j+1}, 'Marker', 'o',...
                        'LineStyle', LineStyle);
                elseif V==2
                    line(xx,yy, 'LineWidth', 2, 'Color', Color{j+1}, 'Marker', 'o', ...
                        'LineStyle', LineStyle); 
                end
                
                if (j == H) && put_labels
                    x = Centroid{j}(1,2*(k-1)+1);
                    y = Centroid{j}(2,2*(k-1)+1);    
                    str = ['[' sprintf('%+i',Hid_Vec{j}(:,2*(k-1)+1)') ']'];
                    if V == 3
                        z = Centroid{j}(3,2*(k-1)+1); 
                        text(x,y,z, str, 'BackgroundColor', 'yellow', 'FontSize', FontSize);                    
                    elseif V==2   
                        text(x,y, str, 'BackgroundColor', 'yellow', 'FontSize', FontSize);  
                    end     
                    
                    x = Centroid{j}(1,2*k);
                    y = Centroid{j}(2,2*k);    
                    str = ['[' sprintf('%+i',Hid_Vec{j}(:,2*k)') ']'];
                    if V == 3
                        z = Centroid{j}(3,2*k); 
                        text(x,y,z, str, 'BackgroundColor', 'yellow', 'FontSize', FontSize);                    
                    elseif V==2   
                        text(x,y, str, 'BackgroundColor', 'yellow', 'FontSize', FontSize);  
                    end 
                end    
            end
        end
    end
end