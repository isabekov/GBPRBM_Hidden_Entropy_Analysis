function GBPRBM_Plot_Decision_Region_for_Visible_Units(opt)
V = opt.V;
H = opt.H;
W = opt.W;
Sinv = opt.Sinv;
b_h = opt.b_h;
XLimLow  = opt.v(1).AxesLimMin;
XLimHigh = opt.v(1).AxesLimMax;

if isfield(opt,'LabelPosition')==1
    LabelPosition = opt.LabelPosition;
else      
    LabelPosition = [0.12, 0.1, 0.8];
end

if isfield(opt,'LineWidth') == 1
    LineWidth = opt.LineWidth;
else      
    LineWidth = 1.5;
end

if H > 3
    LabelPosition = [LabelPosition rand(1,H-3)];
end

if isfield(opt,'FontSize') == 1
    FontSize = opt.FontSize;
else
    FontSize = 10;
end
% Set random color in the RGB format
% H hidden units + visible bias, total (H+1) colors
if isfield(opt,'Color') == 1
    Color = opt.Color;
else      
    Color = {[0, 0.6, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};
end
LC = length(Color);
if LC<=H+1
    Color = [Color, mat2cell(rand(1,3*(H+1-LC)),1,3*ones(1,H+1-LC))];
end

for j=1:H    
    P = Sinv*W(:,j);
    switch V 
        case 1
            X = -b_h(j)/P; 
            YLim = get(gca, 'YLim');
            Y = LabelPosition(j)*diff(YLim);
            line([X,X], YLim, 'Color', Color{j+1}, ...
                  'LineStyle', '--');  
            PM = sign(P);
            text(X,Y, sprintf(' h_%d=%+i ', j,  PM), 'HorizontalAlignment','left',... 
                 'BackgroundColor',Color{j+1}, 'Margin', 1, 'EdgeColor', 'black', ...
                 'FontSize', FontSize, 'FontWeight', 'demi');

            text(X,Y, sprintf(' h_%d=%+i ', j, -PM), 'HorizontalAlignment','right',... 
                 'BackgroundColor',Color{j+1}, 'Margin', 1, 'EdgeColor', 'black', ...
                 'FontSize', FontSize, 'FontWeight', 'demi');            
        case 2
            YLimLow  = opt.v(2).AxesLimMin;
            YLimHigh = opt.v(2).AxesLimMax;
            % Make sure that line does not exceed axes boundaries, i.e.
            % the line is in the box with diagonals (XLimLow, XLimHigh) and 
            % (YLimLow, YLimHigh).
            % Setting the values below will also work as long as axes limits
            % are fixed:
            % X_LB = XLimLow;  % also works, but not advised
            % X_UB = XLimHigh; % also works, but not advised   
            XLowByY  = (sign(P(2)/P(1))==+1)*(-P(2)*YLimHigh - b_h(j))/P(1) + ...
                       (sign(P(2)/P(1))==-1)*(-P(2)*YLimLow  - b_h(j))/P(1);
            X_LB = (XLimLow > XLowByY)*XLimLow + (XLimLow < XLowByY)*XLowByY;

            XHighByY = (sign(P(2)/P(1))==+1)*(-P(2)*YLimLow  - b_h(j))/P(1) + ...
                       (sign(P(2)/P(1))==-1)*(-P(2)*YLimHigh - b_h(j))/P(1);
            X_UB = (XLimHigh > XHighByY)*XHighByY + (XLimHigh < XHighByY)*XLimHigh;           

            Y_LB = -(b_h(j) + X_LB*P(1))/P(2);                 
            Y_UB = -(b_h(j) + X_UB*P(1))/P(2);                 
            line([X_LB, X_UB], [Y_LB, Y_UB], 'Color', Color{j+1}, ...
                  'LineStyle', '--', 'LineWidth', LineWidth);  
            x = X_LB + LabelPosition(j)*diff([X_LB, X_UB]);
            y = -(b_h(j) + x*P(1))/P(2);
            if -P(1)/P(2) > 0 
                Rotation = atan(-P(1)/P(2))/(2*pi)*360 - 90;
                PM = sign([cos(Rotation/360*2*pi) sin(Rotation/360*2*pi)]*P);
                text(x,y, sprintf(' h_%d=%+i ', j, PM), 'HorizontalAlignment','left',... 
                     'BackgroundColor',Color{j+1}, 'Margin', 1, 'EdgeColor', 'black', ...
                     'FontSize', FontSize, 'FontWeight', 'demi', 'Rotation', Rotation);
                
                    text(x,y, sprintf(' h_%d=%+i ', j, -PM), 'HorizontalAlignment','right',... 
                     'BackgroundColor',Color{j+1}, 'Margin', 1, 'EdgeColor', 'black',  ...
                     'FontSize', FontSize, 'FontWeight', 'demi', 'Rotation', Rotation);
            else
                Rotation = atan(-P(1)/P(2))/(2*pi)*360 + 90;
                PM = sign([cos(Rotation/360*2*pi) sin(Rotation/360*2*pi)]*P);
                text(x,y, sprintf(' h_%d=%+i ', j, PM), 'HorizontalAlignment','left',... 
                     'BackgroundColor',Color{j+1}, 'Margin', 1, 'EdgeColor', 'black',  ...
                     'FontSize', FontSize, 'FontWeight', 'demi', 'Rotation', Rotation);                

                text(x,y, sprintf(' h_%d=%+i ', j, -PM), 'HorizontalAlignment','right',... 
                     'BackgroundColor',Color{j+1}, 'Margin', 1, 'EdgeColor', 'black',  ...
                     'FontSize', FontSize, 'FontWeight', 'demi', 'Rotation', Rotation);                  
            end
        case 3
            %% Find a plane passing through a rectangle
            % The plane cutting the parallelepiped yields either
            % a convex polygon with 3-6 vertices.
            % This algorithm resembles the "Marching-Cubes" algorithm
            YLimLow  = opt.v(2).AxesLimMin;
            YLimHigh = opt.v(2).AxesLimMax;                
            ZLimLow  = opt.v(3).AxesLimMin;
            ZLimHigh = opt.v(3).AxesLimMax;    
            %% A dictionary "Limits" for lower and upper limits of every axis
            Key = {'xL', 'xU', 'yL', 'yU', 'zL', 'zU'};
            Val = [XLimLow, XLimHigh, YLimLow, YLimHigh, ZLimLow, ZLimHigh];
            Limits = containers.Map(Key,Val);
            
            %% A dictionary "Edges"
            Edges_Keys = cell(1,V*2^(V-1));
            % To create fields that contain cell arrays, place the cell arrays within a value cell array
            tmp = struct('isused', false, 'CoorValue', NaN, 'Point', zeros(V,1), 'Planes', {cell(1,2)});
            Val = repmat({tmp}, 1, length(Edges_Keys));
            k = 1;
            for d = 'xyz'
                for a='LU'
                    for b='LU'
                        switch d
                            case 'x'
                                Edges_Keys{k} = [d, '._', num2str(a), num2str(b)];
                            case 'y'
                                Edges_Keys{k} = [d, '.', num2str(a), '_', num2str(b)];
                            case 'z'
                                Edges_Keys{k} = [d, '.', num2str(a), num2str(b), '_'];
                        end
                        k = k + 1;
                    end
                end
            end
            Edges = containers.Map(Edges_Keys,Val);

            %% Finding edges which intersect with the plane passing through the parallelepiped
            % Make sure that line does not exceed axes' boundaries, i.e.
            % the line is in the box with diagonals (XLimLow, XLimHigh), 
            % (YLimLow, YLimHigh) and (ZLimLow, ZLimHigh).
            % In order to accomplish this, find intersections of the plane with
            % axes' boundaries. A plane will be plotted as a polygon of 
            % MATLAB type "patch" specified by 3-6 intersection points.
            CH = cell(0);
            % Fix "v3" dimension (Z-axis) to ZLimLow and ZLimHigh first.    
            % Plane equation becomes line equation: 
            % x*P(1) + y*P(2) + (ZLimLow*P(3)  + b_h(j)) = 0, or
            % x*P(1) + y*P(2) + (ZLimHigh*P(3) + b_h(j)) = 0
            for z_lim = {'zL', 'zU'}
                z = Limits(z_lim{1});
                % Find intersection of the line with four axes' boundaries
                % which form a rectangle:
                % x is fixed to XLimLow,   y goes from YLimLow to YLimHigh
                % x is fixed to XLimHigh,  y goes from YLimLow to YLimHigh
                % y is fixed to YLimLow,   y goes from XLimLow to XLimHigh
                % y is fixed to YLimHigh,  y goes from XLimLow to XLimHigh
                % If the line doesn't cross axis boundary, then do not include this
                % point 
                for x_lim = {'xL', 'xU'}
                   key = sprintf('y.%s_%s', x_lim{1}(2), z_lim{1}(2));
                   % Only one level of indexing is supported by a containers.Map
                   Rib = Edges(key);
                   Rib.Planes{1} = [z_lim{1} '.xy'];
                   Rib.Planes{2} = [x_lim{1} '.yz'];
                   x = Limits(x_lim{1});
                   y = -(x*P(1) + (z*P(3) + b_h(j)))/P(2);
                   if (YLimLow<=y) && (y<=YLimHigh)
                      % Add this point to the convex hull set
                      Rib.isused = true;
                      Rib.CoorValue = y;
                      Rib.Point = [x y z]';
                      CH = [CH key];
                   end
                   % Remapping
                   Edges(key) = Rib;                   
                end
                for y_lim = {'yL', 'yU'}
                   key = sprintf('x._%s%s', y_lim{1}(2), z_lim{1}(2));
                   % Only one level of indexing is supported by a containers.Map
                   Rib = Edges(key);
                   Rib.Planes{1} = [z_lim{1} '.xy'];
                   Rib.Planes{2} = [y_lim{1} '.xz'];
                   y = Limits(y_lim{1});
                   x = -(y*P(2) + (z*P(3) + b_h(j)))/P(1);
                   if (XLimLow<=x) && (x<=XLimHigh)
                      % Add this point to the convex hull set
                      Rib.isused = true;
                      Rib.CoorValue = x;
                      Rib.Point = [x y z]';
                      CH = [CH key];              
                   end
                   % Remapping
                   Edges(key) = Rib;  
                end        
            end
            % Fix "v1" dimension (X-axis) to XLimLow and XLimHigh first.    
            % Plane equation becomes line equation: 
            % y*P(2) + z*P(3) + (XLimLow*P(1) + b_h(j)) = 0, or
            % y*P(2) + z*P(3) + (XLimHigh*P(1) + b_h(j)) = 0
            for x_lim = {'xL', 'xU'}
                x = Limits(x_lim{1});
                for y_lim = {'yL', 'yU'}
                   key = sprintf('z.%s%s_', x_lim{1}(2), y_lim{1}(2));
                   % Only one level of indexing is supported by a containers.Map
                   Rib = Edges(key);   
                   Rib.Planes{1} = [x_lim{1} '.yz'];
                   Rib.Planes{2} = [y_lim{1} '.xz'];                   
                   y = Limits(y_lim{1});
                   z = -(y*P(2) + (x*P(1) + b_h(j)))/P(3);
                   if (ZLimLow<=z) && (z<=ZLimHigh)
                      % Add this point to the convex hull set                      
                      Rib.isused = true;
                      Rib.CoorValue = z;
                      Rib.Point = [x y z]';
                      CH = [CH key];
                   end
                   % Remapping
                   Edges(key) = Rib;                      
                end              
            end
            
            %% A dictionary "Position" 
            % which return position of a coordinate label (one letter)
            % for the edge format
            Keys = {'x', 'y', 'z'};
            Values = [1,2,3];
            Position = containers.Map(Keys,Values);
            
            %% A dictionary "Plane_Edges" 
            % Keys: planes 
            % Values: four edges which form the plane (for every plane)
            % A parallelepiped has 6 planes
            Plane_Keys = {'xL.yz', 'xU.yz', 'yL.xz', 'yU.xz', 'zL.xy', 'zU.xy'};
            Plane_Values = cell(1,length(Plane_Keys));
            for k=1:length(Plane_Keys)
                Plane_Borders = cell(1,4);
                str = blanks(3);
                % Constant term
                idx = Position(Plane_Keys{k}(1));
                Constant_Term = Plane_Keys{k}(2);                
                str(idx) = Constant_Term;
                cnt = 1;
                Borders = Plane_Keys{k}(4:5);
                for Pr=1:2
                    Parent = Borders(Pr); 
                    for a='LU'                                     
                       str(Position(Parent)) = '_';
                       idx = Position(Borders((Pr==1)*2 + (Pr==2)*1));
                       str(idx) = a;
                       Plane_Borders{cnt} = [Parent '.' str];
                       cnt = cnt + 1;
                    end
                end
                Plane_Values{k} = Plane_Borders;
            end
            Plane_Edges = containers.Map(Plane_Keys,Plane_Values);
            %% Initial edge and plane
            Rib = CH{1};
            Plane_Start = Edges(Rib).Planes{1};
            CH = setdiff(CH, Rib, 'legacy');
            Convex_Hull = Edges(Rib).Point;
            while(~isempty(CH))      
                % Consequent evaluation 
                % Idea: follow the line on every plane and get
                % edges which intersect the line
                [Rib, Plane_Start] = Get_Edge_Neighbors(Rib, Plane_Start, CH, Edges, Plane_Edges);
                % Remove edge from the CH list
                CH = setdiff(CH, Rib, 'legacy');
                % Convex hull X,Y,Z coordinates
                Convex_Hull = [Convex_Hull Edges(Rib).Point];
            end
            % Draw a filled polygon         
            patch(Convex_Hull(1,:),Convex_Hull(2,:),Convex_Hull(3,:), Color{j+1});
            %% Alternative solution:
            % (not advised, because plane boundaries will exceed axes boundaries)
            %     x = linspace(XLimLow, XLimHigh, 200);
            %     y = linspace(YLimLow, YLimHigh, 200);
            %     [X,Y] = meshgrid(x,y);    
            %     Z=(-b_h(j) - P(1)*X - P(2)*Y)/P(3);
            %     surf(X,Y,Z);
    end
end       

function [Edge, Plane] = Get_Edge_Neighbors(Rib,Plane_Start, CH,  Edges, Plane_Edges)
% Every edge is shared by two planes. Remove the starting edge
Plane = setdiff(Edges(Rib).Planes, Plane_Start, 'legacy');
% Edges of the analyzed plane (remove starting edge)
Analyzed_Edges = setdiff(Plane_Edges(Plane{1}),Rib, 'legacy');
% Among three other edges find the one which intersects with the line
Edge = intersect(Analyzed_Edges, CH);
% Cell to string
Edge = Edge{1};