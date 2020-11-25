function GBPRBM_Plot_High_HEntropy_Regions_3H(opt)
% Inputs:
% opt.A_range = A_range;
% opt.B_range = B_range;
% opt.b_h_a = b_h_a;
% opt.b_h_b = b_h_b;
% Optional inputs:
% opt.Ad_range = Ad_range;
% opt.Bd_range = Bd_range;

if isstruct(opt)
    unpack_struct(opt);   
end
%% ================== Hidden entropy = 1.585 bits (A is fixed) ===================
C = Color_Map(round(log2(3)/log2(4)*size(Color_Map,1)),:);
switch fxd
       case 1
           h = line(b_h_fxd*ones(1,2), b_h_a*ones(1,2), B_range);
       case 2
           h = line(b_h_a*ones(1,2), b_h_fxd*ones(1,2), B_range);
       case 3
           h = line(b_h_a*ones(1,2), B_range, b_h_fxd*ones(1,2));
end
set(h, 'Color', C, 'LineWidth', 3);
hold on;


%% ================== Hidden entropy = 1.585 bits (B is fixed) ===================
C = Color_Map(round(log2(3)/log2(4)*size(Color_Map,1)),:);
switch fxd
       case 1
           h = line(b_h_fxd*ones(1,2), A_range,  b_h_b*ones(1,2));
       case 2
           h = line(A_range,  b_h_fxd*ones(1,2), b_h_b*ones(1,2));
       case 3
           h = line(A_range,  b_h_b*ones(1,2), b_h_fxd*ones(1,2));
end
set(h, 'Color', C, 'LineWidth', 3);
hold on;

%% ============= Hidden entropy = 1.585 bits (A vs B line equation) ==============
if (exist('Ad_range', 'var') == 1) && (exist('Ad_range', 'var') == 1)
    if diff(Bd_range) > 0
        C = Color_Map(round(log2(3)/log2(4)*size(Color_Map,1)),:);
        switch fxd
               case 1
                   h = line(b_h_fxd*ones(1,2), Ad_range, Bd_range);
               case 2
                   h = line(Ad_range,  b_h_fxd*ones(1,2), Bd_range);
               case 3
                   h = line(Ad_range,  Bd_range, b_h_fxd*ones(1,2));
        end
        set(h, 'Color', C, 'LineWidth', 3);
        hold on;       
    end
end

%% =================== Hidden entropy = 1 bit ====================================
% Check the diagonal log2(3)-Entropy range
if diff(Bd_range) > 0
    switch fxd
       case 1
           % In the convex hull there are five points
           x = b_h_fxd*ones(1,5);
           y = [b_h_a*ones(1,2), A_range, (h_const(a)==-1)*A_range(1) + (h_const(a)==+1)*A_range(2)];
           z = [B_range, b_h_b*ones(1,2), (h_const(b)==-1)*B_range(1) + (h_const(b)==+1)*B_range(2)];           
           k = convhull(y,z);
       case 2
           x = [b_h_a*ones(1,2), A_range, (h_const(a)==-1)*A_range(1) + (h_const(a)==+1)*A_range(2)];
           y = b_h_fxd*ones(1,length(x));
           z = [B_range, b_h_b*ones(1,2), (h_const(b)==-1)*B_range(1) + (h_const(b)==+1)*B_range(2)];           
           k = convhull(x,z);
       case 3
           x = [b_h_a*ones(1,2), A_range, (h_const(a)==-1)*A_range(1) + (h_const(a)==+1)*A_range(2)];
           y = [B_range, b_h_b*ones(1,2), (h_const(b)==-1)*B_range(1) + (h_const(b)==+1)*B_range(2)];
           z = b_h_fxd*ones(1,length(x));
           k = convhull(x,y);
    end   
else
    switch fxd
       case 1
           % In the convex hull there are five points
           x = b_h_fxd*ones(1,4);
           y = [A_range(1)*ones(1,2) A_range(2)*ones(1,2)];
           z = repmat(B_range,1,2);  
           k = convhull(y,z);           
       case 2
           x = [A_range(1)*ones(1,2) A_range(2)*ones(1,2)];
           y = b_h_fxd*ones(1,length(x));
           z = repmat(B_range,1,2); 
           k = convhull(x,z);           
       case 3
           x = [A_range(1)*ones(1,2) A_range(2)*ones(1,2)];
           y = repmat(B_range,1,2);
           z = b_h_fxd*ones(1,length(x));
           k = convhull(x,y);           
    end    
end
C = Color_Map(round(log2(2)/log2(4)*size(Color_Map,1)),:);
h = fill3(x(k), y(k), z(k), C);

%% =================== Hidden entropy = 2 bits =========================  
C = Color_Map(round(log2(4)/log2(4)*size(Color_Map,1)),:);
% Check the diagonal log2(3)-Entropy range
if diff(Bd_range) > 0
    % Mark two dots with 2-bit entropy
    for k=1:2
        switch fxd
           case 1
               h = line(b_h_fxd*ones(1,2), Ad_range(k)*[1 1], Bd_range(k)*[1 1]);                    
           case 2
               h = line(Ad_range(k)*[1 1], b_h_fxd*ones(1,2), Bd_range(k)*[1 1]);  
           case 3
               h = line(Ad_range(k)*[1 1], Bd_range(k)*[1 1], b_h_fxd*ones(1,2));  
        end     
        set(h, 'LineWidth', 2, 'Marker', 'sq', ...
               'MarkerEdgeColor',C, 'MarkerFaceColor',C);  
        % uistack(h, 'top');
    end
else
    switch fxd
       case 1
           h = line(b_h_fxd*ones(1,2), b_h_a*[1 1], b_h_b*[1 1]);     
       case 2
           h = line(b_h_a*[1 1], b_h_fxd*ones(1,2),  b_h_b*[1 1]);   
       case 3
           h = line(b_h_a*[1 1], b_h_b*[1 1], b_h_fxd*ones(1,2));   
    end    
   set(h, 'LineWidth', 2, 'Marker', 'sq', ...
          'MarkerEdgeColor',C, 'MarkerFaceColor',C); 
   % uistack(h, 'top');
end
 