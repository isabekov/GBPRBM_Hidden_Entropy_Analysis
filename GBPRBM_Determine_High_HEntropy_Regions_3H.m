function opt = GBPRBM_Determine_High_HEntropy_Regions_3H(opt)
% Input variables:
% opt.W
% opt.b_v
% opt.Sinv
% opt.h_const
% opt.fxd
% opt.a
% opt.b
% opt.b_h_lim

if isstruct(opt)
    unpack_struct(opt);   
end
%% ================== Hidden entropy = 1.585 bits (A is fixed) ===================
fprintf('\n=========================================================\n');
fprintf('Hidden entropy vs. hidden bias (H=3), with one bias fixed:\n')
fprintf('b_h(%i) = %f\n', fxd, b_h_fxd);

fprintf('--------------------------------\n');
fprintf('Hidden entropy = log2(3) bits if\n');
fprintf('h_const(b)=%i, b=%i\n', h_const(b), b);

h_fxd = -sign(4*h_const(a)*W(:,a)'*Sinv*W(:,fxd));
fprintf('    h_j=%+i corresponding to\n',h_fxd);
b_h_a = -W(:,a)'*Sinv*(b_v + h_fxd*W(:,fxd) + h_const(b)*W(:,b));

%fprintf('b_h_a = -W(:,a)''*Sinv*(b_v + h_fxd*W(:,fxd) + h_const(b)*W(:,b))');
fprintf('    b_h(%i) = %f,\n', a, b_h_a);
fprintf('and only for\n'); 
B_val = [-W(:,b)'*Sinv*(b_v + h_fxd*W(:,fxd) + h_const(a)*W(:,a)), ...
        -W(:,b)'*Sinv*(b_v - h_fxd*W(:,fxd) + h_const(a)*W(:,a)), ...
        -W(:,b)'*Sinv*(b_v + h_fxd*W(:,fxd) - h_const(a)*W(:,a)), ...
        -W(:,b)'*Sinv*(b_v - h_fxd*W(:,fxd) - h_const(a)*W(:,a)) + 2*h_const(a)*h_fxd/h_const(b)*W(:,a)'*Sinv*W(:,fxd)]; 
if h_const(b) == +1
   B_lob = max(B_val);             
   fprintf('    b_h(%i) > %f.\n', b, B_lob);
   B_range = [B_lob, b_h_lim{b}(2)];   
else
   B_upb = min(B_val);                          
   B_range = [b_h_lim{b}(1), B_upb]; 
   fprintf('    b_h(%i) < %f.\n', b, B_upb);
end

%% ================== Hidden entropy = 1.585 bits (B is fixed) ===================
fprintf('--------------------------------\n');
fprintf('Hidden entropy = log2(3) bits if\n');
fprintf('h_const(a)=%+i, a=%i\n', h_const(a),a);

h_fxd = -sign(4*h_const(b)*W(:,b)'*Sinv*W(:,fxd));
fprintf('    h_j=%+i corresponding to\n',h_fxd);
b_h_b = -W(:,b)'*Sinv*(b_v + h_fxd*W(:,fxd) + h_const(a)*W(:,a));

fprintf('    b_h(%i) = %f,\n', b, b_h_b);
fprintf('and only for\n'); 
A_val = [-W(:,a)'*Sinv*(b_v + h_fxd*W(:,fxd) + h_const(b)*W(:,b)), ...
        -W(:,a)'*Sinv*(b_v - h_fxd*W(:,fxd) + h_const(b)*W(:,b)), ...
        -W(:,a)'*Sinv*(b_v + h_fxd*W(:,fxd) - h_const(b)*W(:,b)), ...
        -W(:,a)'*Sinv*(b_v - h_fxd*W(:,fxd) - h_const(b)*W(:,b)) + 2*h_const(b)*h_fxd/h_const(a)*W(:,b)'*Sinv*W(:,fxd)];  
if h_const(a) == +1
   A_lob = max(A_val); 
   fprintf('    b_h(%i) > %f.\n', a, A_lob);
   A_range = [A_lob, b_h_lim{a}(2)];   
else
   A_upb = min(A_val);                 
   A_range = [b_h_lim{a}(1), A_upb]; 
   fprintf('    b_h(%i) < %f.\n', a, A_upb);
end
%% ============= Hidden entropy = 1.585 bits (A vs B line equation) ==============
fprintf('------------------------------\n');
fprintf('Hidden entropy = log2(3) bits if\n');
h_fxd = -sign(4*W(:,fxd)'*Sinv*(h_const(a)*W(:,a) + h_const(b)*W(:,b)));
fprintf('    h_j=%+i corresponding to\n',h_fxd);
B_intercept = -(b_v + h_fxd*W(:,fxd))'*Sinv*(h_const(b)/h_const(a)*W(:,b) + W(:,a));
slope =  -h_const(b)/h_const(a);
% Binary dictionary: char(44-1) = "+" and char(44+1) = "-"
fprintf('    b_h(%i) = %sb_h(%i)%s,\n', a, char(44 - slope), b, sprintf('%+f', B_intercept));
fprintf('and only for\n'); 
B_pos = [-W(:,b)'*Sinv*(b_v + h_fxd*W(:,fxd) - h_const(a)*W(:,a)), ...
         -W(:,b)'*Sinv*(b_v + h_fxd*W(:,fxd) - h_const(a)*W(:,a)) - 2*h_const(a)*h_fxd/h_const(b)*W(:,a)'*Sinv*W(:,fxd)];
      
B_neg = [-W(:,b)'*Sinv*(b_v + h_fxd*W(:,fxd) + h_const(a)*W(:,a)), ...
         -W(:,b)'*Sinv*(b_v - h_fxd*W(:,fxd) + h_const(a)*W(:,a))];  
if h_const(b) == +1
   Bd_upb = min(B_pos);                    
   fprintf('    b_h(%i) < %f,\n', b, Bd_upb);
   Bd_lob = max(B_neg); 
   fprintf('    b_h(%i) > %f.\n', b, Bd_lob);   
else
   Bd_lob = max(B_pos); 
   fprintf('    b_h(%i) > %f,\n', b, Bd_lob);   
   Bd_upb = min(B_neg);                 
   fprintf('    b_h(%i) < %f.\n', b, Bd_upb);
end
Bd_range = [Bd_lob, Bd_upb];  
if diff(Bd_range) > 0
    % Equation
    Ad_srt = slope*Bd_lob + B_intercept;
    Ad_stp = slope*Bd_upb + B_intercept;    
    % Attention! The range is not 
    Ad_range = [Ad_srt, Ad_stp];
else
    fprintf('There is no solution! The range of b_h(%i) is invalid!\n', b);
end

%% =================== Hidden entropy = 1 bit ====================================
fprintf('-------------------------\n');
fprintf('Hidden entropy = 1 bit if\n');
% Binary dictionary: char(61-1) = "<" and char(61+1) = ">"
fprintf('    b_h(%i) %s %f,\n', a, char(61 + h_const(a)), b_h_a);
% Binary dictionary: char(61-1) = "<" and char(61+1) = ">"
fprintf('    b_h(%i) %s %f,\n', b, char(61 + h_const(b)), b_h_b);
tmp = -(h_const(a)*W(:,a) + h_const(b)*W(:,b))'*Sinv*(b_v + h_fxd*W(:,fxd));
fprintf('    %sb_h(%i)%sb_h(%i) > %+f.\n', char(44 - h_const(a)), a,...
        char(44 -h_const(b)), b, tmp);       
    
%% =================== Hidden entropy = 2 bits =========================    
fprintf('--------------------------\n');
fprintf('Hidden entropy = 2 bits if\n'); 
% Check the diagonal log2(3)-Entropy range
if diff(Bd_range) > 0
    % Mark two dots with 2-bit entropy
    for k=1:2
        fprintf('    b_h(%i) = %f; ', a, Ad_range(k));
        % ASCII codes: 44: comma, 46: dot.
        fprintf('b_h(%i) = %f%s\n', b, Bd_range(k), char(42+2*k)); 
    end             
else
    fprintf('    b_h(%i) = %f; ', a, b_h_a);
    fprintf('b_h(%i) = %f;\n', b, b_h_b);              
end
% Output argument
opt.A_range = A_range;
opt.B_range = B_range;
opt.b_h_a = b_h_a;
opt.b_h_b = b_h_b;
if exist('Ad_range', 'var') == 1
   opt.Ad_range = Ad_range;
end
if exist('Bd_range', 'var') == 1
   opt.Bd_range = Bd_range;
end
