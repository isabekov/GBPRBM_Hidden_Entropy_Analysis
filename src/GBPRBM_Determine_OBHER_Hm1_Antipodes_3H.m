function opt = GBPRBM_Determine_OBHER_Hm1_Antipodes_3H(opt)
% This function finds a one-bit hidden entropy region in the space of hidden
% bias for a given Gaussian-Bipolar Restricted Boltzmann Machine model and
% (H-1) antipode hidden units.
% Basically, it solves: p(h_j,h_a,h_b) = p(h_j,-h_a,-h_b)
% Input variables:
% opt.W       : model weights
% opt.b_v     : visible bias vector
% opt.Sinv    : inverse of the "covariance matrix"
% opt.fxd     : index of the fixed hidden unit
% opt.a       : index of the first antipode hidden unit
% opt.b       : index of the second antipode hidden unit
% opt.b_h_lim

if isstruct(opt)
    unpack_struct(opt);   
end
%% ================== Hidden entropy = 1 bit ===================
fprintf('\n===============================================================\n');
fprintf('One-Bit Hidden Entropy Region, p(h_j,h_a,h_b) = p(h_j,-h_a,-h_b)\n');
fprintf('----------------------------------------------------------------\n');
fprintf('======================= j=%i, h_j=%+i ========================\n', fxd, h_fxd);
fprintf('Hidden entropy = 1 bit if\n');
h_s = sign(W(:,a)'*Sinv*W(:,b));

B_range = [-W(:,b)'*Sinv*(b_v + h_fxd*W(:,fxd) + h_s*W(:,a)), ...
           -W(:,b)'*Sinv*(b_v + h_fxd*W(:,fxd) - h_s*W(:,a))]; 

fprintf('    %f < b_h(%i) < %f.\n', B_range(1), b, B_range(2));

fprintf('and only for\n'); 

h_jab = sign(W(:,fxd)'*Sinv*(h_s*W(:,a) + W(:,b)));
if h_fxd == +1
   %J_lob = max(J_val); 
   J_lob = -W(:,fxd)'*Sinv*(b_v - h_jab*(h_s*W(:,a) + W(:,b)));
   fprintf('    b_h(%i) > %f.\n', fxd, J_lob);
   J_range = [J_lob, b_h_lim{fxd}(2)];   
else
   %J_upb = min(J_val);                 
   J_upb = -W(:,fxd)'*Sinv*(b_v + h_jab*(h_s*W(:,a) + W(:,b)));
   J_range = [b_h_lim{fxd}(1), J_upb]; 
   fprintf('    b_h(%i) < %f.\n', fxd, J_upb);
end

B_intercept = -(b_v + h_fxd*W(:,fxd))'*Sinv*(h_s*W(:,b) + W(:,a));
slope =  -h_s;

% Binary dictionary: char(44-1) = "+" and char(44+1) = "-"
fprintf('    b_h(%i) = %sb_h(%i)%s,\n', a, char(44 - slope), b, sprintf('%+f', B_intercept));


A_srt = slope*B_range(1) + B_intercept;
A_stp = slope*B_range(2) + B_intercept;    

A_range = [A_srt, A_stp];
% ======== Mapping =========
% fxd = 1 => a = 2; b = 3;   
% fxd = 2 => a = 1; b = 3;   
% fxd = 3 => a = 1; b = 2;   

%                   j = fxd
%                 
%                          C point
%     ^	+------------------*                x <----  right corner
%     |	| 	               \-
%     |	|	                 \-
%     |	|	                   \-
%     |	|	                     \-
%     |	|	                       \-  
%     |	|                            \-
%     |	|                              \-
% b_h(j)|	                             \- 
%     |	|                                   * D point
%     |	|                                   |
%     |	|                                   |
%     |	|                                   |
%     |	|                                   |
%     |	|                                   |
%     |	|                                   |
%     |	|                                   |
%     v	+-----------------------------------+
% 	<-------------  b_h(b) ------------>

 if ~(h_s*(1 + h_jab)*W(:,fxd)'*Sinv*W(:,a) > (1 - h_jab)*W(:,fxd)'*Sinv*W(:,b))
      if ~(h_s*(1 + h_jab)*W(:,fxd)'*Sinv*W(:,a) + 2*h_s*W(:,b)'*Sinv*W(:,a) > (1 - h_jab)*W(:,fxd)'*Sinv*W(:,b))
         % Two corners are cut
         if (h_fxd == -1)
            % ======= Point K ======= (left)
            % Bias "B" is set to: b_h(b) = -W(:,b)'*Sinv*(b_v + h_s*W(:,a) - W(:,fxd));
            b_h_fxd_K = -W(:,fxd)'*Sinv*(b_v - h_s*W(:,a) + W(:,b)) + 2*h_s*W(:,a)'*Sinv*W(:,b);
            % ======= Point L ======= (right)
            % Bias "B" is set to: b_h(b) = -W(:,b)'*Sinv*(b_v - h_s*W(:,a) - W(:,fxd));
            b_h_fxd_L = -W(:,fxd)'*Sinv*(b_v - h_s*W(:,a) + W(:,b));
            x = [J_range(1), b_h_fxd_K, b_h_fxd_L, J_range(1)];
            y = [repmat(A_range(1),1,2), repmat(A_range(2),1,2)];
            z = [repmat(B_range(1),1,2), repmat(B_range(2),1,2)];
         elseif (h_fxd == +1)
            % ======= Point M ======= (left)
            % Bias "B" is set to: b_h(b) = -W(:,b)'*Sinv*(b_v + h_s*W(:,a) + W(:,fxd));
            b_h_fxd_M = -W(:,fxd)'*Sinv*(b_v + h_s*W(:,a) - W(:,b));
            % ======= Point N ======= (right)
            % Bias "B" is set to: b_h(b) = -W(:,b)'*Sinv*(b_v - h_s*W(:,a) + W(:,fxd));
            b_h_fxd_N = -W(:,fxd)'*Sinv*(b_v + h_s*W(:,a) - W(:,b)) - 2*h_s*W(:,a)'*Sinv*W(:,b);
            x = [b_h_fxd_M, J_range(2), J_range(2), b_h_fxd_N];
            y = [repmat(A_range(1),1,2), repmat(A_range(2),1,2)];
            z = [repmat(B_range(1),1,2), repmat(B_range(2),1,2)];
          end 
      else
         % One corner is cut
         if (h_fxd == -1)
            % ======= Point C =======
            b_h_b_C = -W(:,b)'*Sinv*(b_v - h_s*W(:,a) - h_jab*W(:,fxd)) + h_s*(1 + h_jab)*W(:,fxd)'*Sinv*W(:,a);
            b_h_a_C= -W(:,a)'*Sinv*(b_v + W(:,b) + h_jab*W(:,fxd)) + h_s*(1-h_jab)*W(:,fxd)'*Sinv*W(:,b);
            % ======= Point D =======
            b_h_fxd_D = -W(:,fxd)'*Sinv*(b_v - h_s*W(:,a) + W(:,b));
            x = [J_range, J_range(2), b_h_fxd_D, J_range(1)];
            y = [repmat(A_range(1),1,2), b_h_a_C, repmat(A_range(2),1,2)];
            z = [repmat(B_range(1),1,2), b_h_b_C, repmat(B_range(2),1,2)];
         elseif (h_fxd == +1)
            % ======= Point E =======
            b_h_b_E = -W(:,b)'*Sinv*(b_v + h_s*W(:,a) + h_jab*W(:,fxd)) - h_s*(1 + h_jab)*W(:,fxd)'*Sinv*W(:,a);
            b_h_a_E = -W(:,a)'*Sinv*(b_v - W(:,b) - h_jab*W(:,fxd)) - h_s*(1-h_jab)*W(:,fxd)'*Sinv*W(:,b);
            % ======= Point F =======
            b_h_fxd_F = -W(:,fxd)'*Sinv*(b_v + h_s*W(:,a) - W(:,b));

            x = [J_range(1), b_h_fxd_F, repmat(J_range(2),1,2), J_range(1)];
            y = [b_h_a_E, repmat(A_range(1),1,2),  repmat(A_range(2),1,2)];
            z = [b_h_b_E, repmat(B_range(1),1,2),  repmat(B_range(2),1,2)];
          end      
      end    
 elseif ~(h_s*(1 - h_jab)*W(:,fxd)'*Sinv*W(:,a) < (1 + h_jab)*W(:,fxd)'*Sinv*W(:,b))
      if ~(h_s*(1 - h_jab)*W(:,fxd)'*Sinv*W(:,a) - 2*h_s*W(:,b)'*Sinv*W(:,a) < (1 + h_jab)*W(:,fxd)'*Sinv*W(:,b))
         % Two corners are cut
         if (h_fxd == -1)
            % ======= Point O ======= (left)
            % Bias "B" is set to: b_h(b) = -W(:,b)'*Sinv*(b_v + h_s*W(:,a) - W(:,fxd));
            b_h_fxd_O = -W(:,fxd)'*Sinv*(b_v + h_s*W(:,a) - W(:,b));
            % ======= Point P ======= (right)
            % Bias "B" is set to: b_h(b) = -W(:,b)'*Sinv*(b_v - h_s*W(:,a) - W(:,fxd));
            b_h_fxd_P = -W(:,fxd)'*Sinv*(b_v + h_s*W(:,a) - W(:,b)) + 2*h_s*W(:,a)'*Sinv*W(:,b);
            x = [J_range(1), b_h_fxd_O, b_h_fxd_P, J_range(1)];
            y = [repmat(A_range(1),1,2), repmat(A_range(2),1,2)];
            z = [repmat(B_range(1),1,2), repmat(B_range(2),1,2)];
         elseif (h_fxd == +1)
            % ======= Point Q ======= (left)
            % Bias "B" is set to: b_h(b) = -W(:,b)'*Sinv*(b_v + h_s*W(:,a) + W(:,fxd));
            b_h_fxd_Q = -W(:,fxd)'*Sinv*(b_v - h_s*W(:,a) + W(:,b)) - 2*h_s*W(:,a)'*Sinv*W(:,b);
            % ======= Point R ======= (right)
            % Bias "B" is set to: b_h(b) = -W(:,b)'*Sinv*(b_v - h_s*W(:,a) + W(:,fxd));
            b_h_fxd_R = -W(:,fxd)'*Sinv*(b_v - h_s*W(:,a) + W(:,b)) ;
            x = [b_h_fxd_Q, J_range(2), J_range(2), b_h_fxd_R];
            y = [repmat(A_range(1),1,2), repmat(A_range(2),1,2)];
            z = [repmat(B_range(1),1,2), repmat(B_range(2),1,2)];
          end 
      else
         % One corner is cut
         if (h_fxd == -1)
            % ======= Point G =======
            b_h_b_G = -W(:,b)'*Sinv*(b_v + h_s*W(:,a) + h_jab*W(:,fxd)) + h_s*(1 - h_jab)*W(:,fxd)'*Sinv*W(:,a);
            b_h_a_G = -W(:,a)'*Sinv*(b_v - W(:,b) - h_jab*W(:,fxd)) + h_s*(1+h_jab)*W(:,fxd)'*Sinv*W(:,b);
            % ======= Point H =======
            b_h_fxd_H = -W(:,fxd)'*Sinv*(b_v + h_s*W(:,a) - W(:,b));

            x = [J_range(1), b_h_fxd_H, repmat(J_range(2),1,2), J_range(1)];
            y = [repmat(A_range(1),1,2), b_h_a_G, repmat(A_range(2),1,2)];
            z = [repmat(B_range(1),1,2), b_h_b_G, repmat(B_range(2),1,2)];
         elseif (h_fxd == +1)
            % ======= Point I =======
            b_h_b_I = -W(:,b)'*Sinv*(b_v - h_s*W(:,a) - h_jab*W(:,fxd)) - h_s*(1 - h_jab)*W(:,fxd)'*Sinv*W(:,a);
            b_h_a_I = -W(:,a)'*Sinv*(b_v + W(:,b) + h_jab*W(:,fxd)) - h_s*(1+h_jab)*W(:,fxd)'*Sinv*W(:,b);
            % ======= Point J =======
            b_h_fxd_J = -W(:,fxd)'*Sinv*(b_v - h_s*W(:,a) + W(:,b));

            x = [J_range, J_range(2), b_h_fxd_J, J_range(1)];
            y = [repmat(A_range(1),1,2), repmat(A_range(2),1,2), b_h_a_I];
            z = [repmat(B_range(1),1,2), repmat(B_range(2),1,2), b_h_b_I];
        end
      end
 else
    % A rectangular region. No corners are cut. 
    x = [J_range fliplr(J_range)];
    y = [repmat(A_range(1),1,2) repmat(A_range(2),1,2)];
    z = [repmat(B_range(1),1,2) repmat(B_range(2),1,2)];     
 end 

C = Color_Map(round(log2(2)/log2(4)*size(Color_Map,1)),:);
switch fxd 
   case 1
    fill3(x, y, z, C);    
   case 2
    fill3(y, x, z, C);   
   case 3     
    fill3(y, z, x, C);        
end