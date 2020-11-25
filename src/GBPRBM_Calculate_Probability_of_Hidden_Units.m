function [p_h TT] = Calculate_Probability_of_Hidden_Units(opt)
H = opt.H;
W = opt.W;
Sinv = opt.Sinv;
b_v = opt.b_v;
b_h = opt.b_h;
Nm = zeros(2^H,1);
Zn  = 0;
% Index terms
r = (2.^fliplr([0:(H-1)]))/2;
s = sum(r) + 1;
TT = (double(dec2bin(0:(2^H)-1)-'1') +  double(dec2bin(0:(2^H)-1)-'0'))';
Sum = -Inf;
% Find the largest element
for h=TT
   a = 1/2*(2*b_v + W*h)'*Sinv*W*h + b_h'*h;
   if a > Sum
      Sum = a; 
   end   
end

for h=TT 
   % Decimal index, e.g. for h=[-1 -1]', idx = 1;  for h=[+1 +1]', idx = 4; 
   idx = r*h + s;
   Nm(idx) = exp(1/2*(2*b_v + W*h)'*Sinv*W*h + b_h'*h - Sum);   
   Zn = Zn + Nm(idx);
end
% Normalization
p_h = Nm'/Zn;