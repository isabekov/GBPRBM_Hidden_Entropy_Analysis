function v2 = GBPRBM_Inverse_CDF_Sampling_Y_given_X(opt)
v1 = opt.v1;
YLimLow  = opt.YLimLow;
YLimHigh = opt.YLimHigh;
N_Samples = opt.N_Samples;
Nv2 = opt.Nv2;
p_h = opt.p_h;
TT = opt.TT;

% Memory preallocation for Y-coordinates of the dots to be scattered
v2 = zeros(1,N_Samples);
% Realizations of uniformly [0,1] distributed random variable
u = rand(1,N_Samples);
% Quantized values of the Y-axis (Nv2 possible values)
y = linspace(YLimLow, YLimHigh, Nv2);
Qstep = abs(YLimLow - YLimHigh)/(Nv2);
% For every sample "k", do inverse transform sampling for conditional
% distribution F(v2|v1)
for k=1:N_Samples  
    % Realization of the uniformly [0,1] distributed random variable
    Value = u(k);
    % Conditional cumulative distribution function F(v2|v1) 
    [Fv2] = CDF_Y_given_X(y,v1(k),p_h,TT,opt);
    % Find "v2" for which "Fv2" is equal to "Value" but
    % firstly find index of "Fv2" at which Fv2(idx) = Value
    % Instead of "equal" use "greater or equal" due to quantization
    % limitations
    idc = find(Fv2 >= Value);
    if isempty(idc)
        % Then "Value" is very close to 1 which corresponds to the last
        % element of the CDF. Map this value to the random variable axes
        v2(k) = y(end); 
    else
        % Since many values of Fv1 can be greater than "Value" use the 
        % closest one, i.e. the first index
        n = idc(1);
        % Piece-wise linearization
        if n==1            
            % Line equation
            v2(k) = Value*Qstep/Fv2(n);             
        else
            d = Value - Fv2(n-1);
            % Line equation
            v2(k) = y(n-1) + d*Qstep/(Fv2(n) - Fv2(n-1)); 
        end
    end
end

function [Fv2] = CDF_Y_given_X(v2,v1,p_h,TT,opt)
H = opt.H;
W = opt.W;
sigma_v = opt.sigma_v;
b_v = opt.b_v;
mu  = repmat(b_v,1,2^H) + W*TT;
a = p_h.*normpdf(v1, mu(1,:), sigma_v(1));
Z = sum(a);
Fv2=zeros(size(v2));
for k=1:length(a)    
    Fv2 = Fv2 + a(k)*normcdf(v2,mu(2,k),sigma_v(2));
end
Fv2 = Fv2/Z;