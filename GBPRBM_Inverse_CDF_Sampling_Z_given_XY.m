function v3 = GBPRBM_Inverse_CDF_Sampling_Z_given_XY(opt)
v1 = opt.v1;
v2 = opt.v2;
ZLimLow  = opt.ZLimLow;
ZLimHigh = opt.ZLimHigh;
N_Samples = opt.N_Samples;
Nv3 = opt.Nv3;
p_h = opt.p_h;
TT = opt.TT;

% Memory preallocation for Z-coordinates of the dots to be scattered
v3 = zeros(1,N_Samples);
% Realizations of uniformly [0,1] distributed random variable
u = rand(1,N_Samples);
% Quantized values of the Z-axis (Nv3 possible values)
z = linspace(ZLimLow, ZLimHigh, Nv3);
Qstep = abs(ZLimLow - ZLimHigh)/(Nv3);
% For every sample "k", do inverse transform sampling for conditional
% distribution F(v3|v2,v1)
for k=1:N_Samples   
    % Realization of the uniformly [0,1] distributed random variable
    Value = u(k);
    % Conditional cumulative distribution function F(v3|v2,v1) 
    [Fv3] = CDF_Z_given_XY(z,v1(k),v2(k),p_h,TT,opt);
    % Find "v3" for which "Fv3" is equal to "Value" but
    % firstly find index of "Fv3" at which Fv3(idx) = Value
    % Instead of "equal" use "greater or equal" due to quantization
    % limitations
    idc = find(Fv3 >= Value);
    if isempty(idc)
        % Then "Value" is very close to 1 which corresponds to the last
        % element of the CDF. Map this value to the random variable axes
        v3(k) = z(end);           
    else
        % Since many values of Fv1 can be greater than "Value" use the 
        % closest one, i.e. the first index
        n = idc(1);
        % Piece-wise linearization
        if n==1
            % Line equation
            v3(k) = Value*Qstep/Fv3(n);            
        else        
            d = Value - Fv3(n-1);
            % Line equation
            v3(k) = z(n-1) + d*Qstep/(Fv3(n) - Fv3(n-1));   
        end
    end
end

function [Fv3] = CDF_Z_given_XY(v3,v1,v2,p_h,TT,opt)
H = opt.H;
W = opt.W;
sigma_v = opt.sigma_v;
b_v = opt.b_v;
mu  = repmat(b_v,1,2^H) + W*TT;
a = p_h.*normpdf(v1, mu(1,:), sigma_v(1)).*normpdf(v2, mu(2,:), sigma_v(2));
Z = sum(a);
Fv3=zeros(size(v3));
for k=1:length(a)    
    Fv3 = Fv3 + a(k)*normcdf(v3,mu(3,k),sigma_v(3));
end
Fv3 = Fv3/Z;