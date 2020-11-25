function v1 = GBPRBM_Inverse_CDF_Sampling_X(opt)
XLimLow  = opt.XLimLow;
XLimHigh = opt.XLimHigh;
N_Samples = opt.N_Samples;
Nv1 = opt.Nv1;
p_h = opt.p_h;
TT = opt.TT;
% Memory preallocation for X-coordinates of the dots to be scattered
v1 = zeros(1,N_Samples);
% Realizations of uniformly [0,1] distributed random variable
u = rand(1,N_Samples);
% Quantized values of the X-axis (Nv1 possible values)
x = linspace(XLimLow, XLimHigh, Nv1);
Qstep = abs(XLimLow - XLimHigh)/(Nv1);
% Cumulative distribution function of F(v1)
Fv1 = CDF_X(x,p_h,TT,opt);
% For every sample "k", do inverse transform sampling for F(v1)
for k=1:N_Samples
    % Realization of the uniformly [0,1] distributed random variable
    Value = u(k);
    % Find "v1" for which "Fv1" is equal to "Value" but
    % firstly find index of "Fv1" at which Fv1(idx) = Value.
    % Instead of "equal" use "greater or equal" due to quantization
    % limitations
    idc = find(Fv1 >= Value);    
    if isempty(idc)
        % Then "Value" is very close to 1 which corresponds to the last
        % element of the CDF. Map this value to the random variable axes
        v1(k) = x(end); 
    else
        % Since many values of Fv1 can be greater than "Value" use the 
        % closest one, i.e. the first index
        n = idc(1);
        % Piece-wise linearization
        if n==1
            % Line equation
            v1(k) = Value*Qstep/Fv1(n);            
        else
            d = Value - Fv1(n-1);
            % Line equation
            v1(k) = x(n-1) + d*Qstep/(Fv1(n) - Fv1(n-1)); 
        end
    end   
end


function [Fv1] = CDF_X(v1,p_h,TT,opt)
Fv1=zeros(size(v1));
for k=1:length(p_h)
    % Same operation as if normcdf were used
    %x = (v - (opt.b_v(1) + opt.W(1,:)*TT(:,k)))/(sqrt(2)*opt.sigma_v(1));
    %Fv1 = Fv1 + p_h(k)*(1/2*(1 + erf(x)));
    mu = opt.b_v(1) + opt.W(1,:)*TT(:,k);    
    Fv1 = Fv1 + p_h(k)*normcdf(v1,mu,opt.sigma_v(1));
end