function GBPRBM_Most_Distant_Configurations_3V
% Number of visible units
opt.V = 3;
% Number of hidden units
opt.H = 3;
% Model geometry
opt.W = [9 2  3;
        2 -6 4; 
        -3 4 5];
% Visible bias
opt.b_v = [8 12 10]';
% Standard deviations for visible units
opt.sigma_v = [2 3 4]';
% Covariance matrix
opt.S = diag(opt.sigma_v.^2);
% Inverse of the "covariance matrix"
opt.Sinv = diag(1./(opt.sigma_v.^2));
% Span ratio for all axes
opt.VRatio = [0.4 0.4 0.4];
% Number of dots to be scattered
opt.N_Samples = 2000;
% Number of samples used to quantize X-axis
opt.Nv1 = 190;
% Number of samples used to quantize Y-axis
opt.Nv2 = 180;
% Number of samples used to quantize Y-axis
opt.Nv3 = 170;
opt.View_Elevation = 60;
opt.View_Azimuth = 233;
opt.Color = {'red', 'green', 'blue'};
% Unpack structure
unpack_struct(opt);   
   
    
%% Determine span of the plotted are according to the widest model (H=3)
opt = GBPRBM_Visible_Units_Span(opt);
opt.XLimLow  = opt.v(1).AxesLimMin;
opt.XLimHigh = opt.v(1).AxesLimMax;
opt.YLimLow  = opt.v(2).AxesLimMin;
opt.YLimHigh = opt.v(2).AxesLimMax;
opt.ZLimLow  = opt.v(3).AxesLimMin;
opt.ZLimHigh = opt.v(3).AxesLimMax;

figure('Name', 'GBPRBM: Probability of visible units for V=3', ...
       'Units', 'normalized', 'Position', [0.0198, 0.3021, 0.9488, 0.4036], ...
       'NumberTitle', 'Off');
for H=1:3
    opt.H = H;
    opt.W = W(:,1:H);
    opt.b_h =  (-opt.b_v'*Sinv*opt.W)';
    opt.h_axes = subplot(1,3,H);
    GBPRBM_Scatter_PDF_of_Visible_Units_3V(opt);   
    alpha(0.5);
    title(sprintf('p(v_1,v_2,v_3) with H=%i', H));
end
