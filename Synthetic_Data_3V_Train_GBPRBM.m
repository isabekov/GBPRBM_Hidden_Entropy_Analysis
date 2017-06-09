function opt = Synthetic_Data_3V_Train_GBPRBM(varargin)
% If "Enable_Visual_Debugging" is set to "true", then this function will
% create a directory called "ContrastiveDivergence" and store figures with
% model geometry evaluated after every iteration as PNG images. 
% If opt.Make_Video is set to one, then on Linux platform it will also
% call "ffmpeg" program to merge these images into an animation (MKV movie)
% to visualize evolution of the model parameters.

if ~isempty(varargin)
    if length(varargin) ==1
        LBG = varargin{1};
    end
else
    LBG = 0;
end
V = 3;
% Number of hidden units
H = 9;

% Standard deviations for visible units
sigma_v = 0.01*ones(V,1);
% Covariance matrix
S = diag(sigma_v.^2);
% Inverse of the "covariance matrix"
Sinv = diag(1./(sigma_v.^2));

Nv1 = 120;
Nv2 = 121;
Nv3 = 122;
Color = {[0, 0.6, 0.6], 'red', 'magenta', 'cyan', 'green'};

N_SpC = 2000;
NC = 8;
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
seed = 1;
data = GenerateData(V,NC,sigma_v(1)*4,N_SpC,seed);
data = data(randperm(NC*N_SpC),:);

features = data(1:ceil(NC*N_SpC/2),:);
N_Samples = size(features,1);
b_v = (mean(features))';
features_test = data(ceil(NC*N_SpC/2)+1:end,:);


if LBG ==0 
    W = 0.02*rand(V,H)-0.01;
    b_h = 0.02*rand(H,1);   
else
    load('GeometryLBG.mat');
end
% Mini-batch size (samples)
mBatch_Size = 100;
% 3D or imagesc?
plot_type = 'imagesc';
VRatio = [0.5 0.5 0.5];
FontSize = 10; %#ok<*NASGU>


Enable_Visual_Debugging = true;
vars = {'V', 'H', 'W', 'mBatch_Size', 'b_v', 'sigma_v', 'S'...
        'Sinv', 'plot_type', 'VRatio', 'b_h', 'Nv1', 'Nv2', 'Nv3', 'Color', ...
        'N_Samples', 'FontSize', 'Enable_Visual_Debugging'};
    % Pack structure
pack_struct('opt', vars); 

opt.Make_Video = true;

h_fig = figure('Name', 'Contrastive Divergence Training',  ...
       'Units', 'normalized', 'NumberTitle', 'Off', ...     
       'Position', [0.1947    0.2917    0.6911    0.4232]);
h_axes = subplot(121);
opt.h_axes = h_axes;
opt.LabelPosition = [0.2, 0.8, 0.7];
%opt = GBPRBM_Scatter_PDF_of_Visible_Units_3V(opt);
if opt.Enable_Visual_Debugging==true
    opt.h_fig = h_fig;    
    title('p(v_1,v_2,v_3) of the data');
else
    close(h_fig);
    pause(1);
end
 scatter3(data(:,1),data(:,2),data(:,3),'.y');
 xlabel('v_1');
 ylabel('v_2');
 zlabel('v_3');
set(gca, 'Projection', 'perspective');
set(gca, 'XLim', [0 1]);
set(gca, 'YLim', [0 1]);
set(gca, 'ZLim', [0 1]);

clear data;
% Original model
opt.model.W = W;
opt.model.b_v = b_v;
opt.model.b_h = b_h;
opt.model.sigma_v = sigma_v;
opt.model.Sinv = Sinv;
opt.model.S = S;
opt.XLimLow = 0;
opt.XLimHigh = 1;
opt.YLimLow = 0;
opt.YLimHigh = 1;
opt.ZLimLow = 0;
opt.ZLimHigh = 1;
%opt.model.Centroid = opt.Centroid;
%opt.model.Color = opt.Color;
opt.model.XLimLow  = opt.XLimLow;
opt.model.XLimHigh = opt.XLimHigh;
opt.model.YLimLow  = opt.YLimLow;
opt.model.YLimHigh = opt.YLimHigh;
opt.model.ZLimLow  = opt.ZLimLow;
opt.model.ZLimHigh = opt.ZLimHigh;
%features = [opt.v1' opt.v2' opt.v3'];     


% Contrastive divergence order
opt.CD_order = 5;
% Learning rate
opt.nu = 1e-5;
% Regularization constant
opt.lambda = 0.00000000;
% Momentum
opt.mu = 0.5;

opt.init.b_v = opt.model.b_v;
opt.init.W = opt.model.W;




opt.optim_HE = false;
if opt.optim_HE == true
    opt.init.b_h =  (-opt.init.b_v'*opt.Sinv*opt.init.W)';
else
    opt.init.b_h = rand(H,1);
end

opt.init.z_v = log(opt.model.sigma_v.^2);


% Number of mini-batches
M = ceil(N_Samples/mBatch_Size);
N_Epochs = 2;
opt.N_Epochs = N_Epochs;


% Number of decimal digits for printing iterations
N_Iter_Digits = log10(M);
if rem(N_Iter_Digits,1)==0
    % Consider logarithm of 1,10,100,1000 etc.
    N_Iter_Digits = N_Iter_Digits + 1;
else
    N_Iter_Digits = ceil(N_Iter_Digits);
end
opt.N_Iter_Digits = N_Iter_Digits;

% Number of decimal digits for printing epochs
N_Epoch_Digits = log10(N_Epochs);
if rem(N_Epoch_Digits,1)==0
    % Consider logarithm of 1,10,100,1000 etc.
    N_Epoch_Digits = N_Epoch_Digits + 1;
else
    N_Epoch_Digits = ceil(N_Epoch_Digits);
end
opt.N_Epoch_Digits = N_Epoch_Digits;


RMSE_Evolution = zeros(1,N_Epochs*M);

% Initial
opt.b_v = opt.init.b_v;
opt.W   = opt.init.W;
opt.b_h = opt.init.b_h;
opt.z_v = opt.init.z_v;

% Each cell Centroid{j} will contain a matrix of size (V x 2^j) with the
% centroid coordinates.
Centroid_Evolution = cell(H,1);
for j=1:H
    Centroid_Evolution{j} = zeros(V, 2^j, N_Epochs*M+1);
end
opt = GBPRBM_Visible_Units_Span(opt);
for j=1:H
    Centroid_Evolution{j}(:,:,1) = opt.Centroid{j};
end
b_v_Evolution = zeros(V,N_Epochs*M+1);
b_v_Evolution(:,1) = opt.b_v;

opt.DirSave = 'ContrastiveDivergence';
if exist(opt.DirSave,'dir') ~=7   
   mkdir(opt.DirSave); 
end

if opt.Enable_Visual_Debugging==true
    delete(fullfile('ContrastiveDivergence','*.png'));
    opt.plot_decision_regions = false;
    set(0, 'CurrentFigure', opt.h_fig);
    opt.h_axes = subplot(122);
    set(opt.h_axes, 'Projection', 'perspective');
end
    
for x = 1:N_Epochs
    opt.epoch = x;
    idx = randperm(size(features,1));
    data.features = features(idx,:);   

    % Train the GRBM model using Contrastive Divergence
    [opt] = GBPRBM_Train_3V(data, opt);   
    for j=1:H
        Centroid_Evolution{j}(:,:,(x-1)*M+2:x*M+1) = opt.Centroid_Evolution{j};
    end
    b_v_Evolution = zeros(V,N_Epochs*M+1);
    RMSE_Evolution((x-1)*M+1:x*M) = opt.RMSE_Evolution;
end
opt.final.b_v = opt.init.b_v;
opt.final.W = opt.init.W;
opt.final.b_h = opt.init.b_h;

N_Samples_Test = size(features_test,1);
H_Test = zeros(H, N_Samples_Test);
for i=1:N_Samples_Test
    [H_Test(:,i), P_h_m1_v, P_h_p1_v] = GBPRBM_Gibbs_Sampling_of_Hidden_Units(opt, features_test(i,:)');
end

figure('Name','Hidden Units', 'NumberTitle', 'Off');
subplot(131);
imagesc(H_Test);
xlabel('Test sample');
ylabel('Hidden vector');

subplot(132);
[h_count, pos] = hist(H_Test',2);
p_h_m = h_count(2,:)/N_Samples_Test;
p_h_p = h_count(1,:)/N_Samples_Test;
bar(p_h_m);
set(gca, 'XLim', [0.5, opt.H+0.5]);
set(gca, 'YLim',[0 1]);
set(gca, 'view', [90 90]);
ylabel('Probability');
Entropy = zeros(H,1);
for j=1:H
    if p_h_m(j)==1
       p_h_m(j)=1-eps; 
    end
    if p_h_p(j)==0
       p_h_p(j)=1+eps; 
    end
    Entropy(j) = -(  p_h_p(j).*log2(p_h_p(j)) + p_h_m(j).*log2(p_h_m(j))  );
end

subplot(133);
bar(Entropy);
set(gca, 'view', [90 90]);
set(gca, 'XLim', [0.5, opt.H+0.5]);
set(gca, 'YLim',[0 1]);
title(sprintf('Sum of Entropy of indiv. units=%f, Max is = %d,\nNormalized Entropy=%0.3f', sum(Entropy), opt.H, sum(Entropy)/opt.H));

%% Post-training procedures
figure('Name', 'Contrastive Divergence Results', 'NumberTitle', 'Off', ...
    'Units', 'normalized',...
   'Position', [ 0.2006    0.2148    0.5666    0.6523]); 

%% Subplot #2
opt.h_axes = subplot(222);
opt.LineStyle = '-';
% Plot final values
GBPRBM_Visible_Units_Plot_Geometry(opt);
set(gca, 'Box', 'on');
set(gca, 'Position', [0.5042    0.5788    0.3743    0.3545]);
hold on;    
% Plot initial values
opt.b_v = opt.init.b_v;
opt.W = opt.init.W;
opt.b_h = opt.init.b_h;
% Draw initial geometry in different color
opt.Color = {'black', 'green', 'cyan', 'blue'};
opt.LineStyle = '--';
%[opt] = GBPRBM_Visible_Units_Span(opt);
%GBPRBM_Visible_Units_Plot_Geometry(opt);

%set(gca, 'XLim', [opt.model.XLimLow opt.model.XLimHigh]);
%set(gca, 'YLim', [opt.model.YLimLow opt.model.YLimHigh]);
set(gca, 'XLim', [0 1]);
set(gca, 'YLim', [0 1]);
set(gca, 'ZLim', [0 1]);
plot(opt.b_v_Evolution(1,:), opt.b_v_Evolution(2,:), ...
        'LineStyle', ':');
hold on;
Centroid_Trajectories = zeros(N_Epochs*M+1,V,2^H);
for j=1:2^H
    Centroid_Trajectories(:,:,j) = shiftdim(Centroid_Evolution{H}(:,j,:),2);
end
for j=1:2^H
    plot3(Centroid_Trajectories(:,1,j), Centroid_Trajectories(:,2,j), Centroid_Trajectories(:,3,j), ...
        'LineStyle', ':');
    hold on;
end
daspect([1 1 1]);
title('(b) Evolution of the model geometry')
xlabel('v_1');
ylabel('v_2');
view( -37.5000,30);
%set(gca, 'Projection', 'perspective');

%% Subplot #1
opt.h_axes = subplot(221);
opt.LabelPosition = [0.2, 0.8, 0.7];
opt.v(1).AxesLimMax = opt.model.XLimHigh;
opt.v(1).AxesLimMin = opt.model.XLimLow;
opt.v(2).AxesLimMax = opt.model.YLimHigh;
opt.v(2).AxesLimMin = opt.model.YLimLow;
% Original model
opt.W   = opt.model.W;
opt.b_v = opt.model.b_v;
opt.b_h = opt.model.b_h;
%opt.Centroid = opt.model.Centroid;
% Area of each marker
S = 1.5*ones(1,N_Samples);
scatter3(features(:,1),features(:,2),features(:,3), 'filled', 'Marker', '.', 'MarkerFaceColor', 'blue', ...
                     'MarkerEdgeColor', 'blue');
xlabel('v_1');
ylabel('v_2');
zlabel('v_3');
%opt.Color = opt.model.Color;
%[opt] = GBPRBM_Visible_Units_Span(opt);
%GBPRBM_Visible_Units_Plot_Geometry(opt);
%GBPRBM_Plot_Decision_Region_for_Visible_Units(opt);
set(gca, 'XLim', [opt.model.XLimLow opt.model.XLimHigh]);
set(gca, 'YLim', [opt.model.YLimLow opt.model.YLimHigh]);
set(gca, 'ZLim', [opt.model.ZLimLow opt.model.ZLimHigh]);
set(gca, 'Box', 'on');
daspect([1 1 1]);
title('(a) p(v_1,v_2,v_3) of the data');
% set(gca, 'Position', [0.0552    0.5788    0.3743    0.3545]);
%set(gca, 'Projection', 'perspective');

%% Subplot #3
subplot(212);
plot(1:M*N_Epochs,RMSE_Evolution);
%set(gca, 'Position', [0.0552    0.1098    0.3743    0.3545]);
axis tight;
title('(c) Evolution of RMSE (training)');
xlabel('Iteration');
ylabel('RMSE');


if opt.Make_Video == true    
    if strncmp(computer('arch'),'glnx',4)        
        VideoFileTemplate = ['ContrastiveDivergence_Epoch_%0' num2str(N_Epoch_Digits)  'd.mkv'];
        for epoch=1:N_Epochs            
            FileName = fullfile(sprintf(VideoFileTemplate, epoch));
            delete(FileName);
            str = ['cd ' opt.DirSave ' && ffmpeg -y -nostdin -framerate 5 -start_number 0 -i CD_' sprintf(['%0' num2str(N_Epoch_Digits) 'd'],epoch) 'e_%0' num2str(N_Iter_Digits) 'di.png -c:v libx264 -r 3 ' FileName];       
            disp('================= Executing =================');
            disp(str);
            disp('=============================================');
            unix(str,'-echo');
        end
        % List of files to be merged
        MergeList = 'MergeList.txt';
        fid = fopen(fullfile(opt.DirSave, MergeList), 'w+');
        for epoch=1:N_Epochs
            FileName = sprintf(VideoFileTemplate, epoch);
            fprintf(fid, ['file ' FileName '\n']);
        end
        fclose(fid);
        str = ['cd ' opt.DirSave ' && ffmpeg -y -nostdin -f concat -i ' MergeList ' -c copy ContrastiveDivergence_Evolution.mkv'];
        disp('================= Executing =================');
        disp(str);
        disp('=============================================');
        unix(str,'-echo');
        delete('ContrastiveDivergence_Epoch_*.mp4')
    end
end    