function Synthetic_Data_Generate
seed = 1;
sigma_v = 0.1;
DirSave = 'Synth_Data';
if exist(DirSave,'dir') ~=7
   mkdir(DirSave);
end
N_Samples = 60000;
disp('Generating synthetic data');
for V= 2.^(3:7)
    fprintf('V = %d\n',V);
    % Number of centroids
    NC = round(50*1.05^V);
    % Number of samples per centroid
    N_SpC = round(N_Samples/NC);
    data = GenerateData(V, NC, sigma_v, N_SpC, seed);
    % Shuffle data samples
    data = data(randperm(NC*N_SpC),:);
    % Divide data into training set (70%) and test set (30%)
    N_Samples_Train = round(0.7*NC*N_SpC);
    trainData = data(1:N_Samples_Train,:); %#ok<*NASGU>
    testData  = data(N_Samples_Train+1:end,:);
    save(fullfile(DirSave, sprintf('Synth_Data_Train_V=%d.mat',V)), 'trainData');
    save(fullfile(DirSave, sprintf('Synth_Data_Test_V=%d.mat',V)),  'testData');
end