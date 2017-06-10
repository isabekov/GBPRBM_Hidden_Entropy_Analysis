function Download_MNIST_Database
BaseURL = 'http://yann.lecun.com/exdb/mnist/';
fprintf('MNIST training and test sets will be downloaded from %s\n', BaseURL);
% =================== Training Images =======================
FileName = 'train-images-idx3-ubyte.gz';

if(~exist(FileName,'file'))
    urlwrite([BaseURL, FileName], FileName);
end
FileName = gunzip(FileName);
fid = fopen(FileName{1},'r','b');   % big-endian
Magic_Number = fread(fid,1,'int32'); % Magic number
if(Magic_Number~=2051) 
    display('Error: cannot find magic number');
    return;
end
N_Images = fread(fid,1,'int32');
% Image height
N_Rows = fread(fid,1,'int32');  
% Image width
N_Columns = fread(fid,1,'int32');

Image = cell(1,N_Images);
ImageVector = zeros(N_Rows*N_Columns, N_Images);
for k=1:N_Images
    Image{k} = uint8(fread(fid,[N_Rows N_Columns],'uchar')); 
    % Transpose before reshaping
    ImageVector(:,k) =reshape(Image{k}', N_Rows*N_Columns, 1);
end
fclose(fid);

% ================== Training labels ========================
FileName = 'train-labels-idx1-ubyte.gz';

if(~exist(FileName,'file'))
    urlwrite([BaseURL, FileName], FileName);
end
FileName = gunzip(FileName);
fid = fopen(FileName{1},'r','b');   % big-endian
Magic_Number = fread(fid,1,'int32'); 
if(Magic_Number~=2049) 
    display('Error: cannot find magic number');
    return;
end
N_Labels = fread(fid, 1, 'int32');
tmp = fread(fid,N_Labels, 'uint8');
% Load all labels
Label = uint8(tmp');   
% Vector of labels, i.e. integer "2" corresponds to [0 0 1 0 0 0 0 0 0 0]'
% Maximal range is 0 to 9
LabelVector = zeros(10, N_Labels);
for k = 1:N_Labels
   LabelVector(double(Label(k) + 1),k) = 1; 
end
% LabelVector = logical(LabelVector);
fclose(fid);
%save('MNIST_Train_Full.mat', 'Image', 'Label');
trainData = ImageVector'/255;
trainLabels = LabelVector'; 
save('MNIST_Train_Medal_Normalized.mat', 'trainData', 'trainLabels');


% ======================= Test Images ===============================
FileName = 't10k-images-idx3-ubyte.gz';

if(~exist(FileName,'file'))
    urlwrite([BaseURL, FileName], FileName);    
end
FileName = gunzip(FileName);
fid = fopen(FileName{1},'r','b');  
Magic_Number = fread(fid,1,'int32');    
if(Magic_Number~=2051) 
    display('Error: cannot find magic number');
    return;
end
N_Images = fread(fid,1,'int32');  
N_Rows = fread(fid,1,'int32');   
N_Columns = fread(fid,1,'int32');   

Image = cell(1,N_Images);
ImageVector = zeros(N_Rows*N_Columns, N_Images);
for k=1:N_Images
    Image{k} = uint8(fread(fid,[N_Rows N_Columns],'uchar'));
    % Transpose before reshaping
    ImageVector(:,k) = reshape(Image{k}', N_Rows*N_Columns, 1);
end
fclose(fid);

% ======================= Test Labels  =============================
FileName = 't10k-labels-idx1-ubyte.gz';

if(~exist(FileName,'file'))
    urlwrite([BaseURL, FileName], FileName);
end
FileName = gunzip(FileName);
fid = fopen(FileName{1},'r','b');
Magic_Number = fread(fid,1,'int32');
if(Magic_Number~=2049)
    display('Error: cannot find magic number');
    return;
end
N_Labels = fread(fid,1,'int32');
tmp = fread(fid,N_Labels, 'uint8');
Label = uint8(tmp');   
LabelVector = zeros(10, N_Labels);
for k = 1:N_Labels
   LabelVector(double(Label(k) + 1),k) = 1; 
end
testData = ImageVector'/255; %#ok<*NASGU>
testLabels = LabelVector'; 
save('MNIST_Test_Medal_Normalized.mat', 'testData', 'testLabels');
fclose(fid);