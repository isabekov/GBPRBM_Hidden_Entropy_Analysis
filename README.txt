Gaussian-Bipolar Restricted Boltzmann Machines
Training GBPRBM model using Contrastive Divergence alogrithm for three-dimensional data (number of visible units is equal to 3).

Passing one to the function enables loading pretrained weights from LBG-like clustering alogirthm.
Pretrained weights are stored in "GeometryLBG.mat" file.
>> Synthetic_Data_3D_Train_GBPRBM(1)

Invoking the function without any arguments enables random initialization of the weights:
>> Synthetic_Data_3D_Train_GBPRBM

To obtain "GeometryLBG.mat" file, run 
>> GBPRBM_LBG_Pretraining
