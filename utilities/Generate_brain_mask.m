%% Step 3: Motion compensation
load(fullfile(data_path,'data_infor.mat'));
all_sensitivity = reshape(A,dim_maps);
  i_subj = 5;
use_slice = [1:12];
load(fullfile(data_path,'data_infor.mat'));
Temp_test_subj = i_subj;
temp_filename = SIEMENS_Aera_1_494_20_TH{Temp_test_subj,2};
temp_scanner_name = 'SIEMENS_Aera_1_494_20_TH';
temp_filename = temp_filename{1,1};
data_path_h5 = '/home/tanglihong/ML_PI/fastMRI/brain_all';
h5disp(fullfile(data_path_h5,temp_filename))
Groundtruth  = h5read(fullfile(data_path_h5,temp_filename),'/reconstruction_rss');
Groundtruth = Groundtruth(:,:,use_slice);
K_space  = h5read(fullfile(data_path_h5,temp_filename),'/kspace');
K_space = K_space.r+sqrt(-1).* K_space.i;
K_space = K_space(:,1:2:end,:,use_slice);
figure;montagesc(Groundtruth)
Brain_mask_binary = imbinarize(Groundtruth,mean(Groundtruth(:)));
se = strel('disk',1);
Brain_mask = imclose(Brain_mask_binary,se);
se = strel('disk',1);
Brain_mask = imerode(Brain_mask,se);
se = strel('disk',30);
Brain_mask = imclose(Brain_mask,se);
figure;overlay(Brain_mask_binary,Brain_mask)
fprintf('subj %d mask is save\n',i_subj);
brain_mask_all(:,:,:,i_subj) = Brain_mask;