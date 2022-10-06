function resL1SENSE = SparseSENSE3D(sen,underSampleMask,underSampleData,resinitial,lambda)
dim_spatial = size(underSampleData);
% parameters for L1-reconstruction with splitting
nIterCG = 5;       % number of CG iterations for the PI part
nIterSplit = 15;    % number of splitting iterations for CS part
splitWeight = 0.4;  % reasonable value


XOP = Wavelet('Daubechies_TI',4,6);
FT = p3DFT(underSampleMask,dim_spatial);
% zf_recon = FT'*underSampleData;
% figure;montagesc(my_rss(zf_recon));colorbar;
% create and ESPIRiT operator
SENSEOP = SENSE3D(sen);

%% Reconsturctions
% ESPIRiT CG reconstruction with soft-sense and 1 sets of maps



   if lambda==0
       disp('Performing SENSE reconstruction from SOS maps')
   AHA = @(x)vec(SENSEOP'*(FT'*(FT*(SENSEOP*x))));
   AHb = vec(SENSEOP'*( FT'*underSampleData));
   [resL1SENSE,~] = pcg(AHA,AHb,[],[],[],[],resinitial(:));
   resL1SENSE = reshape(resL1SENSE,dim_spatial(1:3));
   % maxit = 500;
   % tol = 1e-6;
   % resL1SENSE = mypcg(AHA,AHb,tol,maxit,[],[],resinitial(:));
   
   % resL1SENSE = reshape(resL1SENSE,dim_spatial(1:3));
%    figure;montagesc(angle(resL1SENSE));colormap('jet')
% figure;montagesc(abs(resL1SENSE));colormap('jet')
   else
disp('Performing L1-SENSE reconstruction from SOS maps')
tic
[resL1SENSE] = cgL1SENSE3D(underSampleData, double(resinitial), FT, SENSEOP, nIterCG,XOP,lambda,splitWeight,nIterSplit);
   
toc
end
