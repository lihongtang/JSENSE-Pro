%% parfor version
function [recon,cmap]=adapt_array_3d_no_phase_corr(yn,bs,rn,norm)

% Reconstruction of array data and computation of coil sensitivities based
% on: a) Adaptive Reconstruction of MRI array data, Walsh et al. Magn Reson
% Med. 2000; 43(5):682-90 and b) Griswold et al. ISMRM 2002: 2410
%-------------------------------------------------------------------------
%   edit by Zekang Ding, 20200722
%   from 2D to 3D
%
%   Input:
%	yn: array data to be combined [nz, ny, nx, nc].
%	rn: data covariance matrix [nc, nc].
%	norm: =1, normalize image intensity
%
%	Output:
%	recon: reconstructed image [nz, ny, nx].
%	cmap: estimated coil sensitivity maps [nz, ny, nx, nc].
%--------------------------------------------------------------------------
% Ricardo Otazo
% CBI, New York University
%--------------------------------------------------------------------------
%

yn=permute(yn,[4,1,2,3]);
[nc,nz,ny,nx]=size(yn);
if nargin<4, norm=0; end
if nargin<3, rn=eye(nc);end%%%eyc(nc) nc unit matrix

% find coil with maximum intensity for phase correction
[mm,maxcoil]=max(sum(sum(sum(permute(abs(yn),[4 3 2 1])))));

bs1=bs;  %x-block size origin 4,10 ok
bs2=bs;  %y-block size
bs3=bs;  %z-block size
st=2;   %increase to set interpolation step size

% wsmall=zeros(nc,floor(nz./st),floor(ny./st),floor(nx./st));
% cmapsmall=zeros(nc,floor(nz./st),floor(ny./st),floor(nx./st));
wsmall=zeros(nc,floor(nz./st)*floor(ny./st)*floor(nx./st));
cmapsmall=zeros(nc,floor(nz./st)*floor(ny./st)*floor(nx./st));

invrn = inv(rn);

%%%%%%%%%%%%preparing loop list%%%%%%%%%%%%%% 
% x_index=[];y_index=[];z_index=[]; %edit by qk 20210304
% for x=st:st:nx
%     for y=st:st:ny
%         for z=st:st:nz
%             x_index=[x_index x];
%             y_index=[y_index y];
%             z_index=[z_index z];
%         end
%     end
% end
x=st:st:nx;
y=st:st:ny;
z=st:st:nz;
z_unit=z;
z_index=repmat(z_unit,1,length(x)*length(y));
y_unit=[];
for i=1:length(y)
    y_unit=[y_unit repmat(y(i),1,length(z))];
end
y_index=repmat(y_unit,1,length(x));
x_index=[];
for i=1:length(x)
    x_index=[x_index repmat(x(i),1,length(z)*length(y))];
end

parfor (i=1:length(x_index),30)
        %Collect block for calculation of blockwise values
        zmin1=max([z_index(i)-bs3./2 1]);
        ymin1=max([y_index(i)-bs1./2 1]);
        xmin1=max([x_index(i)-bs2./2 1]);
        % Cropping edges
        zmax1=min([z_index(i)+bs3./2 nz]);
        ymax1=min([y_index(i)+bs1./2 ny]);
        xmax1=min([x_index(i)+bs2./2 nx]);
        
        lz1=length(zmin1:zmax1);
        ly1=length(ymin1:ymax1);
        lx1=length(xmin1:xmax1);
        m1=reshape(yn(:,zmin1:zmax1,ymin1:ymax1,xmin1:xmax1),nc,lx1*ly1*lz1);
        
        m=m1*m1'; %signal covariance
        
        % eigenvector with max eigenvalue for optimal combination
        [e,v]=eig(invrn*m);
        
        v=diag(v);
        [mv,ind]=max(v);
        
        mf = e(:,ind);
        mf = mf/(mf'*invrn*mf);
        normmf=e(:,ind);
%         cmapsmall_no_corr_phase(:,i) = normmf;
%         wsmall_no_corr(:,i) = mf;
        % Phase correction based on coil with max intensity
        
%         mf = mf.*exp(-1i*angle(mf(maxcoil)));
%         normmf=normmf.*exp(-1i*angle(normmf(maxcoil)));
% %         wsmall(:,z./st,y./st,x./st)=mf;
%         cmapsmall(:,z./st,y./st,x./st)=normmf;
        wsmall(:,i)=mf;
        cmapsmall(:,i)=normmf;
end

% wsmall_temp=zeros(nc,floor(nz./st),floor(ny./st),floor(nx./st)); %edit by qk 20210304
% cmapsmall_temp=zeros(nc,floor(nz./st),floor(ny./st),floor(nx./st));
% for i=1:length(x_index)
%     wsmall_temp(:,z_index(i)./st,y_index(i)./st,x_index(i)./st)=wsmall(:,i);
%     cmapsmall_temp(:,z_index(i)./st,y_index(i)./st,x_index(i)./st)=cmapsmall(:,i);
% end
% wsmall=wsmall_temp;clear wsmall_temp
% cmapsmall=cmapsmall_temp;clear cmapsmall_temp
wsmall=reshape(wsmall,[nc,floor(nz./st),floor(ny./st),floor(nx./st)]);
cmapsmall=reshape(cmapsmall,[nc,floor(nz./st),floor(ny./st),floor(nx./st)]);
% wsmall_no_corr = reshape(wsmall_no_corr,[nc,floor(nz./st),floor(ny./st),floor(nx./st)]);
% cmapsmall_no_corr_phase=reshape(cmapsmall_no_corr_phase,[nc,floor(nz./st),floor(ny./st),floor(nx./st)]);
% cmapsmall_no_corr_phase=permute(cmapsmall_no_corr_phase,[2,3,4,1]);
% wsmall_no_corr = permute(wsmall_no_corr,[2,3,4,1]);
% cmapsmall = permute(cmapsmall,[2,3,4,1]);
% wsmall = permute(wsmall,[2,3,4,1]);
% opt.scale = [-pi,pi];
% figure;montagesc(squeeze(angle(cmapsmall_no_corr_phase(:,:,1,:))),opt);colormap('jet');
% opt.scale = [-pi,pi];
% figure;montagesc(squeeze(angle(cmapsmall(:,:,1,:))),opt);colormap('jet');
% opt.scale = [-pi,pi];
% figure;montagesc(squeeze(angle(wsmall_no_corr(:,:,1,:))),opt);colormap('jet');
% opt.scale = [-pi,pi];
% figure;montagesc(squeeze(angle(wsmall(:,:,1,:))),opt);colormap('jet');
% opt.scale = [];
% figure;montagesc(squeeze(abs(cmapsmall(:,:,1,:))),opt);
% opt.scale = [];
% figure;montagesc(squeeze(abs(wsmall(:,:,1,:))),opt);
% opt.scale = [];
% test = permute(yn,[2,3,4,1]);
% figure;montagesc(squeeze(abs(test(:,:,1,:))),opt);

% Interpolation of weights upto the full resolution
% Done separately for magnitude and phase in order to avoid 0 magnitude
% pixels between +1 and -1 pixels.
wfull=zeros(nc,nz,ny,nx);
cmap=zeros(nc,nz,ny,nx);

for i=1:nc
    wfull(i,:,:,:)=conj(imresize3(squeeze(abs(wsmall(i,:,:,:))),[nz,ny nx],'linear').*exp(1i.*imresize3(angle(squeeze(wsmall(i,:,:,:))),[nz,ny nx],'nearest')));
    cmap(i,:,:,:)=imresize3(squeeze(abs(cmapsmall(i,:,:,:))),[nz,ny nx],'linear').*exp(1i.*imresize3(squeeze(angle(cmapsmall(i,:,:,:))),[nz,ny nx],'nearest'));
end


recon=squeeze(sum(wfull.*yn));   %Combine coil signals.
% normalization proposed in the abstract by Griswold et al.
% if norm, recon=recon.*squeeze(sum(abs(cmap))).^2; end

cmap=permute(cmap,[2,3,4,1]);
% recon=permute(recon,[2,3,4,1]);

%%
% function [recon,cmap]=adapt_array_3d(yn,bs,rn,norm)
% 
% % Reconstruction of array data and computation of coil sensitivities based
% % on: a) Adaptive Reconstruction of MRI array data, Walsh et al. Magn Reson
% % Med. 2000; 43(5):682-90 and b) Griswold et al. ISMRM 2002: 2410
% %-------------------------------------------------------------------------
% %   edit by Zekang Ding, 20200722
% %   from 2D to 3D
% %
% %   Input:
% %	yn: array data to be combined [nz, ny, nx, nc].
% %	rn: data covariance matrix [nc, nc].
% %	norm: =1, normalize image intensity
% %
% %	Output:
% %	recon: reconstructed image [nz, ny, nx].
% %	cmap: estimated coil sensitivity maps [nz, ny, nx, nc].
% %--------------------------------------------------------------------------
% % Ricardo Otazo
% % CBI, New York University
% %--------------------------------------------------------------------------
% %
% 
% yn=permute(yn,[4,1,2,3]);
% [nc,nz,ny,nx]=size(yn);
% if nargin<4, norm=0; end
% if nargin<3, rn=eye(nc);end%%%eyc(nc) nc unit matrix
% 
% % find coil with maximum intensity for phase correction
% [mm,maxcoil]=max(sum(sum(sum(permute(abs(yn),[4 3 2 1])))));
% 
% bs1=bs;  %x-block size origin 4,10 ok
% bs2=bs;  %y-block size
% bs3=bs;  %z-block size
% st=2;   %increase to set interpolation step size
% 
% wsmall=zeros(nc,floor(nz./st),floor(ny./st),floor(nx./st));
% cmapsmall=zeros(nc,floor(nz./st),floor(ny./st),floor(nx./st));
% 
% invrn = inv(rn);
% 
% for x=st:st:nx
%     for y=st:st:ny
%         for z=st:st:nz
%         %Collect block for calculation of blockwise values
%         zmin1=max([z-bs3./2 1]);
%         ymin1=max([y-bs1./2 1]);
%         xmin1=max([x-bs2./2 1]);
%         % Cropping edges
%         zmax1=min([z+bs3./2 nz]);
%         ymax1=min([y+bs1./2 ny]);
%         xmax1=min([x+bs2./2 nx]);
%         
%         lz1=length(zmin1:zmax1);
%         ly1=length(ymin1:ymax1);
%         lx1=length(xmin1:xmax1);
%         m1=reshape(yn(:,zmin1:zmax1,ymin1:ymax1,xmin1:xmax1),nc,lx1*ly1*lz1);
%         
%         m=m1*m1'; %signal covariance
%         
%         % eigenvector with max eigenvalue for optimal combination
%         [e,v]=eig(invrn*m);
%         
%         v=diag(v);
%         [mv,ind]=max(v);
%         
%         mf=e(:,ind);
%         mf=mf/(mf'*invrn*mf);
%         normmf=e(:,ind);
%         
%         % Phase correction based on coil with max intensity
%         mf=mf.*exp(-1i*angle(mf(maxcoil)));
%         normmf=normmf.*exp(-1i*angle(normmf(maxcoil)));
%         
%         wsmall(:,z./st,y./st,x./st)=mf;
%         cmapsmall(:,z./st,y./st,x./st)=normmf;
%         end
%     end
% end
% 
% % Interpolation of weights upto the full resolution
% % Done separately for magnitude and phase in order to avoid 0 magnitude
% % pixels between +1 and -1 pixels.
% wfull=zeros(nc,nz,ny,nx);
% cmap=zeros(nc,nz,ny,nx);
% 
% for i=1:nc
%     wfull(i,:,:,:)=conj(imresize3(squeeze(abs(wsmall(i,:,:,:))),[nz,ny nx],'linear').*exp(1i.*imresize3(angle(squeeze(wsmall(i,:,:,:))),[nz,ny nx],'nearest')));
%     cmap(i,:,:,:)=imresize3(squeeze(abs(cmapsmall(i,:,:,:))),[nz,ny nx],'linear').*exp(1i.*imresize3(squeeze(angle(cmapsmall(i,:,:,:))),[nz,ny nx],'nearest'));
% end
% 
% 
% recon=squeeze(sum(wfull.*yn));   %Combine coil signals.
% % normalization proposed in the abstract by Griswold et al.
% % if norm, recon=recon.*squeeze(sum(abs(cmap))).^2; end
% 
% cmap=permute(cmap,[2,3,4,1]);
% % recon=permute(recon,[2,3,4,1]);
