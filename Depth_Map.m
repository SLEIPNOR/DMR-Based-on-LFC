%% Image format convert from 12 bit [168,4095] to double (0,1)
clear
% white image read 
WhiteImg = LFReadLFP('WI.lfp');
WhiteImg = cast(WhiteImg.RawImg,'double');

% LF image read
LFP = LFReadLFP('raw.lfp');
LensletImg = cast(LFP.RawImg,'double');

% pixel value limits
black = LFP.Metadata.image.rawDetails.pixelFormat.black.gr;
white = LFP.Metadata.image.rawDetails.pixelFormat.white.gr;

% normalize to [0,1]
WhiteImg = (WhiteImg-black)./(white-black);
LensletImg = (LensletImg-black)./(white-black);

WhiteImg = imadjust(WhiteImg);
LensletImg = imadjust(LensletImg);

% demosaic
WhiteImg = cast(WhiteImg.*double(intmax('uint16')), 'uint16');
WhiteImg = demosaic(WhiteImg, LFP.DemosaicOrder);

LensletImg = cast(LensletImg.*double(intmax('uint16')), 'uint16');
LensletImg = demosaic(LensletImg, LFP.DemosaicOrder);

% rgb to gray
WhiteImg = rgb2gray(WhiteImg);
LensletImg = rgb2gray(LensletImg);
WhiteImg = cast(WhiteImg,'double');
LensletImg = cast(LensletImg,'double');

% normalize to [0,1]
WhiteImg = (WhiteImg-0)./(65535-0);
LensletImg = (LensletImg-0)./(65535-0);

%*******************************************************************************
% Grid build

GridPT = GridBuildTools;

GridPT = ProEst (GridPT,WhiteImg);

GridCoords = GridBuild(GridPT);


% Light field process


YSpacing = GridPT.YSpacing;
SpaceRatio = YSpacing/GridPT.YSpacing;

XSpacing = SpaceRatio*GridPT.XSpacing;

YMax = floor((GridPT.IM_size(2)-GridPT.Offset(2))/(cos(GridPT.Rotation)*GridPT.YSpacing))+1 ;
XMax = floor((GridPT.IM_size(1)-GridPT.Offset(1))/(cos(GridPT.Rotation)*GridPT.XSpacing))+1 ; 

u_list = -YSpacing/2:GridPT.YSpacing/(ceil(GridPT.YSpacing)+mod(ceil(GridPT.YSpacing+1),2)-1):YSpacing/2;
v_list = -YSpacing/2:GridPT.YSpacing/(ceil(GridPT.YSpacing)+mod(ceil(GridPT.YSpacing+1),2)-1):YSpacing/2;
s_list = 0:1:YMax-1;
t_list = 0:1:XMax-1;

% LF for lensletimage(LLI)
LF_LLI = zeros(length(t_list), length(s_list), length(v_list), length(u_list));
[t,s,v,u] = ndgrid(t_list, s_list, v_list, u_list);

R = sqrt(u.^2 + v.^2);
outlierID = find(R >= YSpacing/2); % >=

Y_coord = u+s*YSpacing;
X_coord = v+t*XSpacing;

Rot_matrix = [cos(GridPT.Rotation),-sin(GridPT.Rotation);sin(GridPT.Rotation),cos(GridPT.Rotation)];
Y_coord(GridPT.FirstPosShiftRow:2:end,:,:,:) = Y_coord(GridPT.FirstPosShiftRow:2:end,:,:,:) + YSpacing/2;
SpPts = [SpaceRatio,0,0;0,SpaceRatio,0;0,0,1]^-1*[X_coord(:),Y_coord(:),ones(size(X_coord(:)))]';
SpPts = [Rot_matrix,GridPT.Offset']*SpPts;
SpPts = SpPts(1:2,:)';

Interpolant = griddedInterpolant(LensletImg,'linear','none');


SpPts_val = Interpolant(SpPts);
SpPts_val(isnan(SpPts_val)) = 0; % out of 3280*3280
SpPts_val(outlierID) = 0; % out of circle angle space r ~= pitch
LF_LLI = reshape(SpPts_val, size(LF_LLI));

% standardize LF_LLI to LF_LLIstd
fprintf('Standardize Light Field:')
st = size(LF_LLI,1);
ss = size(LF_LLI,2);
sv = size(LF_LLI,3);
su = size(LF_LLI,4);

sz = max(st,ss);
sz = sz+mod(sz-1,2);

LF_LLIstd = zeros(sz, sz, sv, su);

count = 0;

for nv=1:sv
    
    for nu=1:su
        
        img_corrected = zeros(st,2*ss);
        img = squeeze(LF_LLI(:,:,nv,nu));    

        img_shift = img(GridPT.FirstPosShiftRow:2:end,:);
        img_corrected(GridPT.FirstPosShiftRow:2:end,2:2:end) = img_shift;

        if GridPT.FirstPosShiftRow ==1

            img_noshift = img(GridPT.FirstPosShiftRow+1:2:end,:);
            img_corrected(GridPT.FirstPosShiftRow+1:2:end,1:2:end) = img_noshift;

        else

            img_noshift = img(GridPT.FirstPosShiftRow-1:2:end,:);
            img_corrected(GridPT.FirstPosShiftRow-1:2:end,1:2:end) = img_noshift;

        end

% The first few black image cannot be filled becasue all pixels are zero
        if ~all(img_corrected(:) == 0)

            img_corrected = regionfill(img_corrected,img_corrected==0);
            
        end
               
        
        img_corrected = imresize(img_corrected, [sz,sz]);

        LF_LLIstd(:,:,nv,nu) = img_corrected;
        
        fprintf(repmat('\b',1,count))
        count=fprintf('%.2f%%\n',100*nv/sv);
    end
end

LF_LLIstd = LF_LLIstd(:, :, end:-1:1, end:-1:1);


% Fourier slice

% build coordinate main lens diameter is 25e3 um sensor size is 3280*1.4=4592 um
LF = LF_LLIstd;
LF_size = size(LF);
D_Mainlens_org = 25e3;
% D_Mainlens_shrink = 25e3/ShrinkScale; % divide ShrinkScale if use LF_resample
D_Mainlens = D_Mainlens_org;
D_sensor = GridPT.IM_size*1.4;
t_Rg = linspace(-D_sensor(1)/2,D_sensor(1)/2,LF_size(1));
s_Rg = linspace(-D_sensor(1)/2,D_sensor(1)/2,LF_size(2));
v_Rg = linspace(-D_Mainlens/2,D_Mainlens/2,LF_size(3));
u_Rg = linspace(-D_Mainlens/2,D_Mainlens/2,LF_size(4));

t_sr = 1/(t_Rg(2)-t_Rg(1));
s_sr = 1/(s_Rg(2)-s_Rg(1));
v_sr = 1/(v_Rg(2)-v_Rg(1));
u_sr = 1/(u_Rg(2)-u_Rg(1));
% [t_mesh,s_mesh,v_mesh,u_mesh] = ndgrid(t_Rg,s_Rg,v_Rg,u_Rg);

% Fourier domain
LF_FFT = fftshift(fftn(ifftshift(LF)));
LF_FFTABS= log(abs(LF_FFT)+1); 
LF_FFTABS = (LF_FFTABS-min(LF_FFTABS(:)))./(max(LF_FFTABS(:))-min(LF_FFTABS(:)));
% f range (bandwidth)
kt_Rg = linspace(-t_sr/2, t_sr/2, length(t_Rg));
ks_Rg = linspace(-s_sr/2, s_sr/2, length(s_Rg));
kv_Rg = linspace(-v_sr/2, v_sr/2, length(v_Rg));
ku_Rg = linspace(-u_sr/2, u_sr/2, length(u_Rg));

[kt_mesh,ks_mesh,kv_mesh,ku_mesh] = ndgrid(kt_Rg,ks_Rg,kv_Rg,ku_Rg);

% interpolation
FFT_Interpolant = griddedInterpolant(kt_mesh,ks_mesh,kv_mesh,ku_mesh,LF_FFT,'linear','none');

%***********************
% Precision Range      *
%***********************
% Refocus range in t v space
alhpa_max_t = t_sr/(t_sr-v_sr);
alhpa_min_t = t_sr/(t_sr+v_sr);

% Refocus range in s u space
alhpa_max_s = s_sr/(s_sr-u_sr);
alhpa_min_s = s_sr/(s_sr+u_sr);

% Refocus range 
alhpa_max = min(alhpa_max_t,alhpa_max_s);
alhpa_min = max(alhpa_min_t,alhpa_min_s);

% frequency domain refocusing

format long
f_M = 50e3;
F = 53e3;
AN = f_M/D_Mainlens;
count = 0;
scale_factor = 1.007;
alpha = alhpa_min/scale_factor;
% scale_factor_alpha = 1;
% for alpha=scale_factor*alhpa_max:-alpha_unit/scale_factor_alpha:alhpa_min/scale_factor % alpha=scale_factor*alhpa_max:-alpha_unit:alhpa_min/scale_factor
% for alpha = 1.006
%     alpha = 0.995;
while(alpha<=scale_factor *alhpa_max)
    alpha_M(count+1) = alpha;
    beta = [alpha,0,1-alpha,0;0,alpha,0,1-alpha;0,0,1,0;0,0,0,1];
    
    % ktRg_beta = linspace(max(min(kv_Rg)/abs(1-alpha),min(kt_Rg)/alpha),min(max(kv_Rg)/abs(1-alpha),max(kt_Rg)/alpha),length(t_Rg));
    % ksRg_beta = linspace(max(min(ku_Rg)/abs(1-alpha),min(ks_Rg)/alpha),min(max(ku_Rg)/abs(1-alpha),max(ks_Rg)/alpha),length(s_Rg));
    
    ktRg_beta = kt_Rg ;
    ksRg_beta = ks_Rg;
    [ktb_mesh,ksb_mesh,kvb_mesh,kub_mesh] = ndgrid(ktRg_beta,ksRg_beta,0,0);
    
    mesh_slice = beta'*[ktb_mesh(:),ksb_mesh(:),kvb_mesh(:),kub_mesh(:)]';
    
    LF_FFT_Slice = FFT_Interpolant(mesh_slice');
    LF_FFT_Slice(isnan(LF_FFT_Slice)) = 0;
    LF_FFT_Slice = reshape(LF_FFT_Slice,[length(t_Rg),length(s_Rg)]);
    
    LF_inv = fftshift(ifftn(ifftshift(LF_FFT_Slice)));
    LF_inv = (LF_inv-min(LF_inv(:)))./(max(LF_inv(:))-min(LF_inv(:)));
    LF_inv = abs(LF_inv);
    LF_inv = imadjust(LF_inv);
    LF_inv = imadjust(LF_inv,[0,1],[0.12,1]);
    LF_FFT_Slice = log(abs(LF_FFT_Slice)+1);
    LF_FFT_Slice = (LF_FFT_Slice-min(LF_FFT_Slice(:)))./(max(LF_FFT_Slice(:))-min(LF_FFT_Slice(:)));
    
%     % sharp
%     laplacianFilter = [0 -1 0; -1 4 -1; 0 -1 0];
%     LF_inv_filtered = imfilter(LF_inv, laplacianFilter);
%     LF_inv = LF_inv + LF_inv_filtered;

%     figure(1)
%     subplot(2, 2, 1);
%     imshow(LF_FFT_Slice)
%     colormap jet
%     title(sprintf('Fourier Slice (alpha = %.4f)',alpha));
% 
%     subplot(2, 2, 2);
%     imshow(rot90(LF_inv))
%     title('Refocused Image');
% 
%     subplot(2, 2, 3);
%     plot(mesh_slice(1,:),1e2*mesh_slice(3,:))
%     hold on 
%     plot([-1 1], [0 0], 'k--') 
%     plot([0 0], [-1, 1], 'k--')
%     hold off
%     rectangle('Position', [min(kt_Rg), 1e2*min(kv_Rg), t_sr, 1e2*v_sr],'EdgeColor', 'red')
%     xlim([-0.1, 0.1])   
%     ylim([-1e-1, 1e-1])  
%     axis equal  
%     xlabel('kt');
%     ylabel('kv(1e-2)');
%     title('kt,kv Space in Fourier Domain');
% 
%     subplot(2, 2, 4);
%     plot(mesh_slice(2,:),1e2*mesh_slice(4,:))
%     hold on 
%     plot([-1 1], [0 0], 'k--') 
%     plot([0 0], [-1, 1], 'k--')
%     hold off
%     rectangle('Position', [min(ks_Rg), 1e2*min(ku_Rg), s_sr, 1e2*u_sr],'EdgeColor', 'red')
%     xlim([-0.1, 0.1])   
%     ylim([-1e-1, 1e-1]) 
%     axis equal  
%     xlabel('ks');
%     ylabel('ku(1e-2)');
%     title('ks,ku Space in Fourier Domain');

    % crop & save refocus img
    res1 = size(LF_inv,1);
    res2 = size(LF_inv,2);
    cropval = 0.05;
    LF_inv = imcrop(LF_inv,[res1*cropval res2*cropval res1*(1-2*cropval) res2*(1-2*cropval)]);
%     imwrite(LF_inv, sprintf('refocused\\%d.png', count),'png');

    % count depth label
    count = count+1;
%     pause(0.05);

    %***********************
    % Refocus Precision    *
    %***********************
    
    % Precision limit in t v  space
    
    d_alpha_t = (2*(alpha^2))/(t_sr*D_Mainlens-2*alpha);
    
    % Precision limit in s u space
    
    d_alpha_s = (2*(alpha^2))/(s_sr*D_Mainlens-2*alpha);
    
    % DoF

    Ro = (alpha^2* F^2)/(alpha* F - f_M);
    h = (2*f_M^2*AN*(1/t_sr)*Ro^2)/(f_M^4);
    
    
    d_alpha_h = (h*alpha^2*F^2-2*alpha*h*f_M*F+h*f_M^2)/(F*(f_M^2-h*alpha*F+h*f_M));
    
    % Choose minimal unit value as Precision limit
    d_alpha = max(min(d_alpha_t,d_alpha_s),d_alpha_h);

    alpha = alpha+d_alpha;
    IM_mat(:,:,count) = rot90(LF_inv);
end

% template image by Fourier slice

% resample LF_LLIstd
ShrinkScale = 1e10; 

OrigSize = size(LF_LLIstd, 3);
OrigVec = floor((-(OrigSize-1)/2):((OrigSize-1)/2));
NewVec = OrigVec ./ ShrinkScale;


LF_resample = ipermute(LF_LLIstd,[3,2,1,4]);

LF_resample(:,:,:,:) = interp1(OrigVec, LF_resample(:,:,:,:), NewVec);


LF_resample = ipermute(LF_resample,[3,2,1,4]);
LF_resample(isnan(LF_resample)) = 0;


OrigSize = size(LF_LLIstd, 4);
OrigVec = floor((-(OrigSize-1)/2):((OrigSize-1)/2));
NewVec = OrigVec ./ ShrinkScale;


LF_resample = ipermute(LF_resample,[4,2,3,1]);

LF_resample(:,:,:,:) = interp1(OrigVec, LF_resample(:,:,:,:), NewVec);


LF_resample = ipermute(LF_resample,[4,2,3,1]);
LF_resample(isnan(LF_resample)) = 0;


% build coordinate main lens diameter is 25e3 um sensor size is 3280*1.4=4592 um
LF = LF_resample;
LF_size = size(LF);
D_Mainlens_shrink = 25e3/ShrinkScale; % divide ShrinkScale if use LF_resample
D_Mainlens = D_Mainlens_shrink;
D_sensor = GridPT.IM_size*1.4;
t_Rg = linspace(-D_sensor(1)/2,D_sensor(1)/2,LF_size(1));
s_Rg = linspace(-D_sensor(1)/2,D_sensor(1)/2,LF_size(2));
v_Rg = linspace(-D_Mainlens/2,D_Mainlens/2,LF_size(3));
u_Rg = linspace(-D_Mainlens/2,D_Mainlens/2,LF_size(4));

t_sr = 1/(t_Rg(2)-t_Rg(1));
s_sr = 1/(s_Rg(2)-s_Rg(1));
v_sr = 1/(v_Rg(2)-v_Rg(1));
u_sr = 1/(u_Rg(2)-u_Rg(1));

% Fourier domain
LF_FFT = fftshift(fftn(ifftshift(LF)));
LF_FFTABS= log(abs(LF_FFT)+1); 
LF_FFTABS = (LF_FFTABS-min(LF_FFTABS(:)))./(max(LF_FFTABS(:))-min(LF_FFTABS(:)));
% f range (bandwidth)
kt_Rg = linspace(-t_sr/2, t_sr/2, length(t_Rg));
ks_Rg = linspace(-s_sr/2, s_sr/2, length(s_Rg));
kv_Rg = linspace(-v_sr/2, v_sr/2, length(v_Rg));
ku_Rg = linspace(-u_sr/2, u_sr/2, length(u_Rg));

[kt_mesh,ks_mesh,kv_mesh,ku_mesh] = ndgrid(kt_Rg,ks_Rg,kv_Rg,ku_Rg);

% interpolation
FFT_Interpolant = griddedInterpolant(kt_mesh,ks_mesh,kv_mesh,ku_mesh,LF_FFT,'linear','none');

format long

alpha = 1;
beta = [alpha,0,1-alpha,0;0,alpha,0,1-alpha;0,0,1,0;0,0,0,1];

ktRg_beta = kt_Rg ;
ksRg_beta = ks_Rg;
[ktb_mesh,ksb_mesh,kvb_mesh,kub_mesh] = ndgrid(ktRg_beta,ksRg_beta,0,0);

mesh_slice = beta'*[ktb_mesh(:),ksb_mesh(:),kvb_mesh(:),kub_mesh(:)]';

LF_FFT_Slice = FFT_Interpolant(mesh_slice');
LF_FFT_Slice(isnan(LF_FFT_Slice)) = 0;
LF_FFT_Slice = reshape(LF_FFT_Slice,[length(t_Rg),length(s_Rg)]);

LF_inv = fftshift(ifftn(ifftshift(LF_FFT_Slice)));
LF_inv = (LF_inv-min(LF_inv(:)))./(max(LF_inv(:))-min(LF_inv(:)));
LF_inv = abs(LF_inv);
LF_inv = imadjust(LF_inv);
LF_inv = imadjust(LF_inv,[0,1],[0.12,1]);

% figure(2)
% imshow(LF_inv)
LF_inv = imcrop(LF_inv,[res1*cropval res2*cropval res1*(1-2*cropval) res2*(1-2*cropval)]);
% imwrite(LF_inv, sprintf('ImageTemplate.png'),'png');
IM_template = rot90(LF_inv);

% FoV normalization

% reficusing image norm
IM_mat = FoV_normalisation(IM_mat);
% montage({IM_mat(:,:,1),IM_mat(:,:,2),IM_mat(:,:,3),IM_mat(:,:,4),IM_mat(:,:,5)},"Size",[3 2])

% template image normalization
template_fix(:,:,1) =  IM_mat(:,:,round(size(IM_mat,3)/2));
template_fix(:,:,2) =  IM_template;
template_fix = FoV_normalisation(template_fix);

% cropping
IM_mat(:,:,size(IM_mat,3)+1) = template_fix(:,:,2);
IM_mat = ZeroCrop(IM_mat);
% montage({IM_mat(:,:,1),IM_mat(:,:,2),IM_mat(:,:,3),IM_mat(:,:,4),IM_mat(:,:,5)},"Size",[3 2])

IM_template = IM_mat(:,:,end);
IM_mat(:,:,end)=[];
%% view fov normalised image
figure(2)
for i = 1:size(IM_mat,3)
    imshow(IM_mat(:,:,i))
    pause(0.05)
end
%% FV calculation in each patch in FD
laplacianFilter = [0 -1 0; -1 4 -1; 0 -1 0];
hf_ratio=0.5;


% IMT_filtered1 = imfilter(IM_template, laplacianFilter);
% IMT_filtered2 = imfilter(IM_template, laplacianFilter');
% IM_templatehf = IM_template+hf_ratio*(IMT_filtered1+IMT_filtered2);

IM_e = edge(imgaussfilt(IM_template,1.5),'canny',1.3e-1);
figure(2)

imshow(IM_template)
title('Template Image', 'FontSize', 16);
figure(3)

imshow(IM_e);
title('Edge Detection', 'FontSize', 16);


s = size(IM_template);
x = s(2);
y = s(1);

patchsize= 9;
fx_mat_tensor = zeros(y,x,size(IM_mat,3));

numImages = size(IM_mat, 3);

for img_no = 1:numImages
   
    IM = IM_mat(:,:,img_no);
    IM = imgaussfilt(IM,2);
%     IM_filtered = imfilter(IM, laplacianFilter);
%     
%     IM = IM+hf_ratio*IM_filtered;
    IM = imsharpen(IM);
    fx_pixel = LAPD(IM);

 

    
    IM_e_patch = cut_patch(IM_e,patchsize);

    [pos_row,pos_col] = find(IM_e_patch==0);
    fx_mat = cut_patch(fx_pixel,patchsize);

    fx_mat(sub2ind(size(fx_mat), pos_row, pos_col)) = 0;
%     IM_e_patch(sub2ind(size(fx_mat), pos_row, pos_col)) = 0;
    fx_mat_tensor(:,:,img_no) = fx_mat;
    
end

[Max_FV, Mesh]= max(fx_mat_tensor,[],3);
Mesh(Max_FV==0) = 0;

binRange = 2;
outThres = 3;
OC_patch = 4;
Mesh = OutlierCorrectionEdge(OC_patch,binRange,outThres,Mesh,s(1),s(2));

% %% ********************************************************************************
% % zero delete 
[zero_y,zero_x] = find(Max_FV == 0);

for i = 1:length(zero_y) 
         
         Mesh(zero_y(i),zero_x(i)) = 1;
     
end

Mesh = round(Mesh);

%
% zero correction 

% [zero_y,zero_x] = find(Mesh == 0);
% 
% for i = 1:length(zero_y)
% 
%      NZR = find(Max_FV (zero_y(i),:) ~= 0); %% Non-zero for row
% 
%      NZC = find(Max_FV (:,zero_x(i)) ~= 0); %% Non-zero for col
% 
%      if  ~isempty(NZR) && ~isempty(NZC) 
%          % distance between each points and target point in row and col
%          % respectively
%          dist_row = (NZR - zero_x(i)).^2; 
%          dist_col = (NZC - zero_y(i)).^2; 
%         
%          % closest non-zero point in row
%          closeR_col = NZR (dist_row == min(dist_row));
%     
%          % closest non-zero point in col
%          closeC_row = NZC (dist_col == min(dist_col));
% 
%          w1 = min(dist_row)/(min(dist_row)+min(dist_col));
%          w2 = min(dist_col)/(min(dist_row)+min(dist_col));
% 
%          Mesh(zero_y(i),zero_x(i)) = sum([w1*mean(Mesh(zero_y(i),closeR_col)); w2*mean(Mesh(closeC_row,zero_x(i)))],'all');
%      else
%          
%          Mesh(zero_y(i),zero_x(i)) = 1;
%      end
%      
% end
% 
% Mesh = round(Mesh);



%
%
% Depth transfer (cm)
% im_dist indicates F
%
im_dist = 5.3;
focus_depth = 5;
Mesh = alpha_M(Mesh);
Mesh = im_dist*Mesh;

Mesh = (focus_depth*Mesh)./(Mesh-focus_depth); 
%
% binRange = 2;
% outThres = 2;
% OC_patch = 30;
% Mesh = OutlierCorrection(OC_patch,binRange,outThres,Mesh,s(1),s(2));

%********************************************************************************
% Depth Map
% % 
% Mesh = imgaussfilt(Mesh,1); 

image_width = size(IM_template, 2);
image_height = size(IM_template, 1);


figure();
imagesc(Mesh);
h = colorbar;
set(h, 'FontSize', 14); % 14 是字体大小

% 设置标题，并调整字体大小
title('Depth Map (Only Edge) (cm)', 'FontSize', 16); % 16 是标题字体大小
% title('Depth Map (cm)', 'FontSize', 16); % 16 是标题字体大小
colormap('hot')
fig = gcf;
fig.Position(3) = image_width   ;  % 宽度
fig.Position(4) = image_width;  % 高度
axis equal
axis tight;

%% Depth slice

for i = min(Mesh(:)):0.05:max(Mesh(:))
    logic_matrix = false(size(Mesh,1), size(Mesh,2));
    indices = find(Mesh>i-1 & Mesh<i+1);

    logic_matrix(indices) = true;

    figure(6)

%     subplot(6, 3, (i-min(Mesh(:)))+1);
    imshow(logic_matrix)
    title(sprintf('Depth Slice (Depth around %.4f)',i));
%     title(sprintf(' %.4f',i));
end





%% parallel FV
% fx_mat_tensor initialization remains the same

% Start parallel pool (if not already created)

figure(2)
imshow(IM_template)
% figure(3)
IM_e = edge(imgaussfilt(IM_template,0.01),'canny');
% IM_e = edge(IM_template,'canny');
% imshow(IM_e);

s = size(IM_template);
x = s(2);
y = s(1);
n=x; %the photo will be divided into nxn parts
%generate ROI/AOI (region/area of interest) matrix

fx_mat = zeros(1,n*n);
AOI = fROI(n,y,x); % region of interest
laplacianFilter = [0 -1 0; -1 4 -1; 0 -1 0];
hf_ratio=1;

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');
end

numImages = size(IM_mat, 3);
parfor img_no = 1:numImages
    IM = IM_mat(:,:,img_no);
    IM = imgaussfilt(IM, 1);
    IM_filtered = imfilter(IM, laplacianFilter);
    IM = IM+hf_ratio*IM_filtered;


    fx_mat = zeros(1, n * n);

    for k = 1:n * n
        [pos_row, pos_col] = find(imcrop(IM_e, AOI(k,:)) == 1);

        if ~isempty(pos_row)
            fx = 10 * fmeasure(IM, 'LAPD', AOI(k,:));
        else
            fx = 0;
        end

        fx_mat(k) = fx;
    end

    fx_mat = reshape(fx_mat, [n, n])';
    fx_mat_tensor(:,:,img_no) = fx_mat;

    % Display progress
    fprintf('Processing image %d/%d\n', img_no, numImages);
end

fprintf('Processing complete.\n');

% Close the parallel pool when done
delete(gcp('nocreate'));

[Max_FV, Mesh]= max(fx_mat_tensor,[],3);
Mesh(Max_FV==0) = 0;

%%

% 将 IM_e 与 Mesh 叠加显示，并在 IM_e 上应用伪彩色映射
figure;
overlayedImage = imfuse(Mesh, IM_e, 'blend');
imshow(overlayedImage);

% 设置颜色映射
colormap('hot');

% 显示 colorbar，并设置字体大小
h = colorbar;
set(h, 'FontSize', 14); % 14 是字体大小

% 设置标题，并调整字体大小
title('Depth Map with Edge (cm)', 'FontSize', 16); % 16 是标题字体大小