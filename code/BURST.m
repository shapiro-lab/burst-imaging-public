function [BURST_im, Bmode_im, im_timeseries, P] = BURST(directory_name, use_rf, zROI, iV)
%BURST Apply Burst Ultrasound Reconstruction with Signal Templates
%(BURST) to the burst timeseries data stored in directory_name
% INPUTS
%   directory_name -- path to the folder containing the timeseries data
%   use_rf -- whether to reconstruct from the raw RF data (default: false).
%   Use this when you want full transparency in the processing pipeline.
%   For most applications, the images from Verasonics's proprietary
%   reconstruction algorithm look better.
% OUTPUTS
%   BURST_im -- N x M array containing the linear processed BURST image
%   Bmode_im -- N x M array containing the linear B-mode image of the 1st
%   collapse frame
%   x -- N x 1 vector specifying the horizontal axis
%   z -- M x 1 vector specifying the vertical axis
%
% Example Usage
%{
directory_name = '_myData/Sensitivity/20190602_MultiRegimeNisslePhantoms/TMM_HI/2019_06_02@17_01_51_A2CNissleOD1e-1_64ap_L11_11.4MHz_50.0V_SyncBurst/';
[BURST_im, Bmode_im, x, z] = BURST(directory_name);

figure, colormap hot
subplot(121)
imagesc(x,z,20*log10(Bmode_im), [40 inf]), colorbar, axis tight
title('BURST'), xlabel('Lateral position (mm)'), ylabel('Depth (mm)')
subplot(122)
imagesc(x,z,20*log10(BURST_im), [40 inf]), colorbar, axis tight
title('B-mode'), xlabel('Lateral position (mm)'), ylabel('Depth (mm)')
%}

%% Ensure directory_name has trailing (back)slash
archstr = computer('arch');
if strcmpi(archstr,'pcwin64')
    dir_sep = '\';
else
    dir_sep = '/';
end
if ~strcmpi(directory_name(end), dir_sep)
    directory_name = [directory_name dir_sep];
end

%% Load the data from files
if ~exist('use_rf','var')
    use_rf = false;
end
if ~use_rf
    files = dir([directory_name 'image*.mat']);
    n_frames = length(files);
    Im = cell(1,n_frames);
    for i = 1:n_frames
        filename = [directory_name files(i).name];
        load(filename,'P','RData','x','z')
        Im{i} = RData;
    end
else
    files = dir([directory_name 'RF_blocks*.mat']);
    filename = [directory_name files.name];
    load(filename,'P','RFSeries','x','z')
    n_frames = length(RFSeries);
    Im = cell(1,n_frames);
    for i = 1:n_frames
        Im{i} = abs(hilbert(RFSeries{i}));
    end
end


if ~exist('zROI','var')
    zROI = 1:size(Im{1},1);
end

if ~exist('iV','var')
    iV = P.iV;
end
%% Create the template vectors and unmixing matrix
Na = length(Im);
colFrame = P.numPreColFrames + 1;
for k = 1:Na
    Im{k} = Im{k}(zROI,:);
end

wellTemplate = [zeros(colFrame-1,1); 1; zeros(Na-colFrame,1)];
specTemplate = [zeros(colFrame-1,1); ones(Na-colFrame+1,1)];
noiseTemplate = ones(Na,1);

imMat = cat(3,Im{:});
imMat = imMat(:,:,iV);

wellVec = (wellTemplate(iV));
specVec = (specTemplate(iV));
noiseVec = (noiseTemplate(iV));
V = [wellVec specVec noiseVec];
Unmixer = pinv(V'*V)*V';

%% Execute BURST processing
ImW = nan([size(Im{1}),size(V,2)]);
ImD = nan(size(Im{1}));
for i = 1:size(Im{1},1)
    for j = 1:size(Im{1},2)
        imVec = squeeze(imMat(i,j,:));
        
        % Linear unmixing
        weights = Unmixer * imVec;
        weights(weights<0) = 0;
        ImW(i,j,:) = weights;
        
        % Max difference
        ImD(i,j) = max(imVec(colFrame:end)) - min(imVec(colFrame:end));
    end
end

BURST_im = ImW(:,:,1);
Bmode_im = Im{colFrame};
postcol_im = Im{colFrame+1};
im_timeseries = Im;