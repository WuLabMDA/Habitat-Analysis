%	
% Demo: Habitat Analysis (Case Level)
% 
% Reference:
%
%   1. Wu, Jia, Guohong Cao, Xiaoli Sun, Juheon Lee, Daniel L. Rubin, Sandy Napel, Allison W. Kurian, Bruce L. Daniel, and Ruijiang Li. 
%   "Intratumoral spatial heterogeneity at perfusion MR imaging predicts recurrence-free survival
%   in locally advanced breast cancer treated with neoadjuvant chemotherapy." 
%   Radiology 288, no. 1 (2018): 26-35.
%
%   2. Aminu, Muhammad, Divya Yadav, Lingzhi Hong, Elliana Young, Paul Edelkamp Jr, Maliazurina Saad, Morteza Salehjahromi et al. 
%   "Habitat Imaging Biomarkers for Diagnosis and Prognosis in Cancer Patients Infected with COVID-19." Cancers 15, no. 1 (2022): 275.
%
%   version 1.0 --May/2023
%
%   Written by: Muhammad Aminu (muhammadaminu47 AT gmail.com)
%               Jia Wu (jwu11 AT mdanderson.org)

clear
clc

imgDir = 'V:\Amgad\PD-L1 Project\Aminu_Code\Data';
datapath = dir(fullfile(imgDir,'CT'));
segpath = dir(fullfile(imgDir,'tumor_mask'));
datapath(1:2) = [];
segpath(1:2) = [];

rootDir = 'V:\Amgad\PD-L1 Project\Aminu_Code';
Output=[rootDir '\Results\InputImages'];
if ~exist(Output, 'dir')
    mkdir(Output);
end 

SuperpixelsDir = fullfile(rootDir, 'Results', 'Superpixels');
if ~exist(SuperpixelsDir, 'dir')
    mkdir(SuperpixelsDir)
end 

% Windowing (contrast enhancement)
lungs = [1500 -600]; % [WindowLevel WindowWidth]
mediastinum = [350 50];
lungMax = lungs(2) + lungs(1)/2; 
lungMin = lungs(2) - lungs(1)/2; 
medMax = mediastinum(2) + mediastinum(1)/2; 
medMin = mediastinum(2) - mediastinum(1)/2; 

% obtain patient list 
numCases = numel(segpath);
allSuperpixels = cell(numel(datapath), 1);
idxCount = zeros(numel(datapath), 1);

for i = 1:numel(datapath)
    pipath = dir(fullfile(datapath(i).folder, datapath(i).name));
%     pipath(1:2) = [];
    disp(['Processing: ',datapath(i).name]);

    % load CT scans
    ctImg = niftiread(fullfile(pipath.folder, pipath.name));
    ctInfo = niftiinfo(fullfile(pipath.folder, pipath.name));
    imgSize = ctInfo.ImageSize;
    pixelDim = ctInfo.PixelDimensions;
    
    if ctInfo.AdditiveOffset == -1024
        ctImg = ctImg - 1024;
    end

    % load tumor mask
    seg = logical(niftiread(fullfile(segpath(i).folder, segpath(i).name)));

    % smooth the edges of the tumor mask
    radius = 2;
    decomposition = 0;
    se = strel('disk',radius,decomposition);
    maskSeg = imerode(seg, se);

    % Normalize image
    ctLung = mat2gray(ctImg, [lungMin lungMax]);
    ctMed = mat2gray(ctImg, [medMin medMax]);

    % compute bounding box
    [boxBound] = computeBoundingBox(maskSeg);
    maskSeg = maskSeg(boxBound(1,1):boxBound(1,2),...
        boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

    ctLung = ctLung(boxBound(1,1):boxBound(1,2),...
        boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

    ctMed = ctMed(boxBound(1,1):boxBound(1,2),...
        boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

    % compute local entropy
    entLung = entropyfilt(double(ctLung));
    entLung = rescale(entLung);
    entMed = entropyfilt(double(ctMed));
    entMed = rescale(entMed);

    % get region of interest
    ctLung_roi = ctLung(maskSeg);
    ctMed_roi = ctMed(maskSeg);
    entCT_roi = entLung(maskSeg);
    entMed_roi = entMed(maskSeg);

    % Look into tumor region only
    ctLung = ctLung.*maskSeg;
    ctMed = ctMed.*maskSeg;
    entLung = entLung.*maskSeg;
    entMed = entMed.*maskSeg;

    % Fuse different modalities
    imgFused = imadd(imadd(ctLung, entLung), imadd(ctMed, entMed));  % fuse multiple 3D images, imlincomb, imadd, imcompliment, ...
    tempFused = imgFused(maskSeg);

    % SLIC superpixel ooversegmentation
    if sum(maskSeg,'all') < 1000
        numSuperpixels = 30;
    elseif sum(maskSeg,'all') > 1000 && sum(maskSeg,'all') < 10000
        numSuperpixels = 50;
    elseif sum(maskSeg,'all') > 10000 && sum(maskSeg,'all') < 100000
        numSuperpixels = 60;
    else
        numSuperpixels = 80;
    end

    [L, N] = superpixels3(imgFused, numSuperpixels, 'Method', 'slic', 'Compactness', 1e-2);

    indices = find(maskSeg==1);
    [rows, cols, slices] = ind2sub(size(maskSeg), indices); % coordinates

    features = [double(ctLung_roi), double(entCT_roi), double(ctMed_roi), double(entMed_roi),...
        double(tempFused), double(rows), double(cols), double(slices)];

    idxCL = L(maskSeg);

    superpixel = zeros(length(unique(idxCL)), 40);
    label = unique(idxCL);
    for kk = 1:length(label)
        ind = find(idxCL==label(kk));
        clusterData = features(ind,1:end-4);
        Sfea = histFeas(clusterData);
        superpixel(kk,:) = Sfea;
    end
    allSuperpixels{i} = superpixel;
    idxCount(i) = length(unique(idxCL));

    % save outputs
    name = string(extractBetween(segpath(i).name,"RTS_","_P.nii.gz"));
    filename_input = [Output '\' convertStringsToChars(name) '.png'];
    f = figure('visible', 'on');
    %     f.WindowState = 'maximized';

    TumorVox = sum(maskSeg,1);
    TumorVox = sum(TumorVox,2);
    TumorVox = TumorVox(:);
    [~,idx] = max(TumorVox);
    BW = maskSeg(:,:,idx);
    [B1,L1] = bwboundaries(BW);

    subplot(2,2,1)
    imshow(ctLung(:,:,idx),[])
    title('CT Lung')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off

    subplot(2,2,2)
    imshow(ctMed(:,:,idx),[])
    title('CT Med')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off

    subplot(2,2,3);
    imshow(entLung(:,:,idx),[])
    title('entropy CT Lung')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off

    subplot(2,2,4);
    imshow(entMed(:,:,idx),[])
    title('entropy CT Med')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off

    saveas(f,filename_input)
    close(f)

    p_intratumoral_path = fullfile(SuperpixelsDir, name);
    save(p_intratumoral_path, 'features', 'idxCL', 'indices', 'label', 'imgSize', 'pixelDim');
end
path = fullfile(rootDir, 'Results', 'supPixels');
allSuperpixels = cell2mat(allSuperpixels);
save(path, 'allSuperpixels', 'idxCount');
