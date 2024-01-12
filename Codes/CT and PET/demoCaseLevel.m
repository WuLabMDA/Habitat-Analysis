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
%
% Copyright 2023 MAminu. This code may be freely used and distributed, so long as it maintains this copyright line
% 

clearvars
clc

rootDir = 'C:\Projects\ICON';
imgDir = 'V:\ICON\Data';

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
lungMax = lungs(2) + lungs(1)/2; 
lungMin = lungs(2) - lungs(1)/2; 

% obtain patient list 
patientList = dir(fullfile(imgDir));
numPatients = length(patientList)-2;
casesName = {patientList.name};
casesName(1:2) = [];
allSuperpixels = cell(numPatients, 1);
idxCount = zeros(numPatients, 1);

for i = 1:numPatients
    pStruct = patientList(i+2);
    pName = pStruct.name;
    disp(['Processing ', num2str(i), '/', num2str(numPatients), ' ', pName]);

    % load CT scans
    niftiFile = ls(fullfile(imgDir, pName, 'CT.nii'));
    niftiPath = fullfile(imgDir, pName, niftiFile);
    ctImage = niftiread(niftiPath);    
    ctInfo = niftiinfo(niftiPath);
    imgSize = ctInfo.ImageSize;
    pixelDim = ctInfo.PixelDimensions;
    if ctInfo.AdditiveOffset == -1024
        ctImage = ctImage - 1024;
    end

    % load PET scans
    petFile = ls(fullfile(imgDir, pName, 'PET.nii'));
    petPath = fullfile(imgDir, pName, petFile);
    petImage = niftiread(petPath);

    % load tumor mask
    segPath = fullfile(imgDir, pName, 'Seg_tumor.nii');
    segImage = niftiread(segPath);
    seg = logical(segImage);
    
    % Get the largest tumor
    CC = bwconncomp(seg);    
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    filteredSeg = false(size(seg));
    filteredSeg(CC.PixelIdxList{idx}) = true;
    
    % smooth the edges of the tumor mask
    radius = 2;
    decomposition = 0;
    se = strel('disk',radius,decomposition);
    filteredSeg = imerode(filteredSeg, se);
    maskSeg = imfill(filteredSeg,'holes');
    
    % Normalize image
    ctLung = mat2gray(ctImage, [lungMin lungMax]);
    petLung = mat2gray(petImage, [0.01 15.0]);
    
    % compute bounding box
    [boxBound] = computeBoundingBox(maskSeg);
    maskSeg = maskSeg(boxBound(1,1):boxBound(1,2),...
        boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
    
    ctLung = ctLung(boxBound(1,1):boxBound(1,2),...
        boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
    
    petLung = petLung(boxBound(1,1):boxBound(1,2),...
        boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
    
    % compute local entropy
    entLung = entropyfilt(double(ctLung));
    entLung = rescale(entLung);
    entPet = entropyfilt(double(petLung));
    entPet = rescale(entPet);
    
    % get region of interest
    ctLung_roi = ctLung(maskSeg);
    ctPet_roi = petLung(maskSeg);
    entCT_roi = entLung(maskSeg);
    entPet_roi = entPet(maskSeg);
    
    % Look into tumor region only
    ctLung = ctLung.*maskSeg;
    petLung = petLung.*maskSeg;
    entLung = entLung.*maskSeg;
    entPet = entPet.*maskSeg;
    
    % Fuse different modalities
    imgFused = imadd(imadd(ctLung, entLung), imadd(petLung, entPet));  
    tempFused = imgFused(maskSeg);
    
    % SLIC superpixel oversegmentation
    if sum(filteredSeg,'all') < 1000
        numSuperpixels = 30;
    elseif sum(filteredSeg,'all') > 1000 && sum(filteredSeg,'all') < 10000
        numSuperpixels = 50;
    elseif sum(filteredSeg,'all') > 10000 && sum(filteredSeg,'all') < 100000
        numSuperpixels = 60;
    else
        numSuperpixels = 80;
    end
    
    [L, N] = superpixels3(imgFused, numSuperpixels, 'Method', 'slic', 'Compactness', 1e-2);
    
    indices = find(maskSeg==1);
    [rows, cols, slices] = ind2sub(size(maskSeg), indices); 
    
    features = [double(ctLung_roi), double(entCT_roi), double(ctPet_roi), double(entPet_roi),...
        double(tempFused), double(rows), double(cols), double(slices)];
    
    % Get tumor superpixels labels
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
    fInput = [Output '\' pName '.png'];
    f = figure('visible', 'on');
    %     f.WindowState = 'maximized';
    
    TumorVox = sum(sum(maskSeg,1),2);
    TumorVox = TumorVox(:);
    [~,idx] = max(TumorVox);
    BW = maskSeg(:,:,idx);
    [B1,L1] = bwboundaries(BW);
    
    subplot(3,2,1)
    imshow(ctLung(:,:,idx),[])
    title('CT')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off
    
    subplot(3,2,2);
    imshow(entLung(:,:,idx),[])
    title('entropy CT')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off

    subplot(3,2,3)
    imshow(petLung(:,:,idx),[])
    title('PET')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off
    
    subplot(3,2,4);
    imshow(entPet(:,:,idx),[])
    title('entropy PET')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off
    
    imgFused = imgFused.*maskSeg;
    subplot(3,2,5);
    imshow(imgFused(:,:,idx),[])
    title('Fused image')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off

    seg = L.*maskSeg;
    subplot(3,2,6);
    imshow(label2rgb(seg(:,:,idx),'jet','k','shuffle'));
    title('Superpixels')
    hold on
    for k=1:length(B1)
        boundary = B1{k};
        plot(boundary(:,2), boundary(:,1),...
            'r-', 'LineWidth', 1)
    end
    hold off


    saveas(f,fInput)
    close(f)
    
    p_intratumoral_path = fullfile(SuperpixelsDir, pName);
    save(p_intratumoral_path, 'features', 'idxCL', 'indices', 'label', 'imgSize', 'pixelDim');
end
path = fullfile(rootDir, 'Results', 'supPixels');
allSuperpixels = cell2mat(allSuperpixels);
save(path, 'allSuperpixels', 'idxCount', 'casesName');


