%	
% Demo: Habitat Analysis (Population Level)
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
%               Pingjun Chen (pingjunchen AT ieee.org)
%               Jia Wu (jwu11 AT mdanderson.org)

clear
clc

rootDir = '...\DI for immunotherapy response';

SubregionsDir = fullfile(rootDir, 'Results', 'TumorPartitions');
if ~exist(SubregionsDir, 'dir')
    mkdir(SubregionsDir)
end

HabitatPath = fullfile(rootDir, 'Results', 'Habitat maps');
if ~exist(HabitatPath, 'dir')
    mkdir(HabitatPath);
end 

SuperpixelsDir = fullfile(rootDir, 'Results','Superpixels');
supPath = fullfile(rootDir, 'Results', 'supPixels.mat');

caseList = dir(fullfile(SuperpixelsDir, '*.mat'));
numCases = length(caseList);

load(supPath);
allSuperpixels = zscore(allSuperpixels);

rng('default') 
tsneFea = tsne(allSuperpixels,'Algorithm','barneshut'); 
subplot(1, 2, 1);
gscatter(tsneFea(:,1), tsneFea(:,2));

numCluster = 5; 
numNeighbor = size(allSuperpixels, 1);
[idxPL, ~, ~] = spectralcluster(allSuperpixels, numCluster,'NumNeighbors',numNeighbor,'ClusterMethod','kmedoids');
% idxPL = kmeans(allSuperpixels,numCluster,'MaxIter',100,'Replicates',10);
subplot(1, 2, 2);
gscatter(tsneFea(:,1), tsneFea(:,2), idxPL);

gnd = mat2cell(idxPL, idxCount);

% plot and save clustering results
for p = 1:numCases
    caseName = caseList(p).name;
    featuresPath = fullfile(SuperpixelsDir, string(caseName));
    disp(['Processing ', num2str(p), ' of ',num2str(numCases)]);
    load(featuresPath, 'features', 'idxCL', 'indices', 'label', 'imgSize', 'pixelDim');
    
    % update superpixels label
    superpixelsIdx = gnd{p};
    pixelsIdx = idxCL;
    idxPL = idxCL;
    for i = 1:length(label)
        habitatId = superpixelsIdx(i);
        idxPL(pixelsIdx==label(i)) = habitatId;
    end
        
    % plot habitat regions
    rows = features(:,6) - min(features(:,6));
    cols = features(:,7) - min(features(:,7));
    slices = features(:,8) - min(features(:,8)); 
    reconCT = zeros(max(rows)+3, max(cols)+3, max(slices)+3);
    reconFused = reconCT;
    seg = reconCT;
    map = reconCT;
    for j = 1:length(idxCL)
        reconCT(rows(j)+2,cols(j)+2,slices(j)+2) = features(j,1);
        reconFused(rows(j)+2,cols(j)+2,slices(j)+2) = features(j,5);
        seg(rows(j)+2,cols(j)+2,slices(j)+2) = idxCL(j);
        map(rows(j)+2,cols(j)+2,slices(j)+2) = idxPL(j);
    end
    
    habitatname = [HabitatPath '\' char(caseName) '_Habitats.tiff'];
    f = figure('visible', 'on');
%     f.WindowState = 'maximized';
    
    TumorVox = sum(sum(logical(map),1),2);
    TumorVox = TumorVox(:);
    [~,idx] = max(TumorVox);
    
    subplot(2,2,1)
    imshow(reconCT(:,:,idx),[])
    title('Original CT')
    
    subplot(2,2,2)
    imshow(reconFused(:,:,idx),[])
    title('Fused image')

    ax = subplot(2,2,3);
    imshow(label2rgb(seg(:,:,idx),'jet','k','shuffle'));
    title('Superpixels')
    
    bx = subplot(2,2,4);
    imshow(map(:,:,idx),[])
    title('Habitats')
    colormap(bx,jet)
    
    saveas(f,habitatname)
    close(f)    
   
    subregionsPath = fullfile(SubregionsDir, char(caseName));
    save(subregionsPath, 'features', 'idxCL', 'indices', 'label', 'idxPL',...
        'imgSize', 'pixelDim');
end

