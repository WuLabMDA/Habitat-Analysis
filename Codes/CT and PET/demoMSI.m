%	
% Demo: Multiregional spatial interaction (MSI) feature extraction
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
%
%   Copyright 2024 MAminu. This code may be freely used and distributed, so long as it maintains this copyright line
% 

clc
clear

rootDir = '...\All Data';

ClusterNum = 5;
tumorClusterDir = fullfile(rootDir, 'Results', 'TumorPartitions');
SuperpixelsDir = fullfile(rootDir, 'Results', 'Superpixels');

patientList = dir(fullfile(tumorClusterDir, '*.mat'));
numPatient = length(patientList);
msiPath = fullfile(rootDir, 'Results', 'MSI features');
if ~exist(msiPath, 'dir')
    mkdir(msiPath)
end

MSIFeatures = cell(numPatient, 2);
for p = 1:numPatient
    caseName = patientList(p).name;
    [~, pName, ~] = fileparts(caseName);
    disp(['Processing ', num2str(p), '/', num2str(numPatient), ' ', pName]);
    feaPath = fullfile(SuperpixelsDir, caseName);
    tumorPartPath = fullfile(tumorClusterDir, caseName);
    load(feaPath, 'features', 'idxCL', 'imgSize', 'pixelDim');
    load(tumorPartPath);
    
    % reconstruct the tumor cube
    row = features(:,6) - min(features(:,6));
    col = features(:,7) - min(features(:,7));
    slice = features(:,8) - min(features(:,8)); 
    map = zeros(max(row)+3, max(col)+3, max(slice)+3);
    mask = zeros(max(row)+3, max(col)+3, max(slice)+3);
    
    % populate map and mask
    for i = 1:length(idxCL)
        map(row(i)+2,col(i)+2,slice(i)+2) = idxPL(i);
        mask(row(i)+2,col(i)+2,slice(i)+2) = 1;
    end
    
    % build MSI matrix
    ind = find(mask);
    [idx1, idx2, idx3] = ind2sub(size(mask), ind);
    matrixHabitat = zeros(ClusterNum + 1); 
    for pt = 1:length(ind)
        % use 8 connectivity model in 2D
        IDX_1 = map(idx1(pt),idx2(pt),idx3(pt));
        IDX_2 = map(idx1(pt)-1,idx2(pt)-1,idx3(pt));
        
        % summarize into the matrix
        matrixHabitat(IDX_1+1,IDX_2+1) = matrixHabitat(IDX_1+1,IDX_2+1)+1;

        % the habitat label from pre and mid habitat
        IDX_1 = map(idx1(pt),idx2(pt),idx3(pt));
        IDX_2 = map(idx1(pt)-1,idx2(pt),idx3(pt));
        % summarize into the matrix
        matrixHabitat(IDX_1+1,IDX_2+1) = matrixHabitat(IDX_1+1,IDX_2+1)+1;

        % the habitat label from pre and mid habitat
        IDX_1 = map(idx1(pt),idx2(pt),idx3(pt));
        IDX_2 = map(idx1(pt)-1,idx2(pt)+1,idx3(pt));
        % summarize into the matrix
        matrixHabitat(IDX_1+1,IDX_2+1) = matrixHabitat(IDX_1+1,IDX_2+1)+1;

        % the habitat label from pre and mid habitat
        IDX_1 = map(idx1(pt),idx2(pt),idx3(pt));
        IDX_2 = map(idx1(pt),idx2(pt)-1,idx3(pt));
        % summarize into the matrix
        matrixHabitat(IDX_1+1,IDX_2+1) = matrixHabitat(IDX_1+1,IDX_2+1)+1;

        % the habitat label from pre and mid habitat
        IDX_1 = map(idx1(pt),idx2(pt),idx3(pt));
        IDX_2 = map(idx1(pt),idx2(pt)+1,idx3(pt));
        % summarize into the matrix
        matrixHabitat(IDX_1+1,IDX_2+1) = matrixHabitat(IDX_1+1,IDX_2+1)+1;

        % the habitat label from pre and mid habitat
        IDX_1 = map(idx1(pt),idx2(pt),idx3(pt));
        IDX_2 = map(idx1(pt)+1,idx2(pt)-1,idx3(pt));
        % summarize into the matrix
        matrixHabitat(IDX_1+1,IDX_2+1) = matrixHabitat(IDX_1+1,IDX_2+1)+1;

        % the habitat label from pre and mid habitat
        IDX_1 = map(idx1(pt),idx2(pt),idx3(pt));
        IDX_2 = map(idx1(pt)+1,idx2(pt),idx3(pt));
        % summarize into the matrix
        matrixHabitat(IDX_1+1,IDX_2+1) = matrixHabitat(IDX_1+1,IDX_2+1)+1;

        % the habitat label from pre and mid habitat
        IDX_1 = map(idx1(pt),idx2(pt),idx3(pt));
        IDX_2 = map(idx1(pt)+1,idx2(pt)+1,idx3(pt));
        % summarize into the matrix
        matrixHabitat(IDX_1+1,IDX_2+1) = matrixHabitat(IDX_1+1,IDX_2+1)+1;
    end   
    
    matrixHabitat(1,1) = 0;    
    for k = 1:size(matrixHabitat,1)
        matrixHabitat(1,k) = matrixHabitat(k,1);
    end
    
    stats = graycoprops(matrixHabitat,{'contrast','homogeneity','correlation','energy'});
    glcmFeatures = [stats.Contrast,stats.Homogeneity,stats.Correlation,stats.Energy];
    
    matrixHabitat = matrixHabitat.*prod(pixelDim);
    x = diag(matrixHabitat);
    x = x(2:end);
    m = tril(true(size(matrixHabitat)),-1);
    y = matrixHabitat(m);
    
    sum_mth = sum(matrixHabitat(:));
    norm_matrix_habitat = matrixHabitat./sum_mth;
    xx = diag(norm_matrix_habitat);
    xx = xx(2:end);
    n = tril(true(size(norm_matrix_habitat)),-1);
    yy = norm_matrix_habitat(n);
    
    feaMSI = [glcmFeatures x' y' xx' yy'];  % first 4 are GLCM, k habitat diagonal components, off diagnonal , normalized, 
    MSIFeatures{p,1} = pName;
    MSIFeatures{p,2} = feaMSI;
    
end
save(msiPath, 'MSIFeatures');

IDS = MSIFeatures(:,1);
tbl = array2table(cell2mat(MSIFeatures(:,2)));
tbl.Properties.RowNames = IDS;

writetable(tbl,'MSIFeaturesGEMINI.xlsx','WriteRowNames',true)
