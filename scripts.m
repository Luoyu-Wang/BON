%% ========================================================================
%  Brain-Organ Network (BON) Analysis Pipeline 
% ========================================================================

%% 1. CT & PET Coregistration and Bounding Box Extraction
vol_spm_coregister_estwrite('CTBodyImg', 'PETPI2620BodyImg', {'ParametersCTBodyImgSegment'});

subdir = 'PETPI2620BodyImg';
idir = 'ParametersCTBodyImgSegmentCoregisterEstwrite';
odir = 'PETPI2620BrainImg';
md(odir);

sublist = dir(fullfile(subdir, '*.nii'));
ilist = dir(fullfile(idir, '*.nii'));

for i = 1:numel(ilist)
    [data1, hdr1] = vol_spm_read(fullfile(subdir, sublist(i).name));
    [data, hdr]  = vol_spm_read(fullfile(idir, ilist(i).name));
    
    data(data~=90 & data~=91) = 0;

    % Find the largest connected component
    CC = bwconncomp(data);
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [~, idx] = max(numPixels);
    data = false(size(data));
    data(CC.PixelIdxList{idx}) = 1;

    % Get the minimum bounding box containing all true elements
    [x_indices, y_indices, z_indices] = ind2sub(size(data), find(data));
    xmin = min(x_indices)-50; xmax = max(x_indices)+50;
    ymin = min(y_indices)-50; ymax = max(y_indices)+50;
    zmin = min(z_indices)-5;  zmax = max(z_indices)+5;

    rectangular_mask = false(size(data));
    rectangular_mask(xmin:xmax, ymin:ymax, zmin:zmax) = 1;

    % Extract relevant slices based on the bounding box
    slices_to_keep = any(any(rectangular_mask, 1), 2);
    filtered_data = data1(:, :, slices_to_keep, :);

    hdr1.mat(3,4) = -100;
    vol_spm_write(filtered_data, hdr1, fullfile(odir, sublist(i).name));
end


%% 2. Organ Segmentation Extraction (CT & PET)
id = [1,2,3, 5:21, 25:64, 71:74, 77:89, 116];
ctdir = 'CTBodyImg';
idir = 'ParametersCTBodyImgSegment';
odir2 = 'CTOrgansImg';
odir3 = 'ParametersOrgansBodyImgSegment';
md(odir3); md(odir2);

ilist = dir(fullfile(idir, '*.nii'));
ctlist = dir(fullfile(ctdir, '*.nii'));

for i = 1:numel(ilist)
    [ct, hdrct] = vol_spm_read(fullfile(ctdir, ctlist(i).name));
    [data, hdr] = vol_spm_read(fullfile(idir, ilist(i).name));
    
    data2 = data;
    data2(~ismember(data, id)) = 0;
    data = data2;
    
    CC = bwconncomp(data2);
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [~, idx] = max(numPixels);
    data2 = false(size(data2));
    data2(CC.PixelIdxList{idx}) = 1;

    [x_indices, y_indices, z_indices] = ind2sub(size(data2), find(data2));
    xmin = min(x_indices); xmax = max(x_indices);
    ymin = min(y_indices); ymax = max(y_indices);
    zmin = min(z_indices); zmax = max(z_indices);

    rectangular_mask = false(size(data2));
    rectangular_mask(xmin:xmax, ymin:ymax, zmin:zmax) = 1;

    rows_to_keep = any(any(rectangular_mask, 2), 3);
    cols_to_keep = any(any(rectangular_mask, 1), 3);
    slices_to_keep = any(any(rectangular_mask, 1), 2);

    ct = ct(rows_to_keep, cols_to_keep, slices_to_keep);
    seg = data(rows_to_keep, cols_to_keep, slices_to_keep);

    % Update Header matrix transformation
    hdr.mat(1,4) = 0.5*(xmax - xmin) * abs(hdr.mat(1,1));
    hdr.mat(2,4) = 0.5*(-ymax + ymin) * abs(hdr.mat(2,2));
    hdr.mat(3,4) = 0.5*(-zmax + zmin) * abs(hdr.mat(3,3));
    
    hdrct.mat(1,4) = hdr.mat(1,4);
    hdrct.mat(2,4) = hdr.mat(2,4);
    hdrct.mat(3,4) = hdr.mat(3,4);

    vol_spm_write(ct, hdrct, fullfile(odir2, ilist(i).name));
    vol_spm_write(seg, hdr, fullfile(odir3, ilist(i).name));
end


%% 3. T1 & Brain Atlas Processing (Realign & Coregister DK14 to generate DK10)
vol_spm_realign('PETPI2620BrainImg');
% Note: DK14 is kept here because DK10 is derived from DK14 in Step 9.
vol_spm_coregister_estwrite('T1ImgON4Alignment', 'ParametersRealignMeanFun', {'ParametersDK14Structure'});


%% 4. Organ Smoothing, Scaling and Alignment
Parameters.FWHM = [8, 8, 8];
vol_fsl_smooth('PETPI2620OrgansImg', Parameters);
vol_spm_scale_pet_wholebody('PETPI2620OrgansImgSmooth8');

rdirs = {'PETPI2620OrgansImgSmooth8', 'ParametersPETProbability0', 'ParametersPETProbability1'};
vol_ants_realign_nolinear('PETPI2620OrgansImg', rdirs);

% Clean up intermediate files
rmdir('ParametersPETAntsRealignNolinear*', 's');
rmdir('ParametersPETProbability*', 's');
rmdir('PETPI2620OrgansImgR0*', 's');
rmdir('PETPI2620OrgansImgSmooth*', 's');
rmdir('PETPI2620OrgansImgV*', 's');


%% 5. Morphological Erosion on Organs
idir = 'ParametersOrgansBodyImgSegment';
odir = 'ParametersOrgansBodyImgSegmentImerode';
md(odir);
ilist = dir(fullfile(idir, '*.nii'));

for i = 1:numel(ilist)
    [data, hdr] = vol_spm_read(fullfile(idir, ilist(i).name));
    idj = unique(data);
    idj(1) = []; % Remove background (0)
    mask0 = zeros([size(data), length(idj)]);
    
    for j = 1:numel(idj)
        mask = (data == idj(j));
        se = strel('sphere', 1);
        erodedCurrentMask = imerode(mask, se);
        erodedMask = zeros(size(mask));
        erodedMask(erodedCurrentMask) = idj(j);
        mask0(:,:,:,j) = erodedMask;
    end
    mask1 = sum(mask0, 4);
    vol_spm_write(mask1, hdr, fullfile(odir, ilist(i).name));
end


%% 6. CT/PET Preparation & Coregistration
vol_spm_scale_ct_wholebody('CTBodyImgCoregisterEstwriteOrgans');
gunz('PETPI2620OrgansImgR');
vol_spm_mean('PETPI2620OrgansImgR');
vol_spm_coregister_estwrite('CTOrgansImg', 'PETPI2620OrgansImgR3D', {'ParametersOrgansBodyImgSegmentImerode'});


%% 7. Aorta Extraction and Iterative Erosion
% Separate aorta and other organs
idir = 'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwrite';
odir_aorta = 'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwriteAorta';
odir_77 = 'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwrite77';
md(odir_aorta); md(odir_77);

ilist = dir(fullfile(idir, '*.nii'));
for i = 1:numel(ilist)
    [data, hdr] = vol_spm_read(fullfile(idir, ilist(i).name));
    data_aorta = data; data_aorta(data_aorta ~= 52) = 0;
    data_77 = data; data_77(data_77 == 52) = 0;
    
    vol_spm_write(data_aorta, hdr, fullfile(odir_aorta, ilist(i).name));
    vol_spm_write(data_77, hdr, fullfile(odir_77, ilist(i).name));
end

% Merge parcels
vol_spm_parcels_merge(odir_77, 'merge.csv', 'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwrite9');
vol_spm_parcels_merge(odir_77, 'merge10.csv', 'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwrite10');

% Multi-level erosion for aorta (Levels 3, 2, 1)
erosion_levels = [3, 2, 1];
for lvl = erosion_levels
    odir_erode = ['ParametersOrgansBodyImgSegmentImerodeCoregisterEstwriteAortaImerode', num2str(lvl)];
    md(odir_erode);
    for i = 1:numel(ilist)
        [data, hdr] = vol_spm_read(fullfile(odir_aorta, ilist(i).name));
        idj = unique(data); idj(1) = [];
        mask0 = zeros([size(data), length(idj)]);
        for j = 1:numel(idj)
            mask = (data == idj(j));
            se = strel('sphere', lvl); % Dynamically apply erosion radius
            erodedCurrentMask = imerode(mask, se);
            erodedMask = zeros(size(mask));
            erodedMask(erodedCurrentMask) = idj(j);
            mask0(:,:,:,j) = erodedMask;
        end
        vol_spm_write(sum(mask0, 4), hdr, fullfile(odir_erode, ilist(i).name));
    end
end


%% 8. SUVR Calculation & Regressing Nuisance Covariates
vol_spm_pet_suvr('PETPI2620OrgansImgR', 'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwriteAorta');
vol_spm_pet_suvr2('PETPI2620BrainImgR', 'PETPI2620OrgansImgR', 'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwriteAorta');

vol_spm_dvars('PETPI2620OrgansImg', 'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwrite77', 'ParametersRealignDVARS');

% Combine head motion RP and DVARS
movefile('ParametersRealignRP', 'ParametersRealignRP6');
ilist1 = dir(fullfile('ParametersRealignRP6', '*.txt'));
ilist2 = dir(fullfile('ParametersRealignDVARS', '*.txt'));
md('ParametersRealignRP');
for i = 1:numel(ilist1)
    a = load(fullfile('ParametersRealignRP6', ilist1(i).name));
    b = load(fullfile('ParametersRealignDVARS', ilist2(i).name));
    c = [a, b];
    save(fullfile('ParametersRealignRP', ilist1(i).name), 'c', '-ascii');
end

% Regression parameters setup
Parameters.RegressHeadMotion = 3;
% Parameters.RegressWholeBrain = 0;
% Parameters.RegressCSF = 0;
% Parameters.RegressWhiteMatter = 0;
Parameters.AddMeanBack = 1;

vol_spm_regress('PETPI2620BrainImgRRate', Parameters);
Parameters.RegressHeadMotion = 0; % Organs regression params
vol_spm_regress('PETPI2620OrgansImgRRate', Parameters);


%% 9. Generate DK10 Atlas & Extract ROI Signals
% Extract Organ Signals (9 Organs)
vol_extract_signals_privatespace('PETPI2620OrgansImgRRateC0', {'ParametersOrgansBodyImgSegmentImerodeCoregisterEstwrite9'}, 1);

% Generate DK10 from DK14 (Removing nodes 11, 12, 13, 14)
md('ParametersDK10StructureCoregisterEstwrite');
ilist = dir(fullfile('ParametersDK14StructureCoregisterEstwrite', '*.nii'));
for i = 1:numel(ilist)
    [data, hdr] = vol_spm_read(fullfile('ParametersDK14StructureCoregisterEstwrite', ilist(i).name));
    data(ismember(data, [11, 12, 13, 14])) = 0;
    vol_spm_write(data, hdr, fullfile('ParametersDK10StructureCoregisterEstwrite', ilist(i).name));
end

% Extract Brain Signals (DK10)
vol_extract_signals_privatespace('PETPI2620BrainImgRRateC0', {'ParametersDK10StructureCoregisterEstwrite'}, 1);


%% ========================================================================
%  10. Brain-Organ Network Combine, Plot & NBS (DK10 + Organs 9)
% ========================================================================
labels1 = {'Frontal Lobe L', 'Frontal Lobe R', 'Parietal Lobe L', 'Parietal Lobe R', 'Occipital Lobe L', 'Occipital Lobe R', 'Temporal Lobe L', 'Temporal Lobe R','Insular Lobe L', 'Insular Lobe R'};
labels2 = {'Kidney', 'Liver', 'Lung', 'Colon', 'Spleen', 'skeleton', 'Heart', 'Spinal Cord', 'Gluteus Muscle'};
labels3 = [labels1, labels2];

combine_and_plot_networks('TimeSignalParametersDK10StructureCoregisterEstwrite', 'TimeSignalParametersOrgansBodyImgSegmentImerodeCoregisterEstwrite9', 'TimeSignalDK10Organs9Merge', labels3, 1:23, 24:51);


%% ========================================================================
%  Helper Functions (Local functions to streamline repetitive logic)
% ========================================================================

function combine_and_plot_networks(idir1, idir2, out_dir, labels, id_hc, id_ad, do_plot)
    if nargin < 7; do_plot = true; end
    ilist1 = dir(fullfile(idir1, 'ROISignals*.txt'));
    ilist2 = dir(fullfile(idir2, 'ROISignals*.txt'));
    
    out_dir_fz = [out_dir, 'FZ'];
    md(out_dir); md(out_dir_fz);
    
    for i = 1:numel(ilist1)
        sig = [load(fullfile(idir1, ilist1(i).name)), load(fullfile(idir2, ilist2(i).name))];
        ROICorrelation = corrcoef(sig);
        save(fullfile(out_dir, ['ROICorrelation_', ilist1(i).name(1:end-4), '.mat']), 'ROICorrelation');
        
        ROICorrelation_FisherZ = 0.5 * log((1 + ROICorrelation)./(1 - ROICorrelation));
        save(fullfile(out_dir_fz, ['ROICorrelation_FisherZ_', ilist1(i).name(1:end-4), '.mat']), 'ROICorrelation_FisherZ');
    end
    
    if do_plot
        colors = jet(length(labels));
        plot_mcorrmap(out_dir, [out_dir, '_SUVRC_HC'], id_hc, labels, 12, colors);
        plot_mcorrmap(out_dir, [out_dir, '_SUVRC_AD'], id_ad, labels, 12, colors);
        mat_plot_network_minus([out_dir, '_SUVRC'], [out_dir, '_SUVRC_AD_minus_HC'], id_ad, id_hc, labels, 12, colors);
    end
end

