% 2nd Level analysis for WL task (2 CONTRASTS: Animals>Fractals & Animals<Fractals)
clear all; clc;
analysis_dir = 'analysis dir';
patients_dir = 'PATIENTS';
controls_dir = 'CONTROLS';
func_specification = '/func/1stLevel_wl/optcens/03_both/';

patients_list = {'Pseudonyms'};
controls_list = {'Pseudonyms'};

contrasts = {'Patients>Controls: Animals>Fractals', 'Patients>Controls: Fractals>Animals'};
contrast_file = 'con_0001.nii,1';
% NOTE FOR 1sample T-test CALCULATION: 
% the con_0001.nii corresponds to the Animals>Fractals comparison!
% the con_0002.nii corresponds to the Animals<Fractals comparison!
% The contrast vector must be 0100 for Patients>Controls for A>F and
% -1100 for Patients>Controls for F>A, when using
% the con_0001.nii file (one contrast file is always two-sided!)

for p = 1:length(patients_list)
    patients{p} = fullfile(analysis_dir, patients_dir, patients_list{p}, func_specification, contrast_file);
end
patients = patients';

for c = 1:length(controls_list)
    controls{c} = fullfile(analysis_dir, controls_dir, controls_list{c}, func_specification, contrast_file);
end
controls = controls';
    
scans = vertcat(patients, controls);

%% Deal with covariates
% Load covariates - Age, sex (zeros for females, ones for males), 
% handedness (zeros for left-handers, ones for right-handers), and
% the number of volumes as the WL paradigm was self-paced
load([analysis_dir,'SecondLevel_WL','/','covariates.mat']);

% first, distinguish patients from controls
for i = 1:size(covariates,1)
    if contains(covariates.Pseudonym(i), 'P.')
        log_subjects(i) = 0;
    elseif contains(covariates.Pseudonym(i), 'C.')
        log_subjects(i) = 1;
    end
end
log_subjects = log_subjects'; % zeros for patients, ones for controls

Vols_patients = covariates.Volumes_WL(log_subjects == 0); 
Age_patients = covariates.Age(log_subjects == 0); 
Sex_patients = covariates.Sex(log_subjects == 0); 
Hand_patients = covariates.Handedness(log_subjects == 0);

Vols_controls = covariates.Volumes_WL(log_subjects == 1); 
Age_controls = covariates.Age(log_subjects == 1); 
Sex_controls = covariates.Sex(log_subjects == 1); 
Hand_controls = covariates.Handedness(log_subjects == 1);

Vols_subjects = cat(1,Vols_patients,Vols_controls);
Age_subjects = cat(1,Age_patients,Age_controls);
Age_subjects = Age_subjects * 12; % In months
Sex_subjects = cat(1,Sex_patients,Sex_controls);

design1 = log_subjects;
previousOnes = design1 == 1;
design2(design1 == 0) = 1;
design2(previousOnes) = 0;
design2 = design2';

%% Run SPM batch
clear matlabbatch;
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch{1}.spm.stats.factorial_design.dir = {'/home/uni10/nmri/projects/dgarnica/MRI_EEG_PROSPECTIVE/SecondLevel_1sTtest_WL'};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.cov(1).c = design1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Design1';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;

% matlabbatch{1}.spm.stats.factorial_design.cov(2).c = design2;
% matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Design2';
% matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
% matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(2).c = Vols_subjects;
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Volumes';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(3).c = Age_subjects;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'Age';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(4).c = Sex_subjects;
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'Sex';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = contrasts{1};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0,1,0,0,0];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = contrasts{2};
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1,1,0,0,0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec(1).titlestr = '';
matlabbatch{4}.spm.stats.results.conspec(1).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(1).threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec(1).thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec(1).extent = 0;
matlabbatch{4}.spm.stats.results.conspec(1).conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec(1).mask.none = 1;
% matlabbatch{4}.spm.stats.results.conspec(2).titlestr = '';
% matlabbatch{4}.spm.stats.results.conspec(2).contrasts = 1;
% matlabbatch{4}.spm.stats.results.conspec(2).threshdesc = 'FWE';
% matlabbatch{4}.spm.stats.results.conspec(2).thresh = 0.05;
% matlabbatch{4}.spm.stats.results.conspec(2).extent = 0;
% matlabbatch{4}.spm.stats.results.conspec(2).conjunction = 1;
% matlabbatch{4}.spm.stats.results.conspec(2).mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.ps = true;

spm_jobman('run', matlabbatch);