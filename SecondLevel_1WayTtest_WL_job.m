% Created by David Garnica, david.garnica@med.uni-goettingen.de
% 2022, Universitätsmedizin Göttingen, Neurology Department

% 2nd Level analysis (1 SAMPLE T-test) for WL task (2 CONTRASTS: Animals>Fractals & Animals<Fractals)
clear all; clc;
analysis_dir = 'directory';
patients_dir = 'folder_patients';
controls_dir = 'folder_controls';
func_specification = '/func/1stLevel_wl/optcens/03_both/';

patients_list = {'pseudonyms'};
controls_list = {'pseudonyms'};

contrasts = {'Animals<Fractals', 'Animals>Fractals'};
contrast_file = 'con_0001.nii,1';
% NOTE: AS I ASSIGNED 2 CONTRASTS AND THEIR 2 WEIGHTINGS (DOWN), THE con_0002.nii FILES ARE NOT NEEDED!

% STOP HERE AND SELECT THE GROUP
fprintf('Please, first select a GROUP! \n');
prompt = {'1 Patients, 2 Controls, 3 Both groups)'};
dlgtitle = 'Select a group to analyze';
dims = [1 45];
answer = inputdlg(prompt, dlgtitle, dims);
[group] = deal(answer{:});
group = str2num(group);

for p = 1:length(patients_list)
    patients{p} = fullfile(analysis_dir, patients_dir, patients_list{p}, func_specification, contrast_file);
end
patients = patients';

for c = 1:length(controls_list)
    controls{c} = fullfile(analysis_dir, controls_dir, controls_list{c}, func_specification, contrast_file);
end
controls = controls';
    
if group == 1
    scans = patients;
elseif group == 2
    scans = controls;
elseif group == 3
    scans = vertcat(patients, controls);
else
    error('Wrong input. Please repeat')
end

%% Deal with covariates
thr = 0.05; % ***THRESHOLD for p-values, to include covariates (instead of 0.05)
% Load covariates - Age, sex (zeros for females, ones for males), and handedness (zeros for left-handers, ones for right-handers)
load([analysis_dir,'SecondLevel_WL','/','covariates.mat']);
% Check for significant differences between groups in covariates
% first, distinguish patients from controls
for i = 1:size(covariates,1)
    if contains(covariates.Pseudonym(i), 'P.')
        log_subjects(i) = 0;
    elseif contains(covariates.Pseudonym(i), 'C.')
        log_subjects(i) = 1;
    end
end
log_subjects = log_subjects'; % zeros for patients, ones for controls
Age_patients = covariates.Age(log_subjects == 0); Sex_patients = covariates.Sex(log_subjects == 0); Hand_patients = covariates.Handedness(log_subjects == 0);
Age_controls = covariates.Age(log_subjects == 1); Sex_controls = covariates.Sex(log_subjects == 1); Hand_controls = covariates.Handedness(log_subjects == 1);
% Compare, using Wilcoxon rank sum test (equivalent to Mann-Whitney U test)
[p_age,h_age,stats_age] = ranksum(Age_patients, Age_controls);
[p_sex,h_sex,stats_sex] = ranksum(Sex_patients, Sex_controls);
[p_hand,h_hand,stats_hand] = ranksum(Hand_patients, Hand_controls);

%% Run SPM Batch
clear matlabbatch;
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch{1}.spm.stats.factorial_design.dir = {'directory'};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;

% Add covariates to the model if there are significant differences between groups
if p_age < thr && p_sex < thr && p_hand < thr                               % if all are sign. different
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [covariates.Age]; 
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1; 
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = [covariates.Sex];
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Sex';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(3).c = [covariates.Handedness];
    matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'Handedness';
    matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
elseif p_age < thr && p_sex < thr                                           % if age and sex are sign. different
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [covariates.Age]; 
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1; 
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = [covariates.Sex];
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Sex';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
elseif p_age < thr && p_hand < thr                                          % if age and handedness are sign. different
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [covariates.Age]; 
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1; 
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = [covariates.Handedness];
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Handedness';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
elseif p_sex < thr && p_hand < thr                                          % if sex and handedness are sign. different
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [covariates.Sex];
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Sex';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = [covariates.Handedness];
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Handedness';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
elseif p_age < thr                                                          % if only age is sign. different
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [covariates.Age]; 
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1; 
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
elseif p_sex < thr                                                          % if only sex is sign. different
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [covariates.Sex];
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Sex';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
elseif p_hand < thr                                                         % if only handedness is sign. different
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [covariates.Handedness];
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Handedness';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
else                                                                        % if any is sign. different 
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end

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
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = contrasts{2};
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec(1).titlestr = '';
matlabbatch{4}.spm.stats.results.conspec(1).contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec(1).threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec(1).thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec(1).extent = 0;
matlabbatch{4}.spm.stats.results.conspec(1).conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec(1).mask.none = 1;
matlabbatch{4}.spm.stats.results.conspec(2).titlestr = '';
matlabbatch{4}.spm.stats.results.conspec(2).contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec(2).threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec(2).thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec(2).extent = 0;
matlabbatch{4}.spm.stats.results.conspec(2).conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec(2).mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.ps = true;

spm_jobman('run', matlabbatch);