% Created by David Garnica, david.garnica@med.uni-goettingen.de
% 2024, Universitätsmedizin Göttingen, Neurology Department

% 2nd Level analysis for WL task (2 contrasts: Animals>Fractals & Animals<Fractals)
clear all;
analysis_dir = 'analysis dir';
patients_dir = 'PATIENTS';
controls_dir = 'CONTROLS';
func_specification = '/func/1stLevel_wl/optcens/03_both/';

patients_list = {'Pseudonyms'};
controls_list = {'Pseudonyms'};

% NOTE FOR ANOVA CALCULATION: 
% the con_0001.nii corresponds to the Animals>Fractals comparison!
% the con_0002.nii corresponds to the Fractals>Animals comparison!

% The patients>controls contrast designs must be 0100 for A>F and 
% controls>patients must be -1100 for A>F, when using the con_0001.nii file 
% (one contrast file is always two/sided!)

% The main effect over both groups must have the design 1000 for A>F using 
%  the con_0001.nii file and 1000 for A<F using the con_0002.nii file

% STOP HERE AND SELECT THE GROUP & CONTRAST FILE
fprintf('Please, first select a GROUP! \n');
prompt = {'Group (1 Patients/ 2 Controls/ 3 Both groups)', ...
    'Contrast (1 con_0001 Animals>Fractals / 2 con_0002 Fractals>Animals'};
dlgtitle = 'Select group and contrast to analyze';
dims = [1 72];
answer = inputdlg(prompt, dlgtitle, dims);
[group, confile] = deal(answer{:});
group = str2num(group);
confile = str2num(confile);

if group == 1
    fprintf('--> Patients selected \n');
elseif group == 2
    fprintf('--> Controls selected \n'); 
elseif group == 3
    fprintf('--> Both groups selected \n'); 
    fprintf('Please, decide if groups must be compared \n');    
    prompt = {'Compare groups? (1 Yes / 0 No, and for both groups effect)'};
    dlgtitle = 'Groups comparison';
    dims = [1 60];
    answer = inputdlg(prompt, dlgtitle, dims);
    [compare] = deal(answer{:});
    compare = str2num(compare);
    if compare == 0
        fprintf('--> Effect of both groups to be analyzed \n');
    elseif compare == 1
        fprintf('--> Effect of groups to be compared \n');
    end
end

if group == 1 && confile == 1
    contrast_file = 'con_0001.nii';
    contrast_title = 'Effect of Patients: Animals>Fractals';
elseif group == 1 && confile == 2   
    contrast_file = 'con_0002.nii';
    contrast_title = 'Effect of Patients: Fractals>Animals';
elseif group == 2 && confile == 1 
    contrast_file = 'con_0001.nii';
    contrast_title = 'Effect of Controls: Animals>Fractals';
elseif group == 2 && confile == 2  
    contrast_file = 'con_0002.nii';
    contrast_title = 'Effect of Controls: Fractals>Animals';    
elseif group == 3 && confile == 1 && compare == 0
    contrast_file = 'con_0001.nii';
    contrast_title = 'Both groups: Animals>Fractals';
elseif group == 3 && confile == 2 && compare == 0
    contrast_file = 'con_0002.nii';
    contrast_title = 'Both groups: Fractals>Animals';    
elseif group == 3 && confile == 1 && compare == 1
    contrast_file = 'con_0001.nii';
    contrast_title = {'Patients>Controls: Animals>Fractals', 'Controls>Patients: Animals>Fractals'};
elseif group == 3 && confile == 2 && compare == 1
    contrast_file = 'con_0002.nii';
    contrast_title = {'Patients>Controls: Fractals>Animals', 'Controls>Patients: Fractals>Animals'};    
end

% loops to define scans filepaths of every subject
for p = 1:length(patients_list)
    patients{p} = fullfile(analysis_dir, patients_dir, patients_list{p}, func_specification, contrast_file);
end
patients = patients';

for c = 1:length(controls_list)
    controls{c} = fullfile(analysis_dir, controls_dir, controls_list{c}, func_specification, contrast_file);
end
controls = controls';

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
Vols_controls = covariates.Volumes_WL(log_subjects == 1); 
Age_patients = covariates.Age(log_subjects == 0);
Age_controls = covariates.Age(log_subjects == 1); 
Sex_patients = covariates.Sex(log_subjects == 0); 
Sex_controls = covariates.Sex(log_subjects == 1); 

if group == 1
    scans = patients;
    Vols_subjects = Vols_patients;
    Age_subjects = Age_patients * 12; % In months
    Sex_subjects = Sex_patients;
elseif group == 2
    scans = controls;
    Vols_subjects = Vols_controls;
    Age_subjects = Age_controls * 12; % In months
    Sex_subjects = Sex_controls;
elseif group == 3
    scans = vertcat(patients, controls); % patients' scans first, later controls'
    Vols_subjects = cat(1,Vols_patients,Vols_controls); % number of fmri volumes from both patients and controls
    Age_subjects = cat(1,Age_patients,Age_controls);
    Age_subjects = Age_subjects * 12; % In months
    Sex_subjects = cat(1,Sex_patients,Sex_controls);
    Groups_vector = log_subjects;
end

%% Run SPM Batch
clear matlabbatch;
spm('defaults','fmri');
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {'/home/uni10/nmri/projects/dgarnica/MRI_EEG_PROSPECTIVE/SecondLevel_ANOVA_WL'};
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell.scans = scans;
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;

if group == 1 || group == 2
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = Vols_subjects;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Volumes';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = Age_subjects;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Age';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(3).c = Sex_subjects;
    matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'Sex';
    matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
elseif group == 3
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = Groups_vector;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Patients-Controls';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
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
end

matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.review.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.review.display.matrix = 1;
matlabbatch{2}.spm.stats.review.print = 'ps';
matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{4}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

if group == 1
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = contrast_title;
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0];
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
elseif group == 2 
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = contrast_title;
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0];
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
elseif group == 3 && compare == 0
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = contrast_title;
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0]; % FOR BOTH GROUPS EFFECT
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
elseif group == 3 && confile == 1 && compare == 1
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = contrast_title{1};
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [0 1 0 0 0]; % Patients>Controls in A>F
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = contrast_title{2};
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [0 -1 0 0 0]; % Controls>Patients in A>F
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
elseif group == 3 && confile == 2 && compare == 1
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = contrast_title{1};
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0 0 0]; % Patients>Controls in F>A
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = contrast_title{2};
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [1 -1 0 0 0]; % Controls>Patients in F>A
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
end

matlabbatch{4}.spm.stats.con.delete = 0;
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec(1).titlestr = '';
matlabbatch{5}.spm.stats.results.conspec(1).contrasts = Inf;
matlabbatch{5}.spm.stats.results.conspec(1).threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec(1).thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec(1).extent = 0;
matlabbatch{5}.spm.stats.results.conspec(1).conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec(1).mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.ps = true;

spm_jobman('run', matlabbatch);