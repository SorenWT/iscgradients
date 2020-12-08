
basedir = '~/Desktop/masters/NEUR608/iscproject';
cd(basedir)
mkdir glm

worddata = h5read('WordNetFeatures.hdf5','/MOVIE1_CC1');

for i = 1:size(worddata,1)
    batch{1}.spm.stats.fmri_spec.dir = cellstr(fullfile(pwd,'glm'));
    batch{1}.spm.stats.fmri_spec.timing.units = 'scans';
    batch{1}.spm.stats.fmri_spec.timing.RT = 1;
    batch{1}.spm.stats.fmri_spec.sess.scans = {fullfile(pwd,'test_analyze_sub1.img')};
    batch{1}.spm.stats.fmri_spec.sess.cond.name = 'worddata';
    batch{1}.spm.stats.fmri_spec.sess.cond.onset = find(worddata(i,:));
    batch{1}.spm.stats.fmri_spec.sess.cond.duration = 1;
    spm_jobman('run',batch)
    load('glm/SPM.mat')
    regressors(:,i) = SPM.xX.X(:,1);
    cd(basedir)
    system('rm glm/*')
end
save('word_regressors.mat','regressors');

finn_path = '/Volumes/SOREN_SSD_2/';

fmridir = [finn_path 'intersubj_rsa/all_shen268_roi_ts']; 

files = dir(fullfile(fmridir,'*MOVIE1*ts.txt'));% this is the error - this should read *MOVIE1*gsr.txt

% do regression with lasso as in Huth et al. 2012

B = cell(1,length(files)); fitinfo = B; betaweights = B; 
for i = 1:length(files)
    fprintf([num2str(i) ' '])
    betaweights{i} = zeros(size(worddata,1),268);
    subdata = readtable(fullfile(files(i).folder,files(i).name),'ReadVariableNames',false);
    subdata = subdata{:,2:end};
    for q = 1:size(subdata,2)
%         rindx = randperm(size(subdata,1),floor(size(subdata,1)*0.8));
%         train = subdata(:,rindx); test = subdata(:,except(1:size(subdata,1),rindx));
%         B = lasso(regressors(rindx,:),train);
        %cvp = cvpartition(size(subdata,1),'HoldOut',0.2);
        [B{i}{q},fitinfo{i}{q}] = lasso(regressors,subdata(:,q),'CV',5,'NumLambda',10,'LambdaRatio',5e-3,'RelTol',5e-2);
        lamindx = fitinfo{i}{q}.IndexMinMSE;
        betaweights{i}(:,q) = B{i}{q}(:,lamindx);
    end
end

for i = 1:length(files)
    subdata = readtable(fullfile(files(i).folder,files(i).name),'ReadVariableNames',false);
    subdata = subdata{:,2:end};
    for q = 1:size(subdata,2)
        if ~any(betaweights{i}(:,q))
            [~,newlamindx] = min(fitinfo{i}{q}.MSE(any(B{i}{q},1)));
            betaweights{i}(:,q) = B{i}{q}(:,newlamindx);
        end
    end
end

save('sub_word_betaweights.mat','betaweights','-v7.3')
save('sub_lasso_results.mat','B','fitinfo','-v7.3')
