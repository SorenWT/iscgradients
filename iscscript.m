% script for Individual differences in naturalistic viewing revealed through 
% dimensionality reduction of inter-subject correlation matrices
% Soren Wainio-Theberge, 2020

finn_path = '/Volumes/SOREN_SSD_2';

clear allts
files = dir([finn_path '/intersubj_rsa/all_shen268_roi_ts/*MOVIE1*gsr.txt']); % replace the path with the path to that folder on your machine
for i = 1:length(files)
    tmp = readtable(fullfile(files(i).folder,files(i).name),'ReadVariableNames',false);
    allts{i} = tmp{:,2:end};
end
allts = cat(3,allts{:});
allts = permute(allts,[1 3 2]);

for i = 1:268
    iscmats(:,:,i) = corr(allts(:,:,i));
end
%isczmats = rtoz(iscmats): isczmats(find(isinf(isczmats))) = 0;
%cisczmats = mat2cell(isczmats,184,184,ones(1,268));

%load Shen atlas

shatlas = niftiread([finn_path '/intersubj_rsa/shen_2mm_268_parcellation.nii.gz']);
shinfo = niftiinfo([finn_path '/intersubj_rsa/shen_2mm_268_parcellation.nii.gz']);
shtable = readtable('shtable.xlsx'); % taken from Shen et al. 2015

% initial figure of group-level ISC

iscvals = ztor(mean(rtoz(belowDiag(iscmats)),1));

figure
set(gcf,'units','normalized','position',[0.2 0.2 0.6 0.8])
p = panel('no-manage-font');
p.pack('v',{60 40})
p(1).pack('h',{25 50 25})
p(2).pack('h',{1/3 1/3 1/3})
[p,~,ax,cbar] = plot_volonsurf(iscvals,shatlas,shinfo,'conte69','panel',p,'panelindx',{1 2},'colormap',hot);
cbar.FontSize = 14; cbar.Label.String = 'Inter-subject correlation'; cbar.Label.FontSize = 18;

savefig('ISC_grouplevel.fig'); export_fig('ISC_grouplevel.png','-m4')

% adbscan parameter validation

resvals = 50:50:250;
for i = 1:268
    tic
    for k = 1:5
        for res = 1:5
            for q = 1:100
                randset = randperm(184,184/2);
                train = rtoz(iscmats(:,randset,i)); train(train==Inf) = NaN;
                test = rtoz(iscmats(:,except(1:184,randset),i)); test(test==Inf) = NaN;
                cltrain = adbscan(1-corr(train','rows','pairwise'),k,25,resvals(res),0);
                cltest = adbscan(1-corr(test','rows','pairwise'),k,25,resvals(res),0);
                if ~any(isnan(cltrain)) && ~any(isnan(cltest))
                    %                   ri(k,res,q) = rand_index(cltrain,cltest,'adjusted');
                    trainadj = cltrain==cltrain'; testadj = cltest==cltest';
                    stab(k,res,q) = dice(belowDiag(trainadj),belowDiag(testadj));
                else
                    stab(k,res,q) = NaN;
                end
            end
        end
    end
    allstab(i,:,:,:) = permute(stab,[4 1 2 3]);
    toc
end

meanstab = nanmean(allstab,4);
meanstab = permute(meanstab,[2 3 1]);
pks = cell(1,268);
for i = 1:268
    tmp2 = findpeaks2d(meanstab(:,:,i));
    if ~isempty(tmp2)
        pks{i} = tmp2;
    end
end
find(cellfun(@isempty,pks,'uniformoutput',true))

% plotting average stability of each number of clusters over the cortex

shatlas = niftiread('shen_1mm_268_parcellation.nii.gz');
shinfo = niftiinfo('shen_1mm_268_parcellation.nii.gz');

p = panel('no-manage-font');

p.pack('v',{1/2 1/2})
p(1).pack('h',{0.1 0.4 0.4 0.1})
p(2).pack('h',{1/3 1/3 1/3})

plot_volonsurf(squeeze(nanmean(meanstab(1,:,:),2)),shatlas,shinfo,'conte69','panel',p,'panelindx',{1 2})
plot_volonsurf(squeeze(nanmean(meanstab(2,:,:),2)),shatlas,shinfo,'conte69','panel',p,'panelindx',{1 3})

for i = 1:3
    plot_volonsurf(squeeze(nanmean(meanstab(i+2,:,:),2)),shatlas,shinfo,'conte69','panel',p,'panelindx',{2 i})
end

Normalize_Clim(gcf,0)

cbar = findobj('parent',gcf,'type','colorbar');
cbar.Label.String = 'Mean stability'; cbar.Label.FontSize = 20; cbar.FontSize = 16;

savefig('Adbscan_stability.fig'); export_fig('Adbscan_stability.png','-m4')

% determining and plotting silhouette scores

for i = 1:268
    for ii = 1:5
        [~,bestres] = max(meanstab(ii,:,i));
        incorr = rtoz(iscmats(:,:,i)); incorr(incorr==Inf) = NaN;
        %adbclust(:,ii,i) = adbscan(1-corr(incorr,'rows','pairwise'),ii,25,resvals(bestres),0);
        incorr(isnan(incorr)) = 0;
        if ~any(isnan(adbclust(:,ii,i)))
            sil = silhouette(incorr(adbclust(:,ii,i)>0,adbclust(:,ii,i)>0),adbclust(adbclust(:,ii,i)>0,ii,i),'correlation');
            silscore(i,ii) = nanmean(sil);
        else
            silscore(i,ii) = 0;
        end
    end
end

p = panel('no-manage-font');

p.pack('v',{1/2 1/2})
p(1).pack('h',{25 50 25})
%p(1).pack('h',{0.1 0.4 0.4 0.1})
p(2).pack('h',{1/3 1/3 1/3})
set(gcf,'units','normalized','position',[0.15 0.4 0.7 0.6])

%plot_volonsurf(vect2vol(silscore(:,1),shatlas),shinfo,'conte69','panel',p,'panelindx',{1 2})
plot_volonsurf(silscore(:,2),shatlas,shinfo,'conte69','panel',p,'panelindx',{1 2})

for i = 1:3
    plot_volonsurf(silscore(:,i+2)+0.0001,shatlas,shinfo,'conte69','panel',p,'panelindx',{2 i})
end

Normalize_Clim(gcf,0)

cbar = findobj('parent',gcf,'type','colorbar');
delete(cbar(1)); delete(cbar(2)); delete(cbar(3)); cbar(1:3) = [];
cbar.Label.String = 'Silhouette score'; cbar.FontSize = 16; cbar.Label.FontSize = 20; 

savefig('Adbscan_silhouette.fig'); export_fig('Adbscan_silhouette.png','-m4')

% plot example parameter surface

surf(resvals,1:5,squeeze(nanmean(allstab(16,:,:,:),4))); 
xlabel('\epsilon step resolution'); ylabel('Number of clusters'); zlabel('Stability')
FixAxes(gca,14)
savefig('Adbscan_paramsurf_ex1.fig'); export_fig('Adbscan_paramsurf_ex1.png','-m4')
surf(resvals,1:5,squeeze(nanmean(allstab(75,:,:,:),4)))
xlabel('\epsilon step resolution'); ylabel('Number of clusters'); zlabel('Stability')
FixAxes(gca,14)
savefig('Adbscan_paramsurf_ex2.fig'); export_fig('Adbscan_paramsurf_ex2.png','-m4')

% version 2 using regular dbscan

eps = 0.2:0.025:1;
minpts = 25:5:50;

for q = 1:268
    for i = 1:length(eps)
        for ii = 1:length(minpts)
            matin = rtoz(iscmats(:,:,q)); matin(matin==Inf) = NaN;
            cl = dbscan(1-corr(matin,'rows','pairwise'),eps(i),minpts(ii),'distance','precomputed');
            maxclust(i,ii,q) = max(cl);
        end
    end
end

for i = 1:268
    numcl(i) = max(max(maxclust(:,:,i)));
end

% get stability and specificity

dbgenstab = zeros([size(maxclust) 100]);
dbstrictstab = zeros([size(maxclust),100]);
dbsilscore = zeros([size(maxclust)]);

clear trainadj testadj train test genstab strictstab
for i = 1:268
    if numcl(i) > 1
        multiclusts = find(maxclust(:,:,i) > 1);
        [in1,in2] = ind2sub(size(maxclust(:,:,i)),multiclusts);
        for ii = 1:length(in1)
            stabmat = zeros(25,25);
            for iii = 1:100
                cltrain = []; cltest = []; trainadj = []; testadj = [];
                randset = randperm(184,184/2);
                train = rtoz(iscmats(:,randset,i)); train(train==Inf) = NaN;
                test = rtoz(iscmats(:,except(1:184,randset),i)); test(test==Inf) = NaN;
                tmpeps = eps([max(in1(ii)-2,1) max(in1(ii)-1,1) in1(ii) min(in1(ii)+1,length(eps)) min(in1(ii)+2,length(eps))]);
                tmpminpts = minpts([max(in2(ii)-2,1) max(in2(ii)-1,1) in2(ii) min(in2(ii)+1,length(minpts)) min(in2(ii)+2,length(minpts))]);
                for q = 1:5
                    for qq = 1:5
                        cltrain(:,q,qq) = dbscan(1-corr(train','rows','pairwise'),tmpeps(q),tmpminpts(qq),'distance','precomputed');
                        cltest(:,q,qq) = dbscan(1-corr(test','rows','pairwise'),tmpeps(q),tmpminpts(qq),'distance','precomputed');
                        trainadj(:,q,qq) = belowDiag(cltrain(:,q,qq)==cltrain(:,q,qq)');
                        testadj(:,q,qq) = belowDiag(cltest(:,q,qq)==cltest(:,q,qq)');
                    end
                end
                if (max(cltrain(:,3,3))>1) && (max(cltest(:,3,3))>1)
                    strictstab(iii) = dice(trainadj(:,2,2),testadj(:,2,2));
                else
                    strictstab(iii) = 0;
                end
                trainadj = reshape(trainadj,size(trainadj,1),[]);
                testadj = reshape(testadj,size(testadj,1),[]);
                cltrain = reshape(cltrain,size(cltrain,1),[]);
                cltest = reshape(cltest,size(cltest,1),[]);
                
                
                for q = 1:size(trainadj,2)
                    for qq = 1:size(testadj,2)
                        stabmat(q,qq) = dice(trainadj(:,q),testadj(:,qq));
                        hascl(q,qq) = (max(cltrain(:,q))>1)&&(max(cltest(:,qq))>1);
                    end
                end
                genstab(iii) = max(max(stabmat.*hascl));
            end
            dbstrictstab(in1(ii),in2(ii),i,:) = strictstab;
            dbgenstab(in1(ii),in2(ii),i,:) = genstab;
            
            
            matin = rtoz(iscmats(:,:,i)); matin(matin==Inf) = NaN;
            cl = dbscan(1-corr(matin,'rows','pairwise'),eps(in1(ii)),minpts(in2(ii)),'distance','precomputed');
            
            matin(isnan(matin)) = 0;
            sil = silhouette(matin(cl>0,cl>0),cl(cl>0),'correlation'); % compute silhouette without noise points
            dbsilscore(in1(ii),in2(ii),i) = nanmean(sil);
        end
    end
end

for i = 1:268
    tmp = dbstrictstab; tmp(dbstrictstab==0)=NaN;
    regstrictstab(i) = max(max(nanmedian(tmp(:,:,i,:),4)));
    tmp = dbgenstab; tmp(dbgenstab==0)=NaN;
    reggenstab(i) = max(max(nanmedian(tmp(:,:,i,:),4)));
    regsilscore(i) = max(max(dbsilscore(:,:,i)));
end

%save('iscvars_adbscan.mat')

% figures for regular dbscan

mricro_glassbrain(vect2vol(numcl,shatlas),shinfo)
mricro_glassbrain(vect2vol(regstrictstab,shatlas),shinfo)
mricro_glassbrain(vect2vol(reggenstab,shatlas),shinfo)
mricro_glassbrain(vect2vol(regsilscore,shatlas),shinfo)
mricro_glassbrain(vect2vol(double(reggenstab>0.5 & regsilscore>0.5),shatlas),shinfo)

plot_volonsurf(numcl,shatlas,shinfo,'conte69','mask',numcl>1,'CLim',[0 5])
colormap hot;
savefig('numcl_surf.fig'); export_fig('numcl_surf.png','-m4');

regstrictstab(isnan(regstrictstab)) = 0;
plot_volonsurf(regstrictstab+0.001,shatlas,shinfo,'conte69','mask',regstrictstab>0,...
    'CLim',[min(regstrictstab(regstrictstab>0)) max(regstrictstab(regstrictstab>0))])
colormap hot;
cbar = findobj('parent',gcf,'type','colorbar');
cbar.Label.String = 'Mean stability'; cbar.FontSize = 16; cbar.Label.FontSize = 20; 
savefig('regstrictstab_surf.fig'); export_fig('regstrictstab_surf.png','-m4');
regstrictstab(regstrictstab==0) = NaN;

reggenstab(isnan(reggenstab)) = 0;
plot_volonsurf(reggenstab+0.001,shatlas,shinfo,'conte69','mask',reggenstab>0,...
    'CLim',[min(reggenstab(reggenstab>0)) max(reggenstab(reggenstab>0))])
colormap hot;
cbar = findobj('parent',gcf,'type','colorbar');
cbar.Label.String = 'Mean stability'; cbar.FontSize = 16; cbar.Label.FontSize = 20; 
savefig('reggenstab_surf.fig'); export_fig('reggenstab_surf.png','-m4');
reggenstab(reggenstab==0) = NaN;

regsilscore(isnan(regsilscore)) = 0;
plot_volonsurf(regsilscore+0.001,shatlas,shinfo,'conte69','mask',regsilscore>0,...
    'CLim',[min(regsilscore(regsilscore>0)) max(regsilscore(regsilscore>0))])
cbar = findobj('parent',gcf,'type','colorbar');
cbar.Label.String = 'Silhouette score'; cbar.FontSize = 16; cbar.Label.FontSize = 20; 
savefig('regsilscore_surf.fig'); export_fig('regsilscore_surf.png','-m4');
regsilscore(regsilscore==0) = NaN;

plot_volonsurf(double(reggenstab>0.5 & regsilscore>0.5)+1,shatlas,shinfo,'conte69',...
    'mask',(reggenstab>0.5 & regsilscore>0.5),'CLim',[0 5])
colormap hot;
savefig('goodclusts_surf.fig'); export_fig('goodclusts_surf.png','-m4');

% plot example clusterings with ADBSCAN and regular DBSCAN


subplot(1,2,1)
cbar = clustplot(iscmats(:,:,75),adbclust(:,2,75));
cbar.Label.String = 'ISC (Pearson''s r)'; cbar.Label.FontSize = 14; cbar.Label.FontSize = 16;
axis square

subplot(1,2,2)
clustplot(iscmats(:,:,191),adbclust(:,2,191))
cbar.Label.String = 'ISC (Pearson''s r)'; cbar.Label.FontSize = 14; cbar.Label.FontSize = 16;
axis square

set(gcf,'color','w')
savefig('Adbclust_ex.fig'); export_fig('Adbclust_ex.png','-m4')

subplot(1,2,1)
tmp = find(dbsilscore(:,:,147) == max(max(dbsilscore(:,:,147))));
[e1,m1] = ind2sub([33 6],tmp);
inmat = rtoz(iscmats(:,:,147)); inmat(inmat==Inf) = NaN;
cl = dbscan(1-corr(inmat,'rows','pairwise'),eps(e1),minpts(m1),'distance','precomputed');
cbar = clustplot(iscmats(:,:,147),cl);
cbar.Label.String = 'ISC (Pearson''s r)'; cbar.Label.FontSize = 14; cbar.Label.FontSize = 16;
axis square

subplot(1,2,2)
e1 = 15; m1 = 4;
inmat = rtoz(iscmats(:,:,147)); inmat(inmat==Inf) = NaN;
cl = dbscan(1-corr(inmat,'rows','pairwise'),eps(e1),minpts(m1),'distance','precomputed');
cbar = clustplot(iscmats(:,:,147),cl);
cbar.Label.String = 'ISC (Pearson''s r)'; cbar.Label.FontSize = 14; cbar.Label.FontSize = 16;
axis square

savefig('Dbclust_ex.fig'); export_fig('Dbclust_ex.png','-m4')

% PCA and dimensionality reduction

isczmats = rtoz(iscmats);
isczmats(isinf(isczmats)) = 0;

gpca = GradientMaps('kernel','pearson','approach','pca');
gpca = gpca.fit(mat2cell(isczmats,184,184,ones(1,268)),'sparsity',0);

gshuf = GradientMaps('kernel','pearson','approach','pca');

clear eigshuf_pca
parfor i = 1:268
    tic
    iscshuf = cell(1,1000);
    for n = 1:1000
        iscshuf{n} = zeros(size(isczmats(:,:,i)));
        for dim = 1:size(isczmats(:,:,i),2)
            iscshuf{n}(:,dim) = isczmats(randperm(184),dim,i);
        end
    end
    gshufres = gshuf.fit(iscshuf,'sparsity',0);
    toc
    eigshuf_pca{i} = [gshufres.lambda{:}];
end
eigshuf_pca = cat(3,eigshuf_pca{:});

lambda_pca = [gpca.lambda{:}];

for i = 1:268
    lambda_pca_p(:,i) = 1-mean(lambda_pca(:,i) > eigshuf_pca(:,:,i),2);
end

gdm = GradientMaps('kernel','pearson','approach','dm');
gdm = gdm.fit(mat2cell(isczmats,184,184,ones(1,268)),'sparsity',0);

gshuf = GradientMaps('kernel','pearson','approach','dm');

clear eigshuf_dm
eigshuf_dm = cell(1,268);
parfor i = 1:268
    iscshuf = cell(1,1000);
    for n = 1:1000
        iscshuf{n} = zeros(size(isczmats(:,:,i)));
        for dim = 1:size(isczmats(:,:,i),2)
            iscshuf{n}(:,dim) = isczmats(randperm(184),dim,i);
        end
    end
    gshufres = gshuf.fit(iscshuf,'sparsity',0);
    eigshuf_dm{i} = [gshufres.lambda{:}];
end
eigshuf_dm = cat(3,eigshuf_dm{:});

lambda_dm = [gdm.lambda{:}];

for i = 1:268
    lambda_dm_p(:,i) = 1-sum(lambda_dm(:,i) > eigshuf_dm(:,:,i),2)/100;
end

% gradient alignment: joint embedding within a priori networks

cisczmats = mat2cell(isczmats,184,184,ones(1,268));

gpca_proc = GradientMaps('kernel','pearson','approach','pca','alignment','procrustes','n_components',5);
gpca_joint = GradientMaps('kernel','pearson','approach','pca','alignment','joint','n_components',5);

gdm_proc = GradientMaps('kernel','pearson','approach','dm','alignment','procrustes','n_components',5);
gdm_joint = GradientMaps('kernel','pearson','approach','dm','alignment','joint','n_components',5);

pca_allalign = gpca_proc.fit(cisczmats,'niterations',100,'sparsity',0);
pca_allgrad = cat(3,pca_allalign.aligned{:});

dm_allalign = gdm_proc.fit(cisczmats,'niterations',100,'sparsity',0);
dm_allgrad = cat(3,dm_allalign.aligned{:});

pca_allgrad = permute(pca_allgrad,[1 3 2]);
dm_allgrad = permute(dm_allgrad,[1 3 2]);

for i = 1:8
    pca_joint{i} = gpca_joint.fit(cisczmats(find(shtable.Network==i)),'sparsity',0);
    pca_netgrad_joint{i} = cat(3,pca_joint{i}.aligned{:});
    pca_proc{i} = gpca_proc.fit(cisczmats(find(shtable.Network==i)),'niterations',100,'sparsity',0);
    pca_netgrad_proc{i} = cat(3,pca_proc{i}.aligned{:});
    
    dm_joint{i} = gdm_joint.fit(cisczmats(find(shtable.Network==i)),'sparsity',0);
    dm_netgrad_joint{i} = cat(3,dm_joint{i}.aligned{:});
    dm_proc{i} = gdm_proc.fit(cisczmats(find(shtable.Network==i)),'niterations',100,'sparsity',0);
    dm_netgrad_proc{i} = cat(3,dm_proc{i}.aligned{:});
end

% scree plots

%load('lkcmap2')
%monored = lkcmap2(49:end,:);
%monoblue = lkcmap2(1:18,:);

figure
p = panel('no-manage-font');
p.pack('h',{1/2 1/2})
set(gcf,'units','normalized','position',[0.15 0.4 0.7 0.6])

lambda_pca = [pca_allalign.lambda{:}];

p(1).select()
mask = double(lambda_pca_p<0.05); mask(find(mask==0)) = NaN;
plot(lambda_pca,'linewidth',0.5,'linestyle','--')
hold on
plot(lambda_pca.*mask,'linewidth',2)
plot(nanmean(eigshuf_pca,3),'-.ok','linewidth',5)
%stdshade(1:183,nanmean(eigshuf_pca,3)','k',0.3,2,'prctileci');
set(gca,'XLim',[0 10],'XTick',1:10)
xlabel('Component')
ylabel('Percent variance explained')
title('PCA')
FixAxes(gca,18)

lambda_dm = [dm_allalign.lambda{:}];

p(2).select()
mask = double(lambda_dm_p<0.05); mask(find(mask==0)) = NaN;
plot(lambda_dm,'linewidth',0.5,'linestyle','--')
hold on
plot(lambda_dm.*mask,'linewidth',2)
plot(nanmean(eigshuf_dm,3),'-.ok','linewidth',5)
%stdshade(1:183,nanmean(eigshuf_pca,3)','k',0.3,2,'prctileci');
set(gca,'XLim',[0 10],'XTick',1:10)
xlabel('Component')
ylabel('Eigenvalue')
title('Diffusion maps')
FixAxes(gca,18)

p.marginbottom = 20; p.marginleft = 24; p.marginright = 8; p.margintop = 15;
savefig('Scree_pca_dm.fig'); export_fig('Scree_pca_dm.png','-m4')

netnames = {'Medial frontal','Frontoparietal','Default mode','Subcortical-cerebellum',...
    'Motor','Visual I','Visual II','Visual association'}; % for now

% variance explained topoplots

figure
p = panel('no-manage-font');
p.pack(3,2)
set(gcf,'units','normalized','position',[0.2 0 0.6 1])

for i = 1:3
    [p,~,ax1{i},cbar] = plot_volonsurf(lambda_pca(i,:),shatlas,shinfo,'conte69','panel',p,'panelindx',{i 1});
    if i ~=1
        delete(cbar)
    else
        cbar.Label.String = '% Variance Explained'; cbar.FontSize = 14; cbar.Label.FontSize = 16;
    end
    
    [p,~,ax2{i},cbar] = plot_volonsurf(lambda_dm(i,:),shatlas,shinfo,'conte69','panel',p,'panelindx',{i 2});
    if i ~=1
        delete(cbar)
    else
        cbar.Label.String = 'Eigenvalue'; cbar.FontSize = 14; cbar.Label.FontSize = 16;
    end
end
Normalize_Clim([ax1{:}],0)
Normalize_Clim([ax2{:}],0)
p(1).marginbottom = 28; p(2).marginbottom = 28;
p(1,2).marginleft = 40; p(2,2).marginleft = 40; p(3,2).marginleft = 40;
colormap hot

savefig('Component_eigenvalues_topo.fig'); export_fig('Component_eigenvalues_topo.png','-m4')

% correlation of eigenvalues with mean ISC values

W = weights_from_atlas(shatlas,shinfo);
MEM = compute_mem(W);
eig_rand_pca = moran_randomization(lambda_pca(1,:)',W,10000);
eig_rand_dm = moran_randomization(lambda_dm(1,:)',W,10000);

rho_mean_eigpca = corr(lambda_pca(1,:)',iscvals','type','spearman');
pcorr_mean_eigpca = 1-mean(rho_mean_eigpca>corr(iscvals',squeeze(eig_rand_pca),'type','spearman'));

rho_eigpca_eigdm = corr(lambda_pca(1,:)',lambda_dm(1,:)','type','spearman');
pcorr_eigpca_eigdm = 1-mean(rho_eigpca_eigdm>corr(lambda_dm(1,:)',squeeze(eig_rand_pca),'type','spearman'));

rho_mean_eigdm = corr(lambda_dm(1,:)',iscvals','type','spearman');
pcorr_mean_eigdm = 1-mean(rho_mean_eigdm>corr(iscvals',squeeze(eig_rand_dm),'type','spearman'));


% alignment figure

l = flag; l = l([3 1],:);

figure
p = panel('no-manage-font');
set(gcf,'units','normalized','position',[0.05 0 0.9 1])
p.pack('h',{1/2 1/2})
p(1).pack('v',{30 70})
p(2).pack('v',{30 70})
p(1,2).pack(2,4)
p(2,2).pack(2,4)

p(1,1).select()
histogram(belowDiag(corr(pca_allgrad(:,:,1),'type','pearson')));
hold on
yl = get(gca,'YLim');
line([median(belowDiag(corr(pca_allgrad(:,:,1),'type','pearson'))) median(belowDiag(corr(pca_allgrad(:,:,1),'type','pearson')))],...
    yl,'LineWidth',4,'color',l(1,:),'handlevisibility','off','linestyle','--')
set(gca,'YLim',yl)

xlabel('Spearman''s \rho')
ylabel('Count')
title('PCA Whole-brain alignment')
FixAxes(gca,16)

ord = [1 2 3 4 1 2 3 4];
for i = 1:8
    p(1,2,ceil(i/4),ord(i)).select()
    histogram(belowDiag(corr(squeeze(pca_netgrad_proc{i}(:,1,:)),'type','pearson')));
    hold on
    histogram(belowDiag(corr(squeeze(pca_netgrad_joint{i}(:,1,:)),'type','pearson')));
    
    yl = get(gca,'YLim');
    line([median(belowDiag(corr(squeeze(pca_netgrad_proc{i}(:,1,:)),'type','pearson'))) median(belowDiag(corr(squeeze(pca_netgrad_proc{i}(:,1,:)),'type','pearson')))],...
        yl,'LineWidth',4,'color',l(1,:),'handlevisibility','off','linestyle','--')
    line([median(belowDiag(corr(squeeze(pca_netgrad_joint{i}(:,1,:)),'type','pearson'))) median(belowDiag(corr(squeeze(pca_netgrad_joint{i}(:,1,:)),'type','pearson')))],...
        yl,'LineWidth',4,'color',l(2,:),'handlevisibility','off','linestyle','--')
    set(gca,'YLim',yl)
    
    if i == 4
        legend({'Procrustes alignment','Joint embedding'})
    end
    xlabel('Correlation')
    ylabel('Count')
    title(netnames{i})
    FixAxes(gca,14)
end
p.marginleft = 24; p.margintop = 10;
p(1).marginbottom = 28;
p(1,2).de.margin = [22 25 5 5];

p(2,1).select()
histogram(belowDiag(corr(dm_allgrad(:,:,1),'type','pearson')));

hold on
yl = get(gca,'YLim');
line([median(belowDiag(corr(dm_allgrad(:,:,1),'type','pearson'))) median(belowDiag(corr(dm_allgrad(:,:,1),'type','pearson')))],...
    yl,'LineWidth',4,'color',l(1,:),'handlevisibility','off','linestyle','--')
set(gca,'YLim',yl)

xlabel('Spearman''s \rho')
ylabel('Count')
title('Diffusion mapping whole-brain alignment')
FixAxes(gca,16)

ord = [1 2 3 4 1 2 3 4];
for i = 1:8
    p(2,2,ceil(i/4),ord(i)).select()
    histogram(belowDiag(corr(squeeze(dm_netgrad_proc{i}(:,1,:)),'type','pearson')));
    hold on
    histogram(belowDiag(corr(squeeze(dm_netgrad_joint{i}(:,1,:)),'type','pearson')));
    
    yl = get(gca,'YLim');
    line([median(belowDiag(corr(squeeze(dm_netgrad_proc{i}(:,1,:)),'type','pearson'))) median(belowDiag(corr(squeeze(dm_netgrad_proc{i}(:,1,:)),'type','pearson')))],...
        yl,'LineWidth',4,'color',l(1,:),'handlevisibility','off','linestyle','--')
    line([median(belowDiag(corr(squeeze(dm_netgrad_joint{i}(:,1,:)),'type','pearson'))) median(belowDiag(corr(squeeze(dm_netgrad_joint{i}(:,1,:)),'type','pearson')))],...
        yl,'LineWidth',4,'color',l(2,:),'handlevisibility','off','linestyle','--')
    set(gca,'YLim',yl)
    
    if i == 4
        legend({'Procrustes alignment','Joint embedding'})
    end
    xlabel('Correlation')
    ylabel('Count')
    title(netnames{i})
    FixAxes(gca,14)
end
p(2).marginbottom = 28;
p(2,2).de.margin = [22 25 5 5];

savefig('Alignment_networks.fig'); export_fig('Alignment_networks.png','-m4')

% example ISC matrices sorted by gradients

figure 

set(gcf,'units','normalized','position',[0 0 1 1])

p = panel('no-manage-font');
p.pack('v',{30 70})
p(1).pack('h',{25 50 25})
p(2).pack(2,3)
reginds = [191 75];

[p,~,ax,cbar] = plot_volonsurf(vect2vol(ind2log(reginds,268),shatlas),shinfo,'conte69','panel',p,'panelindx',{1 2});
delete(cbar)
 
for i = 1:2
    for q = 1:3
        p(2,i,q).pack('h',{10 90})
        p(2,i,q,1).select()
        imagesc(pca_allalign.gradients{reginds(i)}(:,q))
        set(gca,'colormap',jet,'XTick',[],'YDir','reverse','YLim',[0.5 184.5])
        ylabel(['Gradient ' num2str(q) ' value (unsorted)'])
        p(2,i,q,1).margin = [5 5 5 5];
        
        p(2,i,q,2).select()
        [~,sortindx] = sort(pca_allalign.gradients{reginds(i)}(:,q));
        imagesc(iscmats(sortindx,sortindx,reginds(i)))
        axis square
        set(gca,'CLim',[min(belowDiag(iscmats(:,:,reginds(i)))) max(belowDiag(iscmats(:,:,reginds(i))))],...
            'YDir','reverse','XLim',[0.5 184.5],'YLim',[0.5 184.5]);
        if q == 3
            cbar = colorbar;
            cbar.Label.String = 'ISC (Pearson''s r)'; cbar.Label.FontSize = 16; cbar.FontSize = 14;
        end
    end
end

savefig('PCA_exsort.fig'); export_fig('PCA_exsort.png','-m4')


% pairwise procrustes - pca and dm

pcapairs = cell(268,268); 
proccorrmat_pca = [];
for i = 1:268
    for ii= 1:(i-1)
        pcapairs{i,ii} = gpca_proc.fit([cisczmats(i) cisczmats(ii)],'sparsity',0);
        proccorrmat_pca(i,ii,:) = diag(corr(pcapairs{i,ii}.aligned{1},pcapairs{i,ii}.aligned{2},'type','pearson'));
    end
end
proccorrmat_pca = cat(2,proccorrmat_pca,zeros(268,1,size(proccorrmat_pca,3)));
proccorrmat_pca = proccorrmat_pca+t3d(proccorrmat_pca)+repmat(eye(268),1,1,size(proccorrmat_pca,3));

% took too long
% dmpairs = cell(268,268);
% for i = 1:268
%     for ii= 1:(i-1)
%         dmpairs{i,ii} = gdm_proc.fit([cisczmats(i) cisczmats(ii)],'sparsity',0);
%         proccorrmat_dm(i,ii,:) = diag(corr(dmpairs{i,ii}.aligned{1},dmpairs{i,ii}.aligned{2},'type','pearson'));
%     end
% end
% proccorrmat_dm = cat(2,proccorrmat_dm,zeros(268,1,size(proccorrmat_dm,3)));
% proccorrmat_dm = proccorrmat_dm+t3d(proccorrmat_dm)+repmat(eye(268),1,1,size(proccorrmat_dm,3));

% focus on component 1, pca for now

for i = 1:268
    for ii = 1:268
        iscsubscorr(i,ii) = corr(rtoz(belowDiag(iscmats(:,:,i))),rtoz(belowDiag(iscmats(:,:,ii))));
    end
end

for k = 2:10
    for i = 1:100
        randset = randperm(268,268/2);
        %train = rtoz(iscsubscorr(:,randset)); train(train==Inf) = 0;
        %test = rtoz(iscsubscorr(:,except(1:268,randset))); test(test==Inf) = 0;
        train = rtoz(proccorrmat_pca(:,randset,1)); train(train==Inf) = 0;
        test = rtoz(proccorrmat_pca(:,except(1:268,randset),1)); test(test==Inf) = 0;
        cltrain = kmeans(train,k,'distance','correlation');
        cltest = kmeans(test,k,'distance','correlation');
        %cltrain = spectralcluster(train,k,'distance','precomputed');
        %cltest = spectralcluster(test,k,'distance','precomputed');
        
        trainadj = cltrain==cltrain'; testadj = cltest==cltest';
        stab_pca(k,i) = dice(belowDiag(trainadj),belowDiag(testadj));
    end
    %cl = kmeans(iscsubscorr,k,'distance','correlation');
    %sil = silhouette(iscsubscorr,cl,'correlation');
    %silscore_pca(k) = mean(sil);
    
    tmp = rtoz(proccorrmat_pca(:,:,1)); tmp(tmp==Inf) = 0;
    cl(:,k) = consensus_clust(@kmeans,tmp,k,100,'distance','correlation');
    sil = silhouette(proccorrmat_pca(:,:,1),cl(:,k),'correlation');
    silscore_pca(k) = mean(sil);
end

for k = 8
    %cl1 = consensus_clust(@kmeans,rtoz(proccorrmat_pca(:,:,1),'zeroinf'),k,100,'distance','correlation');
    %cl2 = consensus_clust(@kmeans,rtoz(proccorrmat_pca(:,:,2),'zeroinf'),k,100,'distance','correlation');
    %adj1 = cl1==cl1'; adj2 = cl2==cl2'; meanadj = mean(cat(3,adj1,adj2),3);
    %cl(:,k) = consensus_clust(@kmeans,meanadj,k,100,'distance','correlation');
    
    for i = 1:k
        %pca_cl_joint{i} = gpca_joint.fit(cisczmats(find(cl==i)),'sparsity',0);
        %pca_cl_netgrad_joint{i} = cat(3,pca_cl_joint{i}.aligned{:});
        meancorrmat = ztor(mean(rtoz(proccorrmat_pca(cl(:,k)==i,cl(:,k)==i,1:2)),3));
        [~,bestcentr] = max(nansum(meancorrmat,2));
        tmp = find(cl(:,k)==i);
        pca_cl_proc{k}{i} = gpca_proc.fit(cisczmats(find(cl(:,k)==i)),'niterations',100,'sparsity',0,'reference',pca_allalign.gradients{tmp(bestcentr)}(:,1:5));
        pca_cl_netgrad_proc{k}{i} = cat(3,pca_cl_proc{k}{i}.aligned{:});
        pca_cl_joint{k}{i} = gpca_joint.fit(cisczmats(find(cl(:,k)==i)),'sparsity',0);
        pca_cl_netgrad_joint{k}{i} = cat(3,pca_cl_joint{k}{i}.aligned{:});
    end
    %tmp1 = []; tmp2 = [];
    for i = 1:k
       intraclust_corr_comp1_proc{k}{i} = belowDiag(corr(squeeze(pca_cl_netgrad_proc{k}{i}(:,1,:)),'type','pearson'));
       intraclust_corr_comp2_proc{k}{i} = belowDiag(corr(squeeze(pca_cl_netgrad_proc{k}{i}(:,2,:)),'type','pearson'));
       intraclust_corr_comp1_joint{k}{i} = belowDiag(corr(squeeze(pca_cl_netgrad_joint{k}{i}(:,1,:)),'type','pearson'));
       intraclust_corr_comp2_joint{k}{i} = belowDiag(corr(squeeze(pca_cl_netgrad_joint{k}{i}(:,2,:)),'type','pearson'));
    end
    %intraclust_corr_comp1{k} = ztor(cellfun(@(d) mean(rtoz(d)),tmp1,'uniformoutput',true));
    %intraclust_corr_comp2{k} = ztor(cellfun(@(d) mean(rtoz(d)),tmp2,'uniformoutput',true));
%     intraclust_meancorr(k) = (ztor(mean(rtoz(cat(1,tmp1{:}))))+ztor(mean(rtoz(cat(1,tmp2{:})))))/2;
%     allcomps = cat(3,pca_cl_netgrad_proc{k}{:});
%     sil1(:,k) = silhouette(squeeze(allcomps(:,1,:))',cl(:,k),'correlation');
%     sil2(:,k) = silhouette(squeeze(allcomps(:,2,:))',cl(:,k),'correlation');
%     meansilscore(k) = (mean(sil1(:,k))+mean(sil2(:,k)))/2;
end

% big figure on PCA alignment
figure
set(gcf,'units','normalized','position',[0 0 1 1])
p = panel('no-manage-font');
p.pack('v',{30 70})
p(1).pack('h',{25 50 25})

p(1,2).select()
imagesc(proccorrmat_pca(:,:,1))
title('PCA alignment matrix','FontSize',18)
axis square
set(gca,'YDir','reverse','CLim',...
    [min(belowDiag(proccorrmat_pca(:,:,1))) max(belowDiag(proccorrmat_pca(:,:,1)))],...
    'XLim',[0.5 268.5],'YLim',[0.5 268.5])
cbar = colorbar; cbar.Label.String = 'Pearson''s r'; cbar.Label.FontSize = 18; cbar.FontSize = 14;

p(1,3).pack('h',{1/2 1/2})
p(1,3,1).select()
imagesc(iscsubscorr)
title('Overall ISC similarity matrix')
axis square
set(gca,'YDir','reverse','CLim',...
    [min(belowDiag(iscsubscorr)) max(belowDiag(iscsubscorr))],...
    'XLim',[0.5 268.5],'YLim',[0.5 268.5])
cbar = colorbar; cbar.Label.String = 'Pearson''s r'; cbar.Label.FontSize = 16; cbar.FontSize = 12;


p(2).pack('h',{1/2 1/2})
p(2,1).pack('h',{25 75})
p(2,1,1).pack('v',{25 50 25})
p(2,1,1,2).select()
yyaxis left
plot(mean(stab_pca,2),'LineWidth',2)
ylabel('Stability (mean dice coefficient)')
xlabel('Number of clusters')
hold on
yyaxis right
plot(silscore_pca,'LineWidth',2)
ylabel('Silhouette score')
legend({'Stability','Silhouette score'})
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'FontSize'    , 14        , ...
    'LineWidth'   , 1.5         , ...
    'TickLength',0.5*[0.01 0.025]);

p(2,1,2).pack('v',{1/3 1/3 1/3})
clnums = [3 5 8];

for i = 1:3
    p(2,1,2,i).pack('h',{60 40})
    %inmat = rtoz(proccorrmat_pca(:,:,1)); inmat(inmat ==Inf)= 0;
    %cl = consensus_clust(@kmeans,inmat,clnums(i),100,'distance','correlation');
    [p,~,ax,cbar] = plot_volonsurf(cl(:,clnums(i)),shatlas,shinfo,'conte69','panel',p,'panelindx',{2 1 2 i 1},'colormap',jet);
    delete(cbar)
    
    p(2,1,2,i,2).select()
    clustplot(proccorrmat_pca(:,:,1),cl(:,clnums(i)))
    axis square
    title([num2str(clnums(i)) ' cluster solution'])
end

p(2,2).pack(4,2) % use 8 cluster solution

for i = 1:8
    p(2,2,ceil(i/2),2-mod(i,2)).select()
    histogram(intraclust_corr_comp1_proc{8}{i}); 
    hold on
    histogram(intraclust_corr_comp1_joint{8}{i});
    yl = get(gca,'YLim');
    line([median(intraclust_corr_comp1_proc{8}{i}) median(intraclust_corr_comp1_proc{8}{i})],yl,'LineWidth',4,'color',l(1,:),'handlevisibility','off','linestyle','--')
        line([median(intraclust_corr_comp1_joint{8}{i}) median(intraclust_corr_comp1_joint{8}{i})],yl,'LineWidth',4,'color',l(2,:),'handlevisibility','off','linestyle','--')
    
    if i==1
        legend({'Procrustes alignment','Joint embedding'})
    end
    FixAxes(gca,14)
    xlabel('Correlation'); ylabel('Count');
    title(['Cluster ' num2str(i)])
end

p(1).marginbottom = 40;

p(2,2).marginleft = 40; 
p(2,1,2).marginleft = 25;
p.marginleft = 20; p.margintop = 8;
p(2,2).de.marginbottom = 24;

savefig('PCA_alignment_fig_pearson.fig'); export_fig('PCA_alignment_fig_pearson.png','-m4')

% new way: pca of all significant gradients

allgrads = cat(3,pca_allalign.gradients{:});
tmplambda = lambda_pca(1:5,:); tmpp = lambda_pca_p(1:5,:);
allgrads = permute(allgrads,[2 3 1]);

mask = tmpp<0.05;

allgrads = allgrads.*repmat(mask,1,1,size(allgrads,3));
allgrads = permute(allgrads,[3 1 2]);
allgrads = reshape(allgrads,184,[]);
allgrads(:,~any(allgrads,1)) = [];

[gradweights,gradcomps,~,~,expl] = pca(zscore(allgrads,[],1));


for i = 1:10000
    shuf = zeros(size(allgrads));
    for q = 1:size(allgrads,2)
       shuf(:,q) = allgrads(randperm(size(allgrads,1)),q); 
    end
    [~,~,~,~,shufexpl(:,i)] = pca(zscore(shuf,[],1));
end
gradcomps_p = 1-mean(expl>shufexpl,2);
%kaisercrit = expl>(100/size(allgrads,2));
%gradcomps_p = gradcomps_p.*kaisercrit;
cutoff = find(gradcomps_p>=0.05,1);
gradcomps = gradcomps(:,1:cutoff-1);

allgrads_raw = cat(3,pca_allalign.gradients{:});

clear loadings_z loadings_p loadingsboot loadings_raw
for i = 1:268
    loadings_raw(:,:,i) = corr(gradcomps,allgrads_raw(:,:,i));
    for q = 1:1000
        bootindx = ceil(rand(size(gradcomps,1),1)*184);
        loadingsboot(:,:,q) = corr(gradcomps(bootindx,:),allgrads_raw(bootindx,:,i));
    end
    %loadings_p(:,:,i) = 1-sum(loadings_raw(:,:,i)>loadingshuf,3)/100;
    loadings_z(:,:,i) = loadings_raw(:,:,i)./std(loadingsboot,[],3);
    loadings_p(:,:,i) = ztop(loadings_z(:,:,i));
    %loadings_sig(:,:,i) = (prctile(loadingsboot,2.5,3)>0 & mean(loadingsboot,3)>0) | (prctile(loadingsboot,97.5,3)<0 & mean(loadingsboot,3)<0);
end
%loadings = loadings_z.*(loadings_p<0.05);
%loadings = loadings_raw.*loadings_sig;

% figure to display the new components and plot the loadings across the brain

figure
set(gcf,'units','normalized','position',[0 0 1 1]);

p = panel('no-manage-font');
p.pack('h',{20 60 20})
p(1).pack('v',{25 50 25})
p(2).pack('v',{1/3 1/3 1/3})
p(3).pack('v',{1/3 1/3 1/3})

p(1,2).select()
mask = double(gradcomps_p<0.05 & kaisercrit); mask(mask==0) = NaN;
plot(expl,'linewidth',0.5,'linestyle','--')
hold on
plot(expl.*mask,'linewidth',3)
plot(nanmean(shufexpl,2),'-.ok','linewidth',3)
set(gca,'XLim',[1 40])
FixAxes(gca,20)
xlabel('Component')
ylabel('Variance explained (%)')

clear ax
for i = 1:3
   p(2,i).pack('h',{1/3 1/3 1/3})
   for ii = 1:3
        [p,~,ax{i,ii},cbar] = plot_volonsurf(abs(squeeze(loadings(i,ii,:)))+0.001,shatlas,shinfo,'conte69',...
            'panel',p,'panelindx',{2 i ii},'mask',fdr(loadings_p(i,ii,:))<0.05,'colormap',hot);
        if i~=1 || ii~=1
           delete(cbar) 
        else
            cbar.FontSize = 18; cbar.Label.String = 'Loading (absolute value)'; cbar.Label.FontSize = 20;
        end
        p(2,i,ii).de.margin = [18 18 5 5];
        p(2,i,ii).marginleft = 30;
   end
   p(2,i).marginbottom = 30;
   
   p(3,i).select()
   [~,maxload] = max(max(loadings(i,:,:),[],2),[],3);
   [~,sortindx] = sort(gradcomps(:,i));
   iscproj = isczmats(:,:,maxload); 
   if i > 1
       for q = 1:size(iscproj,2)
            iscproj(:,q) = regresid(gradcomps(:,1:(i-1)),iscproj(:,q)); 
       end
   end
   iscproj = ztor(iscproj); iscproj = iscproj+eye(size(iscproj));
   imagesc(iscproj(sortindx,sortindx))
   set(gca,'YDir','reverse','XLim',[0.5 184.5]','YLim',[0.5 184.5],...
       'CLim',[min(belowDiag(iscproj)) max(belowDiag(iscproj))],'XTick',[],'YTick',[])
   cbar = colorbar; 
   cbar.FontSize = 14; cbar.Label.FontSize = 16; cbar.Label.String = 'ISC'; 
   axis square
   title(['Canonical gradient ' num2str(i)],'FontSize',14)
end
ax = cat(2,ax{1,:},ax{2,:},ax{3,:});
%ax = cat(2,ax{:});
Normalize_Clim(ax,0)
p(2).marginleft = 40; 
p(3).marginleft = 40;
p(3).de.margintop = 20; p(3).de.marginbottom = 20; % make them smaller
p.margintop = 20; p.marginright = 20; p.marginleft = 22;

savefig('PCA_canonicalcomponents_loadings.fig'); export_fig('PCA_canonicalcomponents_loadings.png','-m4');

% load in beta weights and do PCA on these

make_word_regressors

load('sub_word_betaweights.mat');

betaweights = cat(3,betaweights{:});
betaweights = permute(betaweights,[3 1 2]);

cvp = cvpartition(184,'HoldOut',0.25);
tr = training(cvp); tst = test(cvp);
clear rpredlasso ppredlasso gradpredweights
gradpredweights = zeros(859,3,268);
for q = 1:3 % focus on top 3 ISC gradients for now
    tmpgradweights = cell(1,268);
    tmprpredlasso = zeros(1,268); tmpppredlasso = zeros(1,268);
    parfor i = 1:268
        tmpgradweights{i} = zeros(859,1);
        [B,fitinfo] = lasso(betaweights(tr,:,i),gradcomps(tr,q),'CV',10,'NumLambda',20,'RelTol',1e-3);
        tmpgradweights{i}(:,1) = B(:,fitinfo.IndexMinMSE);
        if ~any(tmpgradweights{i}(:,1))
            [~,newlamindx] = min(fitinfo.MSE(any(B,1)));
            tmpgradweights{i}(:,1) = B(:,newlamindx); 
        end
        pred = betaweights(tst,:,i)*tmpgradweights{i}(:,1);
        [tmprpredlasso(i),tmpppredlasso(i)] = corr(pred,gradcomps(tst,q));
    end
    gradpredweights(:,q,:) = cat(3,tmpgradweights{:});
    rpredlasso(q,:) = tmprpredlasso; ppredlasso(q,:) = tmpppredlasso;
end

%now get the words

synsets = h5read('WordNetFeatures.hdf5','/synsets');
%ignore the synset meanings for now unless things become ambiguous
words = extractBefore(synsets,'.'); words = replace(words,'_',' ');

% take all the words in significant regions and make a table summing their
% coefficients
wordcloudtbls = {};
for q = 1:3
    wordcloudtbls{q} = table(words,'VariableNames',{'words'});
    mask = ppredlasso(q,:); tmp = arrayfun(@(grad)fdr(loadings_p(q,grad,:)),1:5,'uniformoutput',false); loadindx = squeeze(any(cat(2,tmp{:}),2));
    mask(loadindx) = fdr(mask(loadindx))<0.05; mask(setdiff(1:length(mask),find(loadindx))) = 0;
    tmp = zscore(gradpredweights(:,q,find(mask)),[],1);
    tmp = squeeze(tmp).*vert(sign(squeeze(median(betaweights(:,:,find(mask)),1)))); % positive means larger absolute value, negative means smaller absolute value
    wordcloudtbls{q}.weights = sum(abs(tmp),2);
    wordcloudtbls{q}.weightspos = sum(tmp.*(tmp>0),2); %sum only positive weights
    wordcloudtbls{q}.weightsneg = sum(abs(tmp.*(tmp<0)),2); %only negative weights
end

% make the figure

figure
p = panel('no-manage-font');
set(gcf,'units','normalized','position',[0 0 1 1])

p.pack('h',{40 60})
p(1).pack('v',{1/3 1/3 1/3}); p(2).pack('v',{1/3 1/3 1/3})

clear ax
for q = 1:3
    mask = ppredlasso(q,:); tmp = arrayfun(@(grad)fdr(loadings_p(q,grad,:)),1:5,'uniformoutput',false); loadindx = squeeze(any(cat(2,tmp{:}),2));
    mask(loadindx) = fdr(mask(loadindx))<0.05; mask(setdiff(1:length(mask),find(loadindx))) = 0;
    [p,~,ax{q},cbar] = plot_volonsurf(rpredlasso(q,:)',shatlas,shinfo,'conte69','mask',mask,'panel',p,...
        'panelindx',{1 q},'colormap',hot,'CLim',[0 max(rpredlasso(q,:).*horz(mask))]);
    if q~=1
        delete(cbar)
    else
        cbar.FontSize = 18; cbar.Label.String = 'Out-of-sample correlation'; cbar.Label.FontSize = 18;
    end
    
    p(2,q).pack('h',{1/2 1/2})

    f = figure;
    wordcloud(wordcloudtbls{q},'words','weightspos','color',palecol(l(2,:)),'HighlightColor','r','title','')
    axw = gca;
    p(2,q,1).select(axw)
    close(f)
    
    f = figure;
    wordcloud(wordcloudtbls{q},'words','weightsneg','color',palecol(l(1,:)),'HighlightColor','b','title','')
    axw = gca;
    p(2,q,2).select(axw)
    close(f)
end
Normalize_Clim([ax{:}],0);

savefig('Canongrads_wordnet_prediction.fig'); export_fig('Canongrads_wordnet_prediction.png','-m4')

% relation to behaviour: PLS

% PLS of second-level components
zgradcomps = zscore(gradcomps,[],1); zneofacts = zscore(neofacts,[],1);
[~,~,XSfacts,YSfacts,~,pctvar_facts,~,stats_facts] = plsregress_perm(zgradcomps,zneofacts,5,10000);
[~,~,XSresps,YSresps,~,pctvar_resps,~,stats_resps] = plsregress_perm(zgradcomps,zscore(neoraw,[],1),5,10000);
stats_facts.pperm = fdr(stats_facts.pperm);
stats_resps.pperm = fdr(stats_resps.pperm);

% ignore version with responses for now since it's not really interpretable
% get loadings for factor version

loads_pls_grads = corr(XSfacts,zgradcomps);
loads_pls_pers = corr(YSfacts,zneofacts);
for i = 1:1000
    bootindx = ceil(rand(size(gradcomps,1),1)*size(gradcomps,1));
    loads_pls_grads_boot(:,:,i) = corr(XSfacts(bootindx,:),zgradcomps(bootindx,:));
    bootindx = ceil(rand(size(gradcomps,1),1)*size(gradcomps,1));
    loads_pls_pers_boot(:,:,i) = corr(YSfacts(bootindx,:),zneofacts(bootindx,:));
end
zloads_pls_grads = loads_pls_grads./std(loads_pls_grads_boot,[],3);
zloads_pls_pers = loads_pls_pers./std(loads_pls_pers_boot,[],3);

for i = 1:5
   ploads_pls_grads(i,:) = fdr(ztop(zloads_pls_grads(i,:)));
   ploads_pls_pers(i,:) = fdr(ztop(zloads_pls_pers(i,:)));
end

% figure for PLS personality results

figure
p = panel('no-manage-font');
set(gcf,'units','normalized','position',[0 0 1 1])
p.pack('v',{1/2 1/2})
p(1).pack('h',{35 45 20})
p(2).pack('h',{35 45 20})

for i = 1:2
    p(i,1).select()
    nicecorrplot(XSfacts(:,i),YSfacts(:,i),{['Gradient component ' num2str(i)],['Personality component ' num2str(i)]})
    FixAxes(gca,16)
    
    p(i,2).select()
    bar(loads_pls_pers(i,:),'facecolor',palecol(l(1,:)),'Edgecolor','none')
    hold on
    mask = ploads_pls_pers(i,:)<0.05;
    bar(find(mask),loads_pls_pers(i,mask),'facecolor','b','edgecolor','none');
    FixAxes(gca,16)
    set(gca,'XTickLabel',{'Agreeableness','Openness','Conscientiousness','Neuroticism','Extraversion'})
    ylabel('Loading')
    errorbar(1:5,loads_pls_pers(i,:),1.96*std(loads_pls_pers_boot(i,:,:),[],3),'LineStyle','none','LineWidth',2,'color','k')
    
    p(i,3).select()
    barh(loads_pls_grads(i,:),'facecolor',palecol(l(2,:)),'Edgecolor','none')
    hold on
    mask = ploads_pls_grads(i,:)<0.05;
    barh(find(mask),loads_pls_grads(i,mask),'facecolor','r','edgecolor','none');
    FixAxes(gca,16)
    xlabel('Loading')
    ylabel('Canonical gradient')
    errorbar(loads_pls_grads(i,:),1:size(loads_pls_grads,2),1.96*std(loads_pls_grads_boot(i,:,:),[],3),'horizontal','LineStyle','none','color','k')
    set(gca,'YDir','reverse')
end

p.de.marginleft = 25; p(1).marginbottom = 25;
p.marginleft = 24; p.marginbottom = 20;

savefig('PLS_personality.fig'); export_fig('PLS_personality.png','-m4')

% supplement: PLS with raw responses

zneoraw = zscore(neoraw,[],1);
resploads_pls_grads = corr(XSresps,zgradcomps);
resploads_pls_pers = corr(YSresps,zneoraw);
for i = 1:1000
    bootindx = ceil(rand(size(gradcomps,1),1)*size(gradcomps,1));
    resploads_pls_grads_boot(:,:,i) = corr(XSresps(bootindx,:),zgradcomps(bootindx,:));
    bootindx = ceil(rand(size(gradcomps,1),1)*size(gradcomps,1));
    resploads_pls_pers_boot(:,:,i) = corr(YSresps(bootindx,:),zneoraw(bootindx,:));
end
respzloads_pls_grads = resploads_pls_grads./std(resploads_pls_grads_boot,[],3);
respzloads_pls_pers = resploads_pls_pers./std(resploads_pls_pers_boot,[],3);

for i = 1:5
   respploads_pls_grads(i,:) = fdr(ztop(respzloads_pls_grads(i,:)));
   respploads_pls_pers(i,:) = fdr(ztop(respzloads_pls_pers(i,:)));
end

% figure for PLS personality results

figure
p = panel('no-manage-font');
set(gcf,'units','normalized','position',[0 0 1 1])
p.pack('v',{1/2 1/2})
p(1).pack('h',{35 45 20})
p(2).pack('h',{35 45 20})

for i = 1:2
    p(i,1).select()
    nicecorrplot(XSresps(:,i),YSresps(:,i),{['Gradient component ' num2str(i)],['Personality component ' num2str(i)]})
    FixAxes(gca,16)
    
    p(i,2).select()
    bar(resploads_pls_pers(i,:),'facecolor',palecol(l(1,:)),'Edgecolor','none')
    hold on
    mask = respploads_pls_pers(i,:)<0.05;
    bar(find(mask),resploads_pls_pers(i,mask),'facecolor','b','edgecolor','none');
    FixAxes(gca,16)
    xlabel('Question')
    %set(gca,'XTickLabel',{'Agreeableness','Openness','Conscientiousness','Neuroticism','Extraversion'})
    ylabel('Loading')
    errorbar(1:60,resploads_pls_pers(i,:),1.96*std(resploads_pls_pers_boot(i,:,:),[],3),'LineStyle','none','LineWidth',2,'color','k')
    
    p(i,3).select()
    barh(resploads_pls_grads(i,:),'facecolor',palecol(l(2,:)),'Edgecolor','none')
    hold on
    mask = respploads_pls_grads(i,:)<0.05;
    barh(find(mask),resploads_pls_grads(i,mask),'facecolor','r','edgecolor','none');
    FixAxes(gca,16)
    xlabel('Loading')
    ylabel('Canonical gradient')
    errorbar(resploads_pls_grads(i,:),1:size(resploads_pls_grads,2),1.96*std(resploads_pls_grads_boot(i,:,:),[],3),'horizontal','LineStyle','none','color','k')
    set(gca,'YDir','reverse')
end

p.de.marginleft = 25; p(1).marginbottom = 25;
p.marginleft = 24; p.marginbottom = 20;

savefig('PLS_personality_resps.fig'); export_fig('PLS_personality_resps.png','-m4')


% working memory

mdl_mem = fitlm(zgradcomps,zscore(behav.ListSort_Unadj));

mem_tbl = mdl_mem.Coefficients;
mem_tbl.Properties.RowNames = cat(1,{'(Intercept)'},cellcat('Gradient ',cellstr(num2str([1:21]')),'',0));
writetable(mem_tbl,'Wrkmem_regression.xlsx','WriteRowNames',true)


% elements for schema

for i = 1:10
subplot(1,10,i)
imagesc(gpca.gradients{75}(:,i))
axis off
end
set(gcf,'color','w')

for i = 1:3
   subplot(3,1,i) 
   histogram(intraclust_corr_comp1{6}{i})
   hold on
   yl = get(gca,'YLim')
    line([median(intraclust_corr_comp1{6}{i}) median(intraclust_corr_comp1{6}{i})],yl,'LineWidth',4,'color',l(1,:),'handlevisibility','off','linestyle','--')
   FixAxes(gca,14)
   title(['Alignment: cluster ' num2str(i)])
   xlabel('Correlation'); ylabel('Counts'); 
   set(gca,'XTick',[],'YTick',[])
end

tmp = cat(3,pca_allalign.gradients{:});
imagesc(squeeze(tmp(:,1,1:10)))
axis off
set(gcf,'color','w')

for i = 1:10
subplot(1,10,i)
imagesc(gradcomps(:,i))
axis off
end
set(gcf,'color','w')

plot_volonsurf(abs(squeeze(loadings(4,1,:)))+0.001,shatlas,shinfo,'conte69',...
    'mask',loadings_p(4,1,:)<0.1,'colormap',hot);

plot_volonsurf(abs(squeeze(loadings(4,2,:)))+0.001,shatlas,shinfo,'conte69',...
    'mask',loadings_p(4,2,:)<0.05,'colormap',hot);

plot_volonsurf(abs(squeeze(loadings(4,3,:)))+0.001,shatlas,shinfo,'conte69',...
    'mask',loadings_p(4,3,:)<0.1,'colormap',hot);

wordcloud(words,exp(rand(859,1)*10),'color',palecol(l(2,:)),'HighlightColor','r','title','')

nicecorrplot(XSresps(:,2),YSresps(:,2),{'Canonical gradients','Personality'})
title('PLS component 1','FontSize',18)
set(gca,'XTick',[],'YTick',[])
