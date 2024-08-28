%% make electrodes plot by region
% adapted from Casey Paquola

clear all

%------set up------%
outDir = '/mica/mica1/03_projects/yigu/summer_iEEG/data/';

% casey's paths
addpath(genpath('/data_/mica1/03_projects/casey/micasoft/matlab/afni_matlab'))
addpath('/data_/mica1/03_projects/casey/micasoft/matlab/useful')
addpath('/data_/mica1/03_projects/casey/sandbox1/deep_embedding/scripts')
addpath('/data_/mica1/03_projects/casey//sandbox1/fieldtrip/')
addpath(genpath('/data_/mica1/03_projects/casey/micaopen/surfstat'))

% for plotting
load(['/mica/mica1/03_projects/yigu/colorbrewer.mat']);
load('/mica/mica1/03_projects/yigu/summer_iEEG/data/fsLR-5k_vertex_match_clean.mat', 'vertex_match_clean');

% make region-wise colors
color_palette = colorbrewer(1).qual.Set2;
% ids - dropping PX004 because scan post-implantation
mr_id = {'PX001', 'PX005', 'PX007', 'PX009', ...
    'PX010', 'PX012', 'PX015', 'PX019', 'PX023', ...
    'PX028', 'PX029', 'PX031', 'PX034', 'PX040', ...
    'PX041', 'PX045', 'PX050', 'PX051', 'PX053', ...
    'PX065'};

%%

% load the surfaces of interest, need to change into .obj for plotting fsLR-5k
fsDir       = '/mica/mica1/03_projects/yigu/surfaces';
f5k_lh = SurfStatReadSurf([fsDir '/fsLR-5k/surf/fsLR-5k.L.surf']);
f5k_rh = SurfStatReadSurf([fsDir '/fsLR-5k/surf/fsLR-5k.R.surf']);

% testing that the surfaces loaded properly
[n, nvtx] = size(f5k_lh.coord);
%tmp = ones([1, nvtx]);
%SurfStatView(tmp, f5k_lh)
%
    
all_vert = [0;0];

for s=1:20
    % count how many electrodes
    nch = size(vertex_match_clean{s}, 2);
    for ch=1:nch
        idx = vertex_match_clean{s}{ch}(:,3); % all the centroid vertices for each channel
        if ~isempty(idx) % append
            all_vert = [all_vert [ones(1,length(idx))*s; idx']];
        end
    end
end
    
% cleaning up - remove 0's and nans
all_vert(:,1) = [];
all_vert(:,isnan(all_vert(2,:))) = [];
all_vert(:,all_vert(2,:)==0) = [];

% separate hemispheres
left = all_vert(2,:) <= nvtx;
lh_vert = all_vert(:,left);
right = all_vert(2,:) > nvtx;
rh_vert = all_vert(:, right);

% mirror rh
rh_vert(2,:) = rh_vert(2,:) - nvtx;

FS_lh = f5k_lh;
FS_rh = f5k_rh;

% make this the first row of the plotting xx_vert
regions_assignment = randi([1 17], 1, size(lh_vert, 2));

% extend the color palette

cmap = vertcat(colorbrewer.qual.Set2{8}, colorbrewer.qual.Set3{1,12});

f = figure('units','centimeters','outerposition',[0 0 30 30]);
    
     a(1) = axes('position', [0.41 0.5 0.2 0.2]);
     trisurf(FS_lh.tri, FS_lh.coord(1,:), FS_lh.coord(2,:), FS_lh.coord(3,:), ones(1,length(FS_lh.coord)))
     view(-90,0); 
     daspect([1 1 1]); axis tight; camlight; axis vis3d off;
     lighting phong; material dull; shading interp;
     colormap(white)
     hold on
     scatter3(FS_lh.coord(1,lh_vert(2,:)), FS_lh.coord(2,lh_vert(2,:)), ...
     FS_lh.coord(3,lh_vert(2,:)), 20, cmap(lh_vert(1,:),:)/255 , ...
     'filled', 'MarkerEdgeColor', 'k')
     


     % left hemisphere medial
     a(2) = axes('position', [0.41 0.33 0.2 0.2]);
     trisurf(FS_lh.tri, FS_lh.coord(1,:), FS_lh.coord(2,:), FS_lh.coord(3,:), ones(1,length(FS_lh.coord)))
     view(90,0); 
     daspect([1 1 1]); axis tight; camlight; axis vis3d off;
     lighting phong; material dull; shading interp;
     colormap(white);
     hold on
     scatter3(FS_lh.coord(1,lh_vert(2,:)), FS_lh.coord(2,lh_vert(2,:)), ...
     FS_rh.coord(3,lh_vert(2,:)), 20, cmap(lh_vert(1,:),:)/255 , ...
     'filled', 'MarkerEdgeColor', 'k')
    
    % right hemisphere lateral
    a(3) = axes('position', [0.7 0.33 0.2 0.2]);
    trisurf(FS_rh.tri, FS_rh.coord(1,:), FS_rh.coord(2,:), FS_rh.coord(3,:), ones(1,length(FS_rh.coord)))
    view(-90,0); 
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading interp;
    colormap(white);
    hold on
    scatter3(FS_rh.coord(1,rh_vert(2,:)), FS_rh.coord(2,rh_vert(2,:)), ...
        FS_rh.coord(3,rh_vert(2,:)), 20, cmap(rh_vert(1,:),:)/255 , ...
        'filled', 'MarkerEdgeColor', 'k')
    
    
    % right hemisphere medial
    a(4) = axes('position', [0.7 0.5 0.2 0.2]);
    trisurf(FS_rh.tri, FS_rh.coord(1,:), FS_rh.coord(2,:), FS_rh.coord(3,:), ones(1,length(FS_rh.coord)))
    view(90,0); 
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading interp;
    colormap(white);
    hold on
    scatter3(FS_rh.coord(1,rh_vert(2,:)), FS_rh.coord(2,rh_vert(2,:)), ...
        FS_rh.coord(3,rh_vert(2,:)), 20, cmap(rh_vert(1,:),:)/255 , ...
        'filled', 'MarkerEdgeColor', 'k')

    hold on
    for k = 1:20
        H(k) = scatter3(1,1,1,1, cmap(k,:)/255, 'filled');
    end
    legend(H, mr_id)
%%
fname = [outDir 'iEEG-BIDS_n20_cortical_electrodes_all_test.png'];
exportfigbo(f, char(fname),'png', 10);    

%% plot hippocampal

% load surface of interest, changed to .obj for plotting
fsDir       = '/mica/mica1/03_projects/yigu/surfaces';
hipp_lh = SurfStatReadSurf([fsDir '/canonical_surf/surf/hipp.surf']);
hipp_rh = SurfStatReadSurf([fsDir '/canonical_surf/surf/hipp.surf']);

[n, nvtx] = size(hipp_lh.coord);
%tmp = ones([1, nvtx]);
%SurfStatView(tmp, hipp_lh)

% read vertex and patient id text
lh_vert = readtable('/mica/mica1/03_projects/yigu/summer_iEEG/data/lh_hipp_vert.txt', 'ReadVariableNames',false);
lh_vert = table2array(lh_vert);
rh_vert = readtable('/mica/mica1/03_projects/yigu/summer_iEEG/data/rh_hipp_vert.txt', 'ReadVariableNames', false);
rh_vert = table2array(rh_vert);

cmap = vertcat(colorbrewer.qual.Set2{8}, colorbrewer.qual.Set3{1,12});

FS_lh = hipp_lh; FS_rh = hipp_rh;
f = figure('units','centimeters','outerposition',[0 0 30 30]);
    
    a(1) = axes('position', [0.41 0.5 0.2 0.2]);
         trisurf(FS_lh.tri, FS_lh.coord(1,:), FS_lh.coord(2,:), FS_lh.coord(3,:), ones(1,length(FS_lh.coord)))
         view(-45,0); 
         daspect([1 1 1]); axis tight; camlight; axis vis3d off;
         lighting phong; material dull; shading interp;
         colormap(white)
         hold on
         scatter3(FS_lh.coord(1,lh_vert(2,:)), FS_lh.coord(2,lh_vert(2,:)), ...
         FS_lh.coord(3,lh_vert(2,:)), 20, cmap(lh_vert(1,:),:)/255 , ...
         'filled', 'MarkerEdgeColor', 'k')
    
    
    a(2) = axes('position', [0.41 0.33 0.2 0.2]);
         trisurf(FS_rh.tri, FS_rh.coord(1,:), FS_rh.coord(2,:), FS_rh.coord(3,:), ones(1,length(FS_rh.coord)))
         view(-45,0); 
         daspect([1 1 1]); axis tight; camlight; axis vis3d off;
         lighting phong; material dull; shading interp;
         colormap(white);
         hold on
         scatter3(FS_rh.coord(1,rh_vert(2,:)), FS_rh.coord(2,rh_vert(2,:)), ...
         FS_rh.coord(3,rh_vert(2,:)), 20, cmap(rh_vert(1,:),:)/255 , ...
         'filled', 'MarkerEdgeColor', 'k')
fname = [outDir 'iEEG-BIDS_n20_hippocampal_electrodes_all_test.png'];
exportfigbo(f, char(fname),'png', 10);   











