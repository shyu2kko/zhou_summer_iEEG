%% mapping electrodes in fsLR5k
% copied from Jess (thank goodness)

addpath('/data_/mica1/02_codes/NIfTI_20140122')
addpath('/data_/mica1/02_codes/matlab_toolboxes/gifti-1.6/');
addpath('/host/fladgate/local_raid/jessica/MICs/scripts/');
addpath('/data_/mica1/02_codes/NIfTI_20140122')
addpath('/data_/mica1/02_codes/matlab_toolboxes/gifti-1.6/');


%-------------------setting up--------------%
bids = '/data_/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0/';
ieeg_bids = '/host/oncilla/local_raid/iEEG-BIDS';
savedir = '/mica/mica1/03_projects/yigu/summer_iEEG/data/';

subjList = ['PX001'; 'PX005'; 'PX007'; 'PX009'; ...
    'PX010'; 'PX012'; 'PX015'; 'PX019'; 'PX023'; ...
    'PX028'; 'PX029'; 'PX031'; 'PX034'; 'PX040'; ...
    'PX041'; 'PX045'; 'PX050'; 'PX051'; 'PX053'; ...
    'PX065']; % omit PX004 because previous implant


%------------------define spherical patch-------------%
nFaces = 49; r = 5; % 5mm radius
%%
% 1. GET CLOSEST VERTEX TO CHANNEL
% initiate data structure; I won't change this
vertex_match = {};
electrode_name_list = {};

% concatenate info from all implanted patients
parfor ii = 1:size(subjList,1)

    % subject
    this_subj = subjList(ii,:); disp(this_subj);

    % paths & filenames
    anat = [ieeg_bids, '/sub-', this_subj, '/anat/nativepro/'];
    surf = [bids, '/sub-', this_subj, '/ses-01/surf/'];
    fnNativepro = ['sub-', this_subj, '_ses-01_space-nativepro_T1w_brain.nii.gz'];
    fnRibbon = ['sub-', this_subj, '_ses-01_space-nativepro_ribbon.nii.gz'];
    fnSurfLh = ['sub-', this_subj, '_ses-01_hemi-L_space-nativepro_surf-fsLR-5k_label-midthickness.surf.gii'];
    fnSurfRh = ['sub-', this_subj, '_ses-01_hemi-R_space-nativepro_surf-fsLR-5k_label-midthickness.surf.gii'];
    electrodeNames = dir([anat,'electrode*']);
    
    % Load scans: need T1w in nativepro space and nativepro ribbon
    nativeprohdr = load_nii([anat,fnNativepro]);nativepro = nativeprohdr.img;
    ribbonhdr = load_nii([anat,fnRibbon]);ribbon = ribbonhdr.img;
    
    % Scan info
    res = nativeprohdr.hdr.dime.pixdim(2:4);
    xOffset = nativeprohdr.hdr.hist.qoffset_x;
    yOffset = nativeprohdr.hdr.hist.qoffset_y;
    zOffset = nativeprohdr.hdr.hist.qoffset_z;
    
    % Load surfaces
    lh_tmp = gifti([surf, fnSurfLh]); lh = struct(); lh.coord = lh_tmp.vertices'; lh.tri = lh_tmp.faces';
    rh_tmp = gifti([surf, fnSurfRh]); rh = struct(); rh.coord = rh_tmp.vertices'; rh.tri = rh_tmp.faces';
    
    % load channel list from timeseries
    px_data = load(['/host/oncilla/local_raid/jessica/microstructure_iEEG/data/ieeg-mni/sub-',this_subj,'_clean_ts.mat'], 'contacts_load', 'dataLabels');
    contacts_names_ts = px_data.contacts_load';
    channel_names_ts = px_data.dataLabels';
    
    % Load electrodes contacts and map to GM
    vertex_match_electrodes = {};
    
    % loop through electrodes
    for jj = 1:size(electrodeNames,1)
        
        % Full name of nii file
        this_name = electrodeNames(jj).name;
        electrode_name_list{ii}{jj} = this_name;
        
        % Name of electrode
        first_split = strsplit(this_name,'_');
        second_split = strsplit(first_split{2},'.');
        shortName = second_split{1};
        
        % Load electrode nifti file
        disp(['Loading electrode ', char(string(jj)), '/', char(string(size(electrodeNames,1))), ': ', this_name])
        tmp = load_nii([anat,electrodeNames(jj).name]);
        this_electrode = tmp.img;
        
        ncontacts = length(unique(this_electrode))-1; % exclude 0
        vertex_match_contacts = zeros(ncontacts,5); % contact1, contact2, vertex match, idx contact, idx channel (in ts)
        
        % Loop through contacts of electrode
        for contact = 1:ncontacts-1 % exclude last contact because bipolar channels
            
            % Compare contact name with names in edf file - some
            % contacts are missing from the recording so we can skip
            % those 
            pat = [lower(shortName),char(string(contact)),'+$'];
            for labels = 1:length(contacts_names_ts)
                this_contact = lower(string(contacts_names_ts(labels)));
                if ~isempty(regexp(this_contact,pat))
                    vertex_match_contacts(contact,4) = labels;
                end
            end
            if vertex_match_contacts(contact,4) == 0
                disp(['contact # ', char(string(contact)), ' not in timeseries... moving to next contact']);
                vertex_match_contacts(contact,4) = NaN;
                continue
            end
            pat = ['\w*',lower(shortName),char(string(contact)),'-\w*'];
            for labels = 1:length(channel_names_ts)
                this_channel = lower(string(channel_names_ts(labels)));
                if ~isempty(regexp(this_channel,pat,'match'))
                    vertex_match_contacts(contact,5) = labels;
                end
            end
            
            % contact numbers constituting bipolar channel
            vertex_match_contacts(contact,1) = contact;
            vertex_match_contacts(contact,2) = contact+1;
            
            % Check if first contact is in gray matter
            c1 = zeros(size(this_electrode));
            c1(this_electrode == contact) = 1;
            
            statsc1 = regionprops(c1, 'centroid'); % centroid of the contact
            statsc1 = [statsc1.Centroid(2), statsc1.Centroid(1), statsc1.Centroid(3)]-1; % fix weird flip x-y
            
            [X,Y,Z] = sphere(nFaces); % Define sphere with 5mm radius
            X = X*(r*res(1)); Y = Y*(r*res(2)); Z = Z*(r*res(3));
            X2 = X + statsc1(1); Y2 = Y + statsc1(2); Z2 = Z + statsc1(3);
            c1_sphere = zeros(size(this_electrode));
            c1_sphere(round(X2(:)), round(Y2(:)), round(Z2(:))) = 1;
            
            overlap = double(ribbon) + (double(c1_sphere)); % check sphere overlap with GM
            is_overlap = find(overlap == 2);
            if (length(is_overlap) / length(find(c1_sphere==1))) < 0.05
                disp(['contact # ', char(string(contact)), ' not in gray matter... moving to next contact']);
                vertex_match_contacts(contact,3) = NaN;
                
                % if no overlap of sphere with GM, continue loop to next contact
                continue
            end
            
            % add second contact of channel to volume
            c1(this_electrode == contact+1) = 1;
            
            % centroid of the channel, made up of 2 contacts
            statsc1 = regionprops(c1, 'centroid');
            statsc1 = [statsc1.Centroid(2), statsc1.Centroid(1), statsc1.Centroid(3)]-1; % fix weird flip x-y
            
            % get centroid x y z in mm
            mid_mm = [statsc1(1) * res(1) + xOffset, ...
                statsc1(2) * res(2) + yOffset, ...
                statsc1(3) * res(3) + zOffset];
            
            % Find closest vertex
            if strcmpi(this_name(11), 'R')
                this_hemi = rh.coord;
                add_hemi = size(this_hemi,2);
            else
                this_hemi = lh.coord;
                add_hemi = 0;
            end
            D = pdist2(mid_mm,this_hemi','euclidean');
            vertex = find(D == min(D));
            vertex_match_contacts(contact,3) = vertex + add_hemi;
        end
        vertex_match_electrodes{jj} = vertex_match_contacts;
    end
    vertex_match{ii} = vertex_match_electrodes;
end
save([savedir, 'fsLR-5k_vertex_match.mat'], 'vertex_match', 'electrode_name_list')

% clean up
vertex_match_clean = {};
for ii = 1:size(subjList,1)
tmp = vertex_match{ii};
for electrode = 1:length(tmp)
    tmpelectrode = tmp{electrode};
    tmpelectrode(tmpelectrode(:,3) == 0 | isnan(tmpelectrode(:,3)) | tmpelectrode(:,5) == 0,:) = [];
    tmp{electrode} = tmpelectrode;
end
vertex_match_clean{ii} = tmp;
end
save([savedir, 'fsLR-5k_vertex_match_clean.mat'], 'vertex_match_clean', 'electrode_name_list')

%done ! don't redo.

%% Generate channel patches: 5mm region on surface
load([savedir, 'fsLR-5k_vertex_match_clean.mat'], 'vertex_match_clean', 'electrode_name_list')
% Save txt files with coordinates and faces to compute geodesic distance in python
% NOTE: in output file 'PX000_hemi-X_mesh.txt', you will have to delete all spaces, and replace all commas by spaces afterwards to run pygeodesic...

hemi = ['L'; 'R'];
surfName = '/data/mica1/01_programs/micapipe-v0.2.0/surfaces/fsLR-5k.R.inflated.surf.gii';
surf = gifti(surfName); nFS = length(surf.vertices');
parfor ii = 1:size(subjList,1)
    
    % subject
    this_subj = subjList(ii,:);
    disp(this_subj);
    
    surfDir = [bids, '/sub-', this_subj, '/ses-01/surf/'];
    
    % output mesh
    for hh = 1:2
        surfName = ['sub-', this_subj, '_ses-01_hemi-', hemi(hh),'_space-nativepro_surf-fsLR-5k_label-midthickness.surf.gii'];
        surf = gifti([surfDir, surfName])
        info = [length(surf.vertices), length(surf.faces)];
        
        dlmwrite([savedir,'/geodist/', this_subj,'_hemi-',hemi(hh),'_mesh0.txt'], info, 'delimiter', ',', 'roffset',0, 'coffset', 0, 'precision', '%7.0f');
        dlmwrite([savedir,'/geodist/',this_subj,'_hemi-',hemi(hh),'_mesh0.txt'], surf.vertices, 'delimiter', ',', 'roffset',0, 'coffset', 0, 'precision', '%7.5f', '-append');
        dlmwrite([savedir,'/geodist/',this_subj,'_hemi-',hemi(hh),'_mesh0.txt'], surf.faces-1, 'delimiter', ',', 'roffset',0, 'coffset', 0, 'precision', '%7.0f', '-append');
    end
end

%% mesh file from fsLR 32k

savedir = '/mica/mica1/03_projects/yigu/surfaces/';

surfName = '/data/mica1/01_programs/micapipe-v0.2.0/surfaces/fsLR-32k.R.inflated.surf.gii';
surf = gifti(surfName); nFS = length(surf.vertices');
info = [length(surf.vertices), length(surf.faces)];
dlmwrite([savedir,'/geodist/fsLR-32k_hemi-R_mesh0.txt'], info, 'delimiter', ',', 'roffset',0, 'coffset', 0, 'precision', '%7.0f');
dlmwrite([savedir,'/geodist/fsLR-32k_hemi-R_mesh0.txt'], surf.vertices, 'delimiter', ',', 'roffset',0, 'coffset', 0, 'precision', '%7.5f', '-append');
dlmwrite([savedir,'/geodist/fsLR-32k_hemi-R_mesh0.txt'], surf.faces-1, 'delimiter', ',', 'roffset',0, 'coffset', 0, 'precision', '%7.0f', '-append');
