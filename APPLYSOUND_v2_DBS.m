function [OUTEEG,TRIALNOISE]=APPLYSOUND_v2_DBSEEG,chanlocs,selecttrials,trialwise,iter,ignorechans)
% AN EXAMPLE MATLAB SCRIPT FOR USING DDWIENER AND SOUND 
%
% This is an EEGLAB-compatible MATLAB script that demonstrates how DDWiener
% and SOUND can be used to clean a partially noisy EEG dataset.
%
% The sample EEG dataset used in this script can be downloaded from the HeadIT
% Data Repository (http://headit.ucsd.edu/) maintained by the Swartz Center
% for Computational Neuroscience at the University of California at San 
% Diego. The used dataset corresponds to the study of Auditory Two-Choice
% Response Task with an Ignored Feature Difference, session 7, recording 1 and can be 
% downloaded from http://headit.ucsd.edu/recording_sessions/99bc255c-a238-11e2-b5e7-0050563f2612.
% In order to use the HeadIT Data Repository, the user must agree to the
% HeadIT Data Repository Terms of use and the Data use agreement.
%
% Since the sample EEG dataset is in Biosemi data format (.bdf), the BioSig
% toolbox is also needed and can be downloaded from http://biosig.sourceforge.net/download.html.
%
% Before running this script, make sure that the script and the functions 
% used in it, EEGLAB toolbox (Version 14.0.0b has been tested), BioSig
% toolbox (Version 3.1.0 has been tested), the used EEG dataset, and the 
% file containing the channel locations are located in your MATLAB path.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, Johanna Metsomaa, Sara Liljander, and Risto Ilmoniemi. 
% Department of Neuroscience and Biomedical Engineering (NBE), School of Science, Aalto University  
% .........................................................................


if ~exist('chanlocs','var')||isempty(chanlocs)
    load /data/rcho/TripolarEEG/DBSstudychanlocs.mat
end
% Set the number of iterations when estimating the channel-specific noise
% estimates. Iter = 5 was found to be sufficient in the datasets studied in
% the original article. 
if ~exist('iter','var')||isempty(iter)
    iter = 10;
end
if ~exist('selecttrials','var')||isempty(selecttrials)
elseif islogical(selecttrials)
    selecttrials=find(selecttrials);     
    EEG=pop_select(EEG,'trial',selecttrials);
elseif ~islogical(selecttrials)
    [trtype,trcount]=occ(selecttrials);
    if median(trcount)>1
        %this routine was added to try to handle passing a condition vector
        %for select trials. The data is returned in the same shape as
        %EEG.data, but the SOUND correction is applied to each trial type
        %individually.
        OUTEEG=zeros(length(chanlocs),size(EEG.data,2),size(EEG.data,3));
        disp('Applying SOUND by conditions given in ''selecttrials''.')
        for i=1:numel(trtype)
            disp(['    Condition==' num2str(trtype(i)) ' with ' num2str(trcount(i)) ' trials.'])
            if trcount(i)>1
                TMP=APPLYSOUND_v2(EEG,chanlocs,selecttrials==trtype(i),trialwise,iter);
            else
                TMP=APPLYSOUND_v2(EEG,chanlocs,selecttrials==trtype(i),0,iter);
            end
            OUTEEG(:,:,selecttrials==trtype(i))=TMP;
        end
        TRIALNOISE=[];
        return
    else
        EEG=pop_select(EEG,'trial',selecttrials);
    end
end
if ~exist('trialwise','var')||isempty(trialwise); trialwise=false; end

endchans=1:length(chanlocs);
keepchans=1:size(EEG.data,1);
missingchans=ismember({chanlocs.labels},{EEG.chanlocs.labels});
refchan=find(mean(mean(EEG.data(:,:,:),3),2)==0);

if ~exist('ignorechans','var')||isempty(ignorechans)
    keepchans(refchan)=[];
else
    keepchans([ignorechans refchan])=[];
end


for i=1:size(chanlocs,2)-1
    chanlocs(i).type='EEG';
end
chanlocs(end).type='EEG';
EEG.chanlocs=chanlocs;
[LFM_sphere] = construct_spherical_lead_field(chanlocs); % full leadfield matrix, (chansXdipoles)
[LFM_sphere_mean] = ref_ave(LFM_sphere); % re-referenced leadfield matrix, for reprojection to channels (this will also project dipole estimates to missing channels)
LFM_sphere_ref = LFM_sphere(missingchans,:); % this is the original leadfield matrix, but reduced down to the input size of the input EEG structure
% quick note about these transformations. The LFM_sphere_ref and LFM_sphere
% matrices will be the same if no channels have been removed from the EEG
% structure prior to sound. If a channel has been removed, or if another
% one has been specified in ignorechans, the computations are completed on
% the existing data, and reprojected to the full set of channels specified
% in chanlocs. This weird roundabout approach is to allow multistage
% channel removals


% Set the SOUND parameters:

% Set the reqularization parameter of the minimum-norm estimate when
% finding the noise-free current estimates. This can be adjusted if clear
% under- or overcorrection is observed. lambda_value = 0.1 was used in the 
% original article. 
lambda_value = 0.1; 


if ~trialwise
    % 4. Taking noisy trials into account with DDWiener
    [EEG_evo,est_noise,tris,scals] = trial_noise_estimator(EEG.data); 
    TRIALNOISE=struct('EEV_evo',EEG_evo,'est_noise',est_noise,'tris',tris,'scals',scals);

    % 5. Use SOUND to clean the channel-specific noise
    % Build the spherical lead field, using the theoretical electrode-locations
    % of the data.

    % Re-reference the data and the lead field to the channel with the least noise
    [~, sigmas] = DDWiener(EEG_evo);
    [~,bestC] = min(sigmas);
    [datatmp] = ref_best(EEG_evo, bestC);
    [LFM_sphere_tmp] = ref_best(LFM_sphere, bestC);

    % Run the SOUND algorithm:

    % (The data contain several channels so SOUND may take a few 
    % minutes.)
    chans = setdiff(1:size(EEG_evo,1),bestC);
    [corrected_data,x,sigmas,dn] = SOUND(datatmp(chans,:), LFM_sphere_tmp(chans,:), iter, lambda_value); 
    TRIALNOISE.sigmas=sigmas;
    TRIALNOISE.dn=dn;
    TRIALNOISE.corrected_data=corrected_data;

    % Re-reference the data and the lead field to the channel average:
    [LFM_sphere_mean] = ref_ave(LFM_sphere);
    OUTEEG = LFM_sphere_mean*x;

elseif trialwise==1
    
    TRIALNOISE=[];
    OUTEEG=EEG.data;

    parfor tr=1:size(EEG.data,3)
    disp(['trial ' int2str(tr)])
    EEGtmp=OUT(:,:,tr);
    [~, sigmas] = DDWiener(EEGtmp);
    [~,bestC] = min(sigmas);
    [datatmp] = ref_best(EEGtmp, bestC);
    [LFM_sphere_tmp] = ref_best(LFM_sphere, bestC); 
    chans = setdiff(1:size(EEGtmp,1),bestC);
    [~,x,~,~] = SOUND(datatmp(chans,:), LFM_sphere_tmp(chans,:), iter, lambda_value); 
    [LFM_sphere_mean] = ref_ave(LFM_sphere);
    tmp=LFM_sphere_mean*x;  
    OUTEEG(:,:,tr) = tmp;
    end
    
elseif trialwise==2
    
    [EEG_evo,est_noise,tris,scals] = trial_noise_estimator(EEG.data(keepchans,:,:)); 
    TRIALNOISE=struct('EEV_evo',EEG_evo,'est_noise',est_noise,'tris',tris,'scals',scals);
    
    % Re-reference the data and the lead field to the channel with the least noise
    [~, sigmas] = DDWiener(EEG_evo);
    [~,bestC] = min(sigmas);
    if ~isempty(refchan); EEG_evo(end+1,:)=0; keepchans(end+1)=keepchans(end)+1; end
    [datatmp] = ref_best(EEG_evo, bestC);
    [LFM_sphere_tmp] = ref_best(LFM_sphere_ref(keepchans,:), bestC);
    
    chans = setdiff(1:size(EEG_evo,1),bestC);
    % Run the SOUND algorithm:
    [~, x, sigmas] = SOUND(datatmp(chans,:), LFM_sphere_tmp(chans,:), iter, lambda_value);
    EEG_tmp_correct = LFM_sphere_mean*x;

    EEG_tmp_trials = zeros(size(LFM_sphere,1),size(EEG.data,2),size(EEG.data,3));
    % (The data contain several channels so SOUND may take a few 
    % minutes.)
     for k = 1:size(EEG.data,3)
        [datatmp] = ref_best(EEG.data(keepchans,:,k), bestC);
        [corrected_data, x] = correct_with_known_noise(datatmp(setdiff(1:size(datatmp,1),bestC),:),...
            LFM_sphere_tmp(setdiff(1:size(datatmp,1),bestC),:), lambda_value,  sigmas);
        EEG_tmp_trials(:,:,k) = LFM_sphere_mean*x;
    end
        OUTEEG=EEG_tmp_trials;
end

% OUTEEG=EEG;
% OUTEEG.data=OUT;
% OUTEEG.chanlocs=chanlocs;
% OUTEEG.nbchan=length(chanlocs);

end



function [corrected_data, x] = correct_with_known_noise(data, LFM,lambda0,  sigmas)
% This function corrects a data segment when the noise distribution is already known.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university  
% .........................................................................
    chanN = length(sigmas);
    W = diag(1./sigmas);
    WL = W*LFM;
    clear LFM;
    WLLW = WL*WL';
    x = WL'*((WLLW + lambda0*trace(WLLW)/chanN*eye(chanN))\(W*data));
    corrected_data = [];

end

function [types,counts]=occ(vector)
%This function parses an input vector to determine if it is a location
%vector (i.e. only one occurance per number; e.g. [1 2 3 4 7 8 9 11 ... ]) 
%or a condition vector (i.e. multiple occurances per number; 
%e.g. [1 1 1 2 2 3 3 3 ... ]).

types=unique(vector); 
if isrow(types); types=types'; end

counts=types;
for i=1:length(types)
    counts(i)=sum(vector==types(i));
end
    
end