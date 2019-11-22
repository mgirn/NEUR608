% clear all

%% Basic experiment info
%(modify this section)

Subjectlist = 1:168;
% sessionlist = {'Rest' 'Evaluation' 'Acceptance'};
Basedir='/Users/manesh/Desktop/DFCOY/'; %Location of FC data
Outputdir='/Users/manesh/Desktop/DFCOY/Community_detection/'; %Where the parcellation results and ROI output will go

NumIterations=1000; %number of times to run the algorithm
gamma=2; %gamma is the resolution parameter that infuences the number of communities that will be detected.
gamma_agreement=2;
% roi_dir='C:\Users\mattd\Documents\Neuroscience_software\ROIs\Yeo\Cortical\Yeo400\'; %location of Yeo ROIs
% numROIs=57; % Number of ROIs
%%
% k=1;%rest
file=[Basedir 'dfcoy_mancovan_results_fnc.mat'];
load(file)
fnc_corrs=icatb_vec2mat(fnc_corrs);
W=fnc_corrs;
% reshape(mean(FC_allsubs,1), 100, 100);
%     Z2=Z(Subjectlist,1:numROIs,1:numROIs);
% W_young=W(young_indices,:,:);
 W_old=W(old_indices,:,:);
W_old=mean(W_old,1); %mean matrix across subjects
W_old=reshape(W_old,45,45);
% W_young(W_young<0)=0;%Remove negative weights
idx = isnan(W); if any(any(idx)); W(idx)=0; end %Remove NaN self-connections


%Iterate community detection algorithm
for iteration=1:NumIterations
    Q0 = -1; Q = 0; % initialize modularity values
    while Q-Q0>1e-5; % while modularity increases
        Q0 = Q;
        [Ci Q]=community_louvain(W_old,gamma);
    end %while
    Ci=Ci';
    communities(iteration,:)=Ci; %Store community assignment vectors for each iteration
    
    %% Create co-classification matrix (probability that pair of nodes belong
    %to same module)
    for SourceNode=1:length(Ci)
        for node=1:length(Ci)
            if Ci(1,SourceNode)==Ci(1,node) %Determine if each pair of nodes are in same community and assign a value of 1 or 0
                Agreement_matrix(SourceNode,node)=1;
            else
                Agreement_matrix(SourceNode,node)=0;
            end
        end
    end
    
    Agreement_matrix(1:length(Agreement_matrix)+1:end)=0; %make diagonal 0
    Agreement_Matrices(:,:,iteration)=Agreement_matrix; %Store values for each iteration
end %iterations

Summed_Agreement_Matrix=sum(Agreement_Matrices,3); %add matrices
Probability_Matrix=Summed_Agreement_Matrix./NumIterations; %determine probability (well actually proportion) of co-classification

%determine community assignments on agreement matrix
Q0 = -1; Q_individual = 0; % initialize modularity values
while Q_individual-Q0>1e-5; % while modularity increases
    Q0 = Q_individual;
    [Ci_Group_Optimal Q_Group_Optimal]=community_louvain(Probability_Matrix,gamma_agreement);
end %while


cd(Outputdir)
save('Community_detection_results.mat','Ci_Group_Optimal','Q_Group_Optimal', 'communities', 'gamma_agreement', 'gamma');


%
%% Create network masks

% cd(roi_dir);
%
% for k = 1:max(Ci_Group_Optimal)
%     current_index=find(Ci_Group_Optimal==k); %set of ROI indices for currrent module
%     ROI_Name=['Network' int2str(k)];
%     clear roilist network_holder
%     for n=1:length(current_index)
%         network_holder{n}=['Yeo_400_' int2str(current_index(n)) '_roi.mat'];  %Grab correspodning set of ROIs
%     end
%
%     roilist=network_holder;
%
%
%     % Combine ROIs
%     clear func j m c a
%     roilist=maroi('load_cell',roilist); % roilist contains the ROIs to be combined
%     [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Combine ROIs');
%     %func = ('r1 | r2 | r3'); %Specify function for combining ROIs. The '|' indicates that ROIs are to be unified.
%
%     for j=1:length(roilist)
%         a{j,:}=['r' int2str(j) '|'];
%     end
%
%     c=[];
%     for m=1:size(a,1)
%         c=[c a{m,:}];
%     end
%     func=c(1:end-1);
%
%
%     for v = 1:length(roilist)
%         eval(sprintf('r%d = roilist{%d};', v, v));
%     end
%     %try
%     %eval(['o=' func ';']);
%     o=eval(func);
%     %catch
%     %   warning(['Hmm, probem with function ' func ': ' lasterr]);
%     %  return
%     %end
%
%
%     % if isa(o, 'maroi')
%     o = label(o, [ROI_Name]); % P contains subject IDs
%     saveroi(o, fullfile(Outputdir, strcat(ROI_Name,'_roi.mat')));  % save as .mat
%     save_as_image(o, fullfile(Outputdir, strcat(ROI_Name, '_roi.nii'))); % Save as image (.nii)
%     %end
%
% end
%
% %}