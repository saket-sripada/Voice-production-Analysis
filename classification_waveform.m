% waveform classifier
% useful command -> lookfor audio
%% read files
% generates complete file name based on serial number of file
 
sBatch = 'ug00'; % common starting string of files names
ftype = '.wav'; %file format
n_subjects = 250; % No. of sujects in study
fnames=[];
%AmPeaks =[];
%tInt = [];
for i=1:n_subjects
    if i<10
        name = strcat(sBatch,'00',num2str(i),ftype);
    elseif i<100
        name = strcat(sBatch,'0',num2str(i),ftype);
    elseif i<1000
        name = strcat(sBatch,num2str(i),ftype);
    end
    fnames = [fnames ; name];
end

%% read and store data
pf = 0.4; % prominence factor for local maxima
wavData = [];
for i=1:size(fnames,1)
    [d,f] = audioread(fnames(i,:));       % read wav file using audioread and store data nd fs
    [pks, locs] = findpeaks(d,pf);        % calculate local maxima and their time intervals/occurances 
    wavData(i).name = d;
    wavData(i).data = d;                  % stores the data from the wave file for convenience  
    wavData(i).fs =f;                     % stores fs value into structure
    wavData(i).Amplitude_Peaks = pks;             % stores vector of amplitude of peaks into structure 
    wavData(i).time_periods = locs;              % stores timepoints at which maxima occur into structure
    % calculates P1 for amplitude of local maxima 
    wavData(i).P1amp = std(wavData(i).Amplitude_Peaks)/mean(wavData(i).Amplitude_Peaks);  
    % calculates P1 for time vactor
    wavData(i).P1_period = std(wavData(i).time_periods)/mean(wavData(i).time_periods); 
end
%% calculate purtubation parametrs
yt=[]; zt=[];
ya=[]; za=[];

for k = 1:size(fnames,1)
    N = size(wavData(k).Amplitude_Peaks,1);
    for i = 1:10
        for j =1:N-i
            ya(i) = wavData(k).Amplitude_Peaks(i+j)-wavData(k).Amplitude_Peaks(j); % calculate x(j+i) - x(j)
            % calculate z(i) for amplitude vector of local maxima
            za(i) = (abs(ya(i))/(N-i))/mean(wavData(k).Amplitude_Peaks);     
            yt(i) = wavData(k).time_periods(i+j)-wavData(k).time_periods(j);
            % calculate z(i) for time vector
            zt(i) = (abs(yt(i))/(N-i))/mean(wavData(k).time_periods);
        end
    end
    wavData(k).P2amp = min(za); % calculate P2 for amplitude vector
    wavData(k).P2_period = min(zt); % calculate P2 for time vector
end
%% exporting from structure to table
X_trn = [];
% for i=1:size(fnames,1)
%     X = [X wavData(i)]
% end
%writetable(struct2table(wavData), '101119A_all_wav_data.xlsx');
%% unsupervised clustering
% k-means
%X_trnP2 = zscore(X_trnP2);
[kmcidx,kmcmeans] = kmeans(X_trnP2,3,'dist','sqeuclidean');

ptsymb = {'bs','g^','md','ro','c+'}; % color scheme
for i = 1:3
    clust = find(kmcidx==i); % identifies clusters
    plot(X_trnP2(clust,1),X_trnP2(clust,2),ptsymb{i}); %plots based on clusters
    hold on
end
%plot(kmcmeans(:,1),kmcmeans(:,2),'ko');
plot(kmcmeans(:,1),kmcmeans(:,2),'kx'); % marks the cluster centers
hold off
title('K-means clustering based on P2_amp and P2_period')
xlabel('P2_ampplitude');
ylabel('P2_period');
grid on
%% hierarchical clustering
% works but need to fix visualization
eucD = pdist(X_trn,'euclidean');
clustTreeEuc = linkage(eucD,'average');
cophenet(clustTreeEuc,eucD);
[h,nodes] = dendrogram(clustTreeEuc,3);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
%
hidx = cluster(clustTreeEuc,'criterion','distance','cutoff',.001);
for i = 1:3
    clust = find(hidx==i);
    plot(X_trn(clust,1),X_trn(clust,2),ptsymb{i});
    hold on
end
hold off
xlabel('P1');
ylabel('P2');
grid on
%% DBSCAN
for eps=0.05:0.005:0.07
    [dbidx,dbcorepts] = dbscan(X_trn,eps,3);
    figure,gscatter(X_trn(:,1),X_trn(:,2),dbidx); 
    title(['eps= ',num2str(eps)])
end
% spectral clustering - spectralcluster() in 2019b

