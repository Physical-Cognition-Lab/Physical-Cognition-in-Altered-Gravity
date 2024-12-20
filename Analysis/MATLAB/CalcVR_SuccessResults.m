clear;clc; close all;
load('FlexGravityVR.mat');
subjects=unique(vr_data(:,1));
sub_success_ratio = zeros(length(subjects),6);
sub_success_ratio_by_distances = zeros(length(subjects),5,2);
percent_of_not_thrown_sub_data_per_subject = [];
success_per_planet_per_distance_per_subject = [];

for s_ix=1:length(subjects)
    sub_num=subjects(s_ix);
    
    % take only the ones that were thrown
    sub_data = vr_data(vr_data(:,1)==sub_num & ...
                       vr_data(:,8)==1 & ...
                       vr_data(:,5)~=-1,[2 5 4 3]); %planet | distance | success 


   percent_of_not_thrown_sub_data_per_subject = [percent_of_not_thrown_sub_data_per_subject; size(vr_data(vr_data(:,1)==sub_num & ...
                       vr_data(:,8)==0 & ...
                       vr_data(:,5)~=-1,[3 2 5 4]),1)/size(sub_data,1)*100]; %how many were not thrown

   planets = [0 1 2 3 4];
   planet_dist_rate = [];
    for p_ix=1:5
       planet = planets(p_ix);
       planet_sub_data = sub_data(sub_data(:,1)==planet,2:3);
       for dist=1:4
            dist_learning_sub_data = planet_sub_data(planet_sub_data(:,1)==dist,2);
            planet_dist_rate = [planet_dist_rate sum(dist_learning_sub_data==1)/size(dist_learning_sub_data,1)*100];
       end 
    end
    success_per_planet_per_distance_per_subject = [success_per_planet_per_distance_per_subject; planet_dist_rate];
end

adapt_data = [];
for s_ix=1:size(success_per_planet_per_distance_per_subject,1)
    subj_adapt_index=[];
    for p_ix=[9 17 1 13] % first is neptune, then earth, then moon, then jupiter, then venus
       adapt_planet_dist1 = ((success_per_planet_per_distance_per_subject(s_ix,p_ix) - success_per_planet_per_distance_per_subject(s_ix,5)) / ...
                            (success_per_planet_per_distance_per_subject(s_ix,p_ix)+success_per_planet_per_distance_per_subject(s_ix,5)));
       adapt_planet_dist2 = ((success_per_planet_per_distance_per_subject(s_ix,p_ix+1) - success_per_planet_per_distance_per_subject(s_ix,6)) / ...
                            (success_per_planet_per_distance_per_subject(s_ix,p_ix+1)+success_per_planet_per_distance_per_subject(s_ix,6)));
       adapt_planet_dist3 = ((success_per_planet_per_distance_per_subject(s_ix,p_ix+2) - success_per_planet_per_distance_per_subject(s_ix,7)) / ...
                            (success_per_planet_per_distance_per_subject(s_ix,p_ix+2)+success_per_planet_per_distance_per_subject(s_ix,7)));
       adapt_planet_dist4 = ((success_per_planet_per_distance_per_subject(s_ix,p_ix+3) - success_per_planet_per_distance_per_subject(s_ix,8)) / ...
                            (success_per_planet_per_distance_per_subject(s_ix,p_ix+3)+success_per_planet_per_distance_per_subject(s_ix,8)));
      if (isnan(adapt_planet_dist1)) adapt_planet_dist1 = 0; end;
      if (isnan(adapt_planet_dist2)) adapt_planet_dist2 = 0; end;
      if (isnan(adapt_planet_dist3)) adapt_planet_dist3 = 0; end;
      if (isnan(adapt_planet_dist4)) adapt_planet_dist4 = 0; end;
      subj_adapt_index = [subj_adapt_index; (adapt_planet_dist1 + adapt_planet_dist2 + adapt_planet_dist3 + adapt_planet_dist4)/4];
    end
    adapt_data = [adapt_data; subj_adapt_index'];
end



% Step 1: Find optimal number of clusters
maxClusters = 25; % You can adjust this based on your expectation of the dataset
optimalRatio = 0;
optimalNumClusters = 0;
for k=2:maxClusters
    [idx,~,sumd] = kmeans(adapt_data, k, 'Display', 'final', 'Replicates', 10);
    withinClusterDist = sum(sumd);
    betweenClusterDist = pdist2(mean(adapt_data,1), adapt_data(idx==mode(idx),:), 'euclidean');
    ratio = withinClusterDist / mean(betweenClusterDist);
    if ratio > optimalRatio
        optimalRatio = ratio;
        optimalNumClusters = k;
    end
end

fprintf('Optimal number of clusters: %d\n', optimalNumClusters);

% Step 2: Perform clustering with the optimal number of clusters
[idx,C] = kmeans(adapt_data, optimalNumClusters, 'Display', 'final', 'Replicates', 10);

% Step 3: Plot the similarity matrix ordered according to clusters
D = pdist2(adapt_data, adapt_data, 'euclidean');
[~, sortIdx] = sort(idx);
sortedD = D(sortIdx, sortIdx);
figure;
colormap('autumn');
imagesc(sortedD);
title('Similarity Matrix Ordered by Clustering');
colorbar;

% Step 4: Print the average patterns for each cluster

% Calculating global y-axis limits
allMeans = [];
allStdDevs = [];
for i=1:optimalNumClusters
    clusterData = adapt_data(idx==i,:);
    meanPattern = mean(clusterData,1);
    stdPattern = std(clusterData, 0, 1);
    allMeans = [allMeans; meanPattern];
    allStdDevs = [allStdDevs; stdPattern];
end
globalYMin = min(allMeans - allStdDevs, [], 'all');
globalYMax = max(allMeans + allStdDevs, [], 'all');

figure;
for i=1:optimalNumClusters
    subplot(1,optimalNumClusters,i);
    clusterData = adapt_data(idx==i,:);
    meanPattern = mean(clusterData,1);
    stdPattern = std(clusterData, 0, 1)/sqrt(25);
    bar(meanPattern); % Plotting the mean pattern as a bar graph
    hold on;
    errorbar(1:length(meanPattern), meanPattern, stdPattern, 'k', 'linestyle', 'none'); % Adding error bars
    hold off;
    title(sprintf('Cluster %d', i));
    ylim([globalYMin globalYMax]); % Ensuring consistent y-axis
end

% Additions to the beginning of the script to initialize storage for the ratios
ratios = zeros(maxClusters-1, 1); % Initialize an array to store the ratio for each number of clusters

% Modify the loop where the optimal number of clusters is determined to store the ratio
for k=2:maxClusters
    [idx,~,sumd] = kmeans(adapt_data, k, 'Display', 'final', 'Replicates', 10);
    withinClusterDist = sum(sumd);
    betweenClusterDist = pdist2(mean(adapt_data,1), adapt_data(idx==mode(idx),:), 'euclidean');
    ratio = withinClusterDist / mean(betweenClusterDist);
    ratios(k-1) = ratio; % Store the ratio
    if ratio > optimalRatio
        optimalRatio = ratio;
        optimalNumClusters = k;
    end
end

% Plot the change in within-cluster to between-cluster ratio
figure;
plot(2:maxClusters, ratios, '-o');
xlabel('Number of Clusters');
ylabel('Within-Cluster to Between-Cluster Ratio');
title('Change in Ratio as a Function of Number of Clusters');
grid on;

% Ensure the rest of the script remains unchanged


%save('SuccessRatios.mat','sub_success_ratio_by_distances','percent_of_not_thrown_sub_data_per_subject','sub_success_ratio','sub_learning_ratio_by_distances','sub_learning_ratio');
