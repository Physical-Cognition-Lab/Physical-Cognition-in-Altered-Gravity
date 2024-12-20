clear; clc; close all;
% Step 1: Read the CSV file
filename = 'VirtualToolsData_VRExperiment.csv';
data = readtable(filename);

% Filter out rows where block type is 'practice'
data = data(~strcmp(data{:, 5}, 'practice'), :);

% Extracting the relevant columns for clarity
participantIDs = data{:, 1};
blockTypes = data{:, 5};
numAttempts = data{:, 7};
games = data{:, 8};
gravityLevels = data{:, 9};
successes = data{:, 14};

% Unique identifiers
participants = unique(participantIDs);
gravities = unique(gravityLevels);
gamesList = unique(games);

% Initialize the results table for gravity analysis
resultsTableGravity_success = array2table(zeros(length(participants), length(gravities)), 'VariableNames', string(gravities), 'RowNames', string(participants));
resultsTableGravity_attempts = array2table(zeros(length(participants), length(gravities)), 'VariableNames', string(gravities), 'RowNames', string(participants));

resultsGamesGravity_success = array2table(zeros(length(gamesList), length(gravities)-1), 'VariableNames', string(gravities([2 3],:)), 'RowNames', string(gamesList));
resultsGamesGravity_attempts = array2table(zeros(length(gamesList), length(gravities)-1), 'VariableNames', string(gravities([2 3],:)), 'RowNames', string(gamesList));

for i = 1:length(participants)
    for g = 1:length(gravities)
        gravity_data = data(participantIDs == participants(i) & gravityLevels == gravities(g),:);
        trials = unique(gravity_data.TRIAL);
        average_trials = [];
        for t=1:length(trials)
            trial_data = gravity_data(gravity_data.TRIAL==trials(t),:); 
            num_of_attempts = trial_data{end,"ATTEMPT"};
            success = trial_data{end,"SUCCESS"};
            average_trials = [average_trials; num_of_attempts success];
        end

        average_trials = mean(average_trials);
        if (isnan(average_trials))
            resultsTableGravity_success(i,g) = {-1};
            resultsTableGravity_attempts(i,g) = {-1};
        else
            resultsTableGravity_success(i,g) = {average_trials(2)};
            resultsTableGravity_attempts(i,g) = {average_trials(1)};
        end
    end
end

low_gravity_adaptability_success = (resultsTableGravity_success{:,1}-resultsTableGravity_success{:,2})./(resultsTableGravity_success{:,1}+resultsTableGravity_success{:,2});
high_gravity_adaptability_success = (resultsTableGravity_success{:,3}-resultsTableGravity_success{:,2})./(resultsTableGravity_success{:,3}+resultsTableGravity_success{:,2});
adapt_data_success = [low_gravity_adaptability_success high_gravity_adaptability_success];

low_gravity_adaptability_attempts = (resultsTableGravity_attempts{:,1}-resultsTableGravity_attempts{:,2})./(resultsTableGravity_attempts{:,1}+resultsTableGravity_attempts{:,2});
high_gravity_adaptability_attempts = (resultsTableGravity_attempts{:,3}-resultsTableGravity_attempts{:,2})./(resultsTableGravity_attempts{:,3}+resultsTableGravity_attempts{:,2});
adapt_data_attempts = [low_gravity_adaptability_attempts high_gravity_adaptability_attempts];



% Step 3: Print a bar graph with error bars comparing the three gravities

avgAttemptsPerGravity = mean(table2array(resultsTableGravity_attempts), 1);
stdDevAttemptsPerGravity = std(table2array(resultsTableGravity_attempts), 0, 1) / sqrt(25); % STD divided by the root of 25

figure;
bar(gravities, avgAttemptsPerGravity);
hold on;
errorbar(gravities, avgAttemptsPerGravity, stdDevAttemptsPerGravity, 'k', 'linestyle', 'none');
hold off;
xlabel('Gravity');
ylabel('Average Attempts');
title('Average Attempts per Gravity with Error Bars');

avgAttemptsPerGravity = mean(table2array(resultsTableGravity_success), 1);
stdDevAttemptsPerGravity = std(table2array(resultsTableGravity_success), 0, 1) / sqrt(25); % STD divided by the root of 25

figure;
bar(gravities, avgAttemptsPerGravity);
hold on;
errorbar(gravities, avgAttemptsPerGravity, stdDevAttemptsPerGravity, 'k', 'linestyle', 'none');
hold off;
xlabel('Gravity');
ylabel('Average Success');
title('Average Success per Gravity with Error Bars');





% Initialize the results table for game analysis
resultsTableGame = array2table(zeros(length(gamesList), length(gravities)), 'VariableNames', string(gravities), 'RowNames', string(gamesList));

adapt_data = adapt_data_success;
% Step 1: Find optimal number of clusters
maxClusters = 10; % You can adjust this based on your expectation of the dataset
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
sum(idx==1)
sum(idx==2)
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















% Initialize the results table for game analysis
resultsTableGame = array2table(zeros(length(gamesList), length(gravities)), 'VariableNames', string(gravities), 'RowNames', string(gamesList));

adapt_data = adapt_data_attempts;
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
sum(idx==1)
sum(idx==2)
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
ylim([0 8]);
grid on;

