function [costs, lifetimes] = capacity2costs(capacity, slope, intercept)

costs = [0];
lifetimes = [];

%% capped at 1
for j = 1:5
capacity1= max(capacity,0);
capacity1(capacity1>1)=1;

capacity1=[0,capacity1,0]; %to fix edges
%count builds and decomm
build  = max(0, capacity1(2:end)-capacity1(1:end-1)); 
decomm =  max(0, capacity1(1:end-1)-capacity1(2:end));

%remove zero values and keep only construction time and decomm time
idx_build = find(build>0);
idx_decomm = find(decomm>0);

new_lifes = idx_decomm-idx_build;
lifetimes = [lifetimes, new_lifes];

for i=1:length(new_lifes)
    costs = [costs, intercept + slope*new_lifes(i)/12];
end

capacity = capacity-1;
end
