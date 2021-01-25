clear all; close all; clc

% Start timer
tic

% RVE desired properties
xmin  = -100;
xmax  = 100;
rmin  = 1.0;
rmax  = 1.0;
N     = 5000;
%radii = 1.0;
radii = rmin + (rmax-rmin)*rand(N,1);
min_dist = 2.25*rmax;
tol   = 1e-10;

% Random circle centers
x = xmin + (xmax-xmin).*rand(N,1);
y = xmin + (xmax-xmin).*rand(N,1);
centers      = [x y];

% %subplot(2,2,1)
% figure(1)
% xlim([xmin xmax]);
% ylim([xmin xmax]);
% rectangle('Position',[xmin xmin xmax-xmin xmax-xmin]);
% hold on;
% for i = 1:N
%     %center = centers(i,:);
%     viscircles(centers(i,:),radii(i));
% end
% hold off;

centers_mod  = centers; 
radii_mod    = radii;

% Check if distance between circles < min_dist
gb_count = 0;
lc_count = zeros(N,2);
for i = 1:N            % ith circle
    x_i = centers(i,1);
    y_i = centers(i,2);
    for j=1:N          % loop over all circles
        x_j = centers(j,1);
        y_j = centers(j,2);
        euc_dist = sqrt((x_j - x_i)^2 + (y_j - y_i)^2);
        if euc_dist < min_dist + tol && i ~= j
          gb_count = gb_count + 1;
          lc_count(i,1) = i;
          lc_count(i,2) = lc_count(i,2) + 1;
          overlap(gb_count) = i;
        end
    end
end

lc_counts = sortrows(lc_count,2,'descend');

k = 1;
for i = 1:length(lc_counts(:,1))
    if lc_counts(i,2) > 0
        to_remove(k) = lc_counts(i,1);
        k = k + 1;
    end
end

%centers_mod = centers;
%to_remove   = unique(overlap);
centers_mod(to_remove,:) = [];
radii_mod(to_remove,:) = [];
% Stop timer
toc

%subplot(2,2,2)
figure
xlim([xmin xmax]);
ylim([xmin xmax]);
rectangle('Position',[xmin xmin xmax-xmin xmax-xmin]);
hold on;
for i = 1:length(centers_mod(:,1))
    %center = centers_mod(i,:);
    viscircles(centers_mod(i,:),radii_mod(i));
end
hold off;

%t = table(centers_mod,radii_mod);
%writetable(t,'masterRVE.txt','Delimiter','\t');



    
