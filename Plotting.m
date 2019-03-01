clear all
clc

% This code simulates a patrol of a grid where expected commute time
% averaged over all pairs of nodes is minimized.


% Optimal edge probabilities
%M = csvread("weightsSmall.csv");
M = xlsread("weightsLarge.xlsx");


% INPUT WIDTH OF SQUARE GRID
%k = 8;
k = 30;

% TIME STEPS TO SIMULATE
time = 1000;

% LENGTH OF TAIL FOR PLOTTING
tail = time;

% WHERE IN THE GRID TO BEGIN
startNode = 465;




% SHOULDN'T NEED TO CHANGE ANYTHING BELOW THIS POINT

n = k^2;

% this is the weighted adjacency matrix
weightedA = zeros(n,n);
for i = 1:length(M)
   if M(i,1) ~= 0
        weightedA(M(i, 1),M(i, 2)) = M(i, 3);    
   end
end


% This alters the weighted adjacency matrix so it is easier to randomly
% sample an edge later. Sometimes the rows don't add exactly to 1 because
% of the rounding. If you get a running error, this is why. Rerun and it
% should be fine.
weightedB = zeros(n,n);
for i = 1:n
   previous = 0;
   for j = 1:n
       if weightedA(i, j) ~= 0
           weightedB(i,j) = previous + weightedA(i, j);
           previous = weightedB(i, j);
       end
   end
end



% path in node form
Path = zeros(time,1);
Path(1) = startNode;
for t = 2:time
    x = rand;
    for j = 1:n
       if x < weightedB(Path(t-1), j)
            Path(t) = j;
            break
       end
    end
end


% path in x and y coordinates
travel = zeros(time,2);
for i = 1:length(Path)
    if mod(Path(i),k) ~= 0
        travel(i,1) = k*(Path(i)/k - floor(Path(i)/k)); 
        travel(i,2) = ceil(Path(i)/k); 
    else
        travel(i,1) = k; 
        travel(i,2) = ceil(Path(i)/k);

    end

end


% plotting
figure
for t = 1:time
    
    if t > tail
        plot(travel(t-tail:t,1), travel(t-tail:t,2), '-');
        hold on
    else
        plot(travel(1:t,1), travel(1:t,2), '-'); 
        hold on
    end
    xlim([0 k+1])
    ylim([0 k+1])
    plot(linspace(1,k), linspace(1,1),'k-');
    hold on
    plot(linspace(1,k), linspace(k,k),'k-');
    hold on
    plot(linspace(1,1), linspace(1,k),'k-');
    hold on
    plot(linspace(k,k), linspace(1,k),'k-');
    hold on
    xlim([0 k+1])
    ylim([0 k+1])
    grid on
    plot(travel(t,1), travel(t,2), '*');
    hold off
    
    pause(0.02);
end

