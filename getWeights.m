clear all
clc

% Boyd & Vandenberghe "Convex Optimization"
% Joëlle Skaf - 09/26/05
%
% The 'fastest mixing Markov chain problem' is to find a transition
% probability matrix P on a graph E that minimizes the mixing rate r, where
% r = max{ lambda_2, -lambda_n } with lambda_1>=...>=lambda_n being the
% eigenvalues of P.

% you need to have the cvx package to run this

% runs in linear time with respect to adding nodes. 

%INPUT WIDTH OF SQUARE GRID
%ONLY THING YOU NEED TO CHANGE
k = 10;

%number of cells
n = k^2;

%adjacency matrix
E = zeros(n,n);

for cell = 1:n
    
    % has all eight neighbors
    if(mod(cell,k) ~= 0 && mod(cell-1,k) ~= 0 && cell < (n-k) && cell > (k+1))
        E(cell, cell - 1) = 1;
        E(cell, cell + 1) = 1;
        E(cell, cell - k) = 1;
        E(cell, cell + k) = 1;

        E(cell, cell - (k+1)) = 1;
        E(cell, cell + (k+1)) = 1;
        E(cell, cell - (k-1)) = 1;
        E(cell, cell + (k-1)) = 1; 
              
    % bottom left corner
    elseif(cell == 1)
        E(cell, cell + 1) = 1;
        E(cell, cell + k) = 1;
        E(cell, cell + (k+1)) = 1;

    
        
        
    % bottom right corner
    elseif(cell == k)
        E(cell, cell - 1) = 1;
        E(cell, cell + k) = 1;
        E(cell, cell + (k-1)) = 1;
        
        
    % top left corner
    elseif(cell == (n-k+1))
        E(cell, cell + 1) = 1;
        E(cell, cell - k) = 1;
        E(cell, cell - (k-1)) = 1;
        
        
    % top right corner
    elseif(cell == n)
        E(cell, cell - 1) = 1;
        E(cell, cell - k) = 1;
        E(cell, cell - (k+1)) = 1;
        
        
    % is on left edge
    elseif(mod(cell-1,k) == 0)
        E(cell, cell + 1) = 1;
        E(cell, cell - k) = 1;
        E(cell, cell + k) = 1;
        E(cell, cell + (k+1)) = 1;
        E(cell, cell - (k-1)) = 1;
    
        
    % is on right edge
    elseif(mod(cell,k) == 0)
        E(cell, cell - 1) = 1;
        E(cell, cell - k) = 1;
        E(cell, cell + k) = 1;

        E(cell, cell - (k+1)) = 1;
        E(cell, cell + (k-1)) = 1;
        
        
    % is on top edge
    elseif(cell > (n-k+1))
        E(cell, cell - 1) = 1;
        E(cell, cell + 1) = 1;
        E(cell, cell - k) = 1;

        E(cell, cell - (k+1)) = 1;
        E(cell, cell - (k-1)) = 1;

    
    % is on bottom edge
    elseif(cell < k)
        E(cell, cell - 1) = 1;
        E(cell, cell + 1) = 1;
        E(cell, cell + k) = 1;
        E(cell, cell + (k+1)) = 1;
        E(cell, cell + (k-1)) = 1;
    end
    
end


% Create and solve model
cvx_begin
    variable P(n,n) symmetric
    minimize(norm(P - (1/n)*ones(n)))
    P*ones(n,1) == ones(n,1);
    P >= 0;
    P(E==0) == 0;
cvx_end
e = flipud(eig(P));
r = max(e(2), -e(n));


% Display results
disp('------------------------------------------------------------------------');
disp('The transition probability matrix of the optimal Markov chain is: ');

disp('The optimal mixing rate is: ');
disp(r);

% rounds probabilities to 4 decimal places
Y = round(P, 4);
disp(Y);







