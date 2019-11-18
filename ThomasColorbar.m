function newmap = ThomasColorbar(minimum, maximum)

coloraxis = linspace(minimum, maximum, 64);
newmap = jet;

% Findind the zero value 
zerocontour = find(abs(coloraxis) == min(abs(coloraxis)));

% Set color of zero only to white  
newmap(zerocontour, :) = [1, 1, 1];
newmap(1:zerocontour, 1) = linspace(0, 1, zerocontour);
newmap(1:zerocontour, 2) = linspace(0, 1, zerocontour);
newmap(1:zerocontour, 3) = linspace(0.5625, 1, zerocontour);
newmap(zerocontour:zerocontour+3, 1) = linspace(1, 0, 4);
newmap(zerocontour:zerocontour+3, 2) = linspace(1, 1, 4);
newmap(zerocontour:zerocontour+3, 3) = linspace(1, 1, 4);

% % Set color of zero and first nearest neighbors to color white
% newmap(zerocontour-1:zerocontour+1, :) = [1, 1, 1; 1, 1, 1; 1, 1, 1];
% newmap(1:zerocontour-1, 1) = linspace(0, 1, zerocontour-1);
% newmap(1:zerocontour-1, 2) = linspace(0, 1, zerocontour-1);
% newmap(1:zerocontour-1, 3) = linspace(0.5625, 1, zerocontour-1);
% newmap(zerocontour+1:26, 1) = linspace(1, 0, 4);
% newmap(zerocontour+1:26, 2) = linspace(1, 1, 4);
% newmap(zerocontour+1:26, 3) = linspace(1, 1, 4);

% newmap(zerocontour, :) = [1, 1, 1];
% newmap(1:zerocontour/2, 1) = linspace(0, 0, zerocontour/2);
% newmap(1:zerocontour/2, 2) = linspace(0, 0, zerocontour/2);
% newmap(1:zerocontour/2, 3) = linspace(143/255, 0.5, zerocontour/2);
% 
% newmap(zerocontour/2+1:zerocontour, 1) = linspace(0, 1, zerocontour/2);
% newmap(zerocontour/2+1:zerocontour, 2) = linspace(0, 1, zerocontour/2);
% newmap(zerocontour/2+1:zerocontour, 3) = linspace(0.5, 1, zerocontour/2);

% newmap(1:zerocontour, 1) = linspace(0, 0, zerocontour);
% newmap(1:zerocontour, 2) = linspace(0, 1, zerocontour);
% newmap(1:zerocontour, 3) = linspace(143/255, 1, zerocontour);

% newmap(1:zerocontour, 1) = linspace(0, 1, zerocontour);
% newmap(1:zerocontour, 2) = linspace(0, 1, zerocontour);
% newmap(1:zerocontour, 3) = linspace(0.5625, 1, zerocontour);
% newmap(zerocontour:zerocontour+3, 1) = linspace(1, 0, 4);
% newmap(zerocontour:zerocontour+3, 2) = linspace(1, 1, 4);
% newmap(zerocontour:zerocontour+3, 3) = linspace(1, 1, 4);