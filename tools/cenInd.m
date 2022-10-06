function ind = cenInd(Nd,Cd)
% returns Cd central indices from total Nd indices
if(nargin<2)
    Cd = 1;
end
ind = floor(Nd/2)+1+ceil(-Cd/2):floor(Nd/2)+ceil(Cd/2);
