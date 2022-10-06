function  res = SENSE3D(sen)
%res = ESPIRiT(EigenVecs, [EigenVals,  imSize)
%
%	Implementation of the ESPIRiT maps projection
%

%
% (c) Michael Lustig 2013
     

res.sen = sen;
res.adjoint = 0;
res = class(res,'SENSE3D');

