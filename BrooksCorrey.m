
 function [swof]=BrooksCorrey(Pd,landau,SwIrr,SoIrr,KrwEnd,KroEnd)
% Input: SI Unit ; Output: Field Unit


%--------------------

% Pd=100;    
% landau=1;    
% SwIrr=0.2;    
% SoIrr=0.2;    
%------------------
pointNumber=10;

Sw=SwIrr:((1-SoIrr-SwIrr)/pointNumber):(1-SoIrr);

SwEffective=((Sw-SwIrr)/(1-SwIrr-SoIrr));

% Pd: Entry Pressure
PcBrCoStatic=Pd*(SwEffective.^(-1/landau));

KrWater=KrwEnd*SwEffective.^(3+2/landau);

KrOil=KroEnd*((1-SwEffective).^2).*(1-(SwEffective.^(1+2/landau)));


% Correction for Pc at SwIrr

slope1=PcBrCoStatic(3)-PcBrCoStatic(4);
slope2=PcBrCoStatic(2)-PcBrCoStatic(3);
slopep=slope2/slope1;
slope3=slope2*slopep;

PcBrCoStatic(1)=PcBrCoStatic(2)+slope3;

% Pascals to Psi:
% PcBrCoStatic=PcBrCoStatic*(14.7/1.01e5);


swof=[Sw' , KrWater' , KrOil' , PcBrCoStatic' ];
swof(1,2)=0;
swof(end,3)=0;

% 
 end