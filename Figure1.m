clear all
syms L A
mu1=0.745; mu2=0.104; delta=0.0607; 
alpha=delta./mu2;   %0.4250
bstar=mu2.*(mu1+delta)./delta; % 1.2207
beta=37.8849; 
k=1000./16.6917 %5.991;  

T=3; barT=14; 
p=fix(barT./T); 
q=barT-p.*T;  
c1star=(alpha.*k.*(sqrt(alpha.*beta)-sqrt(delta+mu1)).^2)./((p+1).*mu1);% 135.8841 
c2star=(alpha.*k.*(sqrt(alpha.*beta)-sqrt(delta+mu1)).^2)./(p.*mu1);% 169.8551 

% 算出平衡点
Lstar=k*delta*(beta-bstar)./(mu1*mu2); 
Astar=k*delta.^2*(beta-bstar)./(mu1*mu2*mu2); 
Estar=[Lstar,Astar];    