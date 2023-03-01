clear all
syms L A
clear all
syms L A
mu1=0.745; mu2=0.104; delta=0.0607; 
alpha=delta./mu2;    
bstar=mu2.*(mu1+delta)./delta;  
beta=37.8849; 
k=1000./16.6917;   
T=3; barT=14; 
p=fix(barT./T); 
q=barT-p.*T;  
c1star=(alpha.*k.*(sqrt(alpha.*beta)-sqrt(delta+mu1)).^2)./((p+1).*mu1);% 135.8841 
c2star=(alpha.*k.*(sqrt(alpha.*beta)-sqrt(delta+mu1)).^2)./(p.*mu1);% 169.8551 
c=146.61; 

n=5000; 

% The unstable periodic solution 
iniL(1)=285.5739; iniA(1)=156.3121;
F=@(t,y) [beta*(y(2)^2/(y(2)+(p+1)*c))-mu1*(1+y(1)/k)*y(1)-delta*y(1); delta*y(1)-mu2*y(2)];                         % The equation for q> t > 0.
G=@(t,y) [beta.*(y(2).^2./(y(2)+p.*c))-mu1.*(1+y(1)./k).*y(1)-delta.*y(1); delta.*y(1)-mu2.*y(2)];                           % The equation for T> t > q.
for i=1:n 
    if mod(i,2)==1   
       [t,y]=ode45(F,[T*fix(i/2),T*fix(i/2)+q],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2); 
       subplot(121); 
       plot3(t,y(:,1),y(:,2),'r','Linewidth',2); hold on 
    elseif mod(i,2)==0 
       [t,y]=ode45(G,[T*(i/2)-(T-q), T*(i/2)],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2);
       subplot(121); 
       plot3(t,y(:,1),y(:,2),'b','Linewidth',2); hold on 
    end 
end
 
% Upper unstable PS 
iniL(1)=293; iniA(1)=160;
F=@(t,y) [beta*(y(2)^2/(y(2)+(p+1)*c))-mu1*(1+y(1)/k)*y(1)-delta*y(1); delta*y(1)-mu2*y(2)];                         % The equation for q> t > 0.
G=@(t,y) [beta.*(y(2).^2./(y(2)+p.*c))-mu1.*(1+y(1)./k).*y(1)-delta.*y(1); delta.*y(1)-mu2.*y(2)];                           % The equation for T> t > q.
for i=1:n 
    if mod(i,2)==1   
       [t,y]=ode45(F,[T*fix(i/2),T*fix(i/2)+q],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2); 
       subplot(121); 
       plot3(t,y(:,1),y(:,2),'Linewidth',1); hold on 
    elseif mod(i,2)==0 
       [t,y]=ode45(G,[T*(i/2)-(T-q), T*(i/2)],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2);
       subplot(121); 
       plot3(t,y(:,1),y(:,2),'Linewidth',1); hold on 
    end 
end

% Lower unstable PS 
iniL(1)=285; iniA(1)=156;
F=@(t,y) [beta*(y(2)^2/(y(2)+(p+1)*c))-mu1*(1+y(1)/k)*y(1)-delta*y(1); delta*y(1)-mu2*y(2)];                         % The equation for q> t > 0.
G=@(t,y) [beta.*(y(2).^2./(y(2)+p.*c))-mu1.*(1+y(1)./k).*y(1)-delta.*y(1); delta.*y(1)-mu2.*y(2)];                           % The equation for T> t > q.
for i=1:n 
    if mod(i,2)==1   
       [t,y]=ode45(F,[T*fix(i/2),T*fix(i/2)+q],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2); 
       subplot(121); 
       plot3(t,y(:,1),y(:,2),'Linewidth',1); hold on 
    elseif mod(i,2)==0 
       [t,y]=ode45(G,[T*(i/2)-(T-q), T*(i/2)],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2);
       subplot(121); 
       plot3(t,y(:,1),y(:,2),'Linewidth',1); hold on 
    end 
end
 
% Stable periodic solution 
n=12000; 
iniL(1)=293.6136; iniA(1)=160.7803;
F=@(t,y) [beta*(y(2)^2/(y(2)+(p+1)*c))-mu1*(1+y(1)/k)*y(1)-delta*y(1); delta*y(1)-mu2*y(2)];                         % The equation for q> t > 0.
G=@(t,y) [beta.*(y(2).^2./(y(2)+p.*c))-mu1.*(1+y(1)./k).*y(1)-delta.*y(1); delta.*y(1)-mu2.*y(2)];                           % The equation for T> t > q.
for i=1:n 
    if mod(i,2)==1   
       [t,y]=ode45(F,[T*fix(i/2),T*fix(i/2)+q],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2); 
       subplot(122); 
       plot3(t,y(:,1),y(:,2),'r','Linewidth',2); hold on 
    elseif mod(i,2)==0 
       [t,y]=ode45(G,[T*(i/2)-(T-q), T*(i/2)],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2);
       subplot(122); 
       plot3(t,y(:,1),y(:,2),'b','Linewidth',2); hold on 
    end 
end
% 
% Upper stable Ps 
iniL(1)=295; iniA(1)=165;
F=@(t,y) [beta*(y(2)^2/(y(2)+(p+1)*c))-mu1*(1+y(1)/k)*y(1)-delta*y(1); delta*y(1)-mu2*y(2)];                         % The equation for q> t > 0.
G=@(t,y) [beta.*(y(2).^2./(y(2)+p.*c))-mu1.*(1+y(1)./k).*y(1)-delta.*y(1); delta.*y(1)-mu2.*y(2)];                           % The equation for T> t > q.
for i=1:n 
    if mod(i,2)==1   
       [t,y]=ode45(F,[T*fix(i/2),T*fix(i/2)+q],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2); 
       subplot(122); 
       plot3(t,y(:,1),y(:,2),'Linewidth',1); hold on 
    elseif mod(i,2)==0 
       [t,y]=ode45(G,[T*(i/2)-(T-q), T*(i/2)],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2);
       subplot(122); 
       plot3(t,y(:,1),y(:,2),'Linewidth',1); hold on 
    end 
end
 
% Lower stable Ps 
iniL(1)=293; iniA(1)=160;
F=@(t,y) [beta*(y(2)^2/(y(2)+(p+1)*c))-mu1*(1+y(1)/k)*y(1)-delta*y(1); delta*y(1)-mu2*y(2)];                         % The equation for q> t > 0.
G=@(t,y) [beta.*(y(2).^2./(y(2)+p.*c))-mu1.*(1+y(1)./k).*y(1)-delta.*y(1); delta.*y(1)-mu2.*y(2)];                           % The equation for T> t > q.
for i=1:n 
    if mod(i,2)==1   
       [t,y]=ode45(F,[T*fix(i/2),T*fix(i/2)+q],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2); 
       subplot(122); 
       plot3(t,y(:,1),y(:,2),'Linewidth',1); hold on 
    elseif mod(i,2)==0 
       [t,y]=ode45(G,[T*(i/2)-(T-q), T*(i/2)],[iniL(i) iniA(i)]); 
       iniL(i+1)=y(end,1);
       iniA(i+1)=y(end,2);
       subplot(122); 
       plot3(t,y(:,1),y(:,2),'Linewidth',1); hold on 
    end 
end
 
subplot(121); 
xlabel('','Interpreter','latex','String','$t$'); 
ylabel('','Interpreter','latex','String','$L(t)$'); 
zlabel('','Interpreter','latex','String','$A(t)$');
title('(a)'); 

subplot(122); 
xlabel('','Interpreter','latex','String','$t$'); 
ylabel('','Interpreter','latex','String','$L(t)$'); 
zlabel('','Interpreter','latex','String','$A(t)$');
title('(b)'); 