%Se cargan las ecuaciones diferenciales que describen cada uno de los 
% estados y también la serie de tiempo que contiene los casos nuevos por 
% esta enfermedad (periodo de tiempo en semanas).

syms  Ms(t) Me(t) Mi(t) Hs(t) He(t) Hi(t) Hr(t) Hit(t)

syms lambda beta_m mu_m theta_m mu_h beta_h theta_h...
    gamma_h

H=Hs+He+Hi+Hr;
M=Ms+Me+Mi;
ode1 = diff(Ms) == lambda - beta_m*Hi*Ms/H - (mu_m)*Ms;
ode2 = diff(Me) == beta_m*Hi*Ms/H - (theta_m+mu_m)*Me;
ode3 = diff(Mi) == theta_m*Me - mu_m*Mi;
ode4 = diff(Hs) == -beta_h*Mi*Hs/M + (He+Hi+Hr)*mu_h;
ode5 = diff(He) == beta_h*Mi*Hs/M - (theta_h+mu_h)*He;
ode6 = diff(Hi) == theta_h*He - (gamma_h+mu_h)*Hi;
ode7 = diff(Hr) == gamma_h*Hi - mu_h*Hr;
ode8 = diff(Hit) == theta_h*He;
odes=[ode1; ode2 ;ode3; ode4; ode5 ;ode6; ode7; ode8];
vars=[Hit Hi Me Hr Hs He Ms Mi];
opts = odeset('NonNegative',1:8); 
load Range7.mat


[T,~]=gsua_dpmat(odes,vars,[0 80],'7m','output',1,'opt',opts,'Range',Range);
load('DataBello_full.mat');
ydata=DataBello.cases(520:600)';
xdata=linspace(0,length(ydata)-1,length(ydata));
ydata2=ydata;
for i=1:length(ydata)
ydata2(i)=sum(ydata(1:i));
end
%% Análisis de Sensibilidad

%Queremos observar la variable Hit (casos acumulados)
T.Properties.CustomProperties.output=1;
ynom= gsua_eval(T.Nominal,T);
M= gsua_dmatrix(T,500);
Tsa= gsua_sa(M,T,'parallel',false, 'SensMethod','Saltelli','ynom',ynom); 
sum(Tsa.Si)/sum(abs(Tsa.Si))
clf 
gsua_plot('Bar',Tsa,Tsa.STi)
savefig('indicesconfi.fig')

%% 
% Las primeras 8 filas de la tabla Range corresponde a cada una de las 8
% vars respectivamente 
% La fila 9: beta_h, fila 10: beta_m, fila 11: gamma_h, fila 2: lambda fila
% 14: mu_m, fila 15: theta_h, fila 16: theta_m

%Según el análisis de sensibilidad se quiere fijar Mi0, Ms0, mu_m, Me0, He0

%Orden de estas variables en la tabla T:
% fila 6: Mi0, fila 5: Ms0, fila 11: mu_m, fila 1: Me0, fila 4: He0
%Orden de estas variables en Range
% fila 8: Mi0, fila 7: Ms0, fila 14: mu_m, fila 3: Me0, fila 6: He0


%Hay que cambiar el orden de las variables y por ende, también su posición
%en la matriz de Ranges

%load Range7
%vars=[Hit Hi Me Hr Hs He Ms Mi];
vars= [Hit Hi Me He Hr Hs Ms Mi];



Range2= Range;
Range(4:6,:)= [Range2(6,:); Range2(4:5,:)];

%Una vez ordenados se puede proceder a fijar las variables y parámetros que
%se definieron anteriormente
Range(3,:)= T.Nominal('Me0');
Range(4,:)= T.Nominal('He0');
Range(15,:)= T.Nominal('theta_h');
Range(end,:)= T.Nominal('theta_m');

[T,~]=gsua_dpmat(odes,vars,[0 80],'7m','output',1,'opt',opts,'Range',Range);

%% 

opt=optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7,res]=gsua_pe(T,xdata,ydata2,'N',100,'opt',opt);
save('Results7.mat','T7','res','xdata','ydata2')
%% 
%load Results7
y7=gsua_eval(T7.Estlsqc,T7,xdata,ydata2); %aca
savefig('familiacurvas.fig')

%% 

%Queremos ver los nuevos casos de contagios
ydata_a= [ydata2(1) diff(ydata2)];
y=y7;
for i=2:size(y7,2)
y(:,i)=y7(:,i)-y7(:,i-1);
end
clf
plot(xdata, ydata_a)
hold on
plot(xdata,y(1,:,1))
legend('Casos nuevos reales', 'estimados')

%% 
T7.Properties.CustomProperties.output=1;
T7.Nominal=T7.Estlsqc(:,1);
T7_2=gsua_ia(T7,T7.Estlsqc,false,true);
%% check
M=gsua_dmatrix(T7,100)
Ua = gsua_ua(M, T7, 'parallel', false, 'ynom',ydata2);

%% 
ynom= gsua_eval(T7.Nominal,T7);
M= gsua_dmatrix(T7,500);
Tsa= gsua_sa(M,T7,'parallel',false, 'SensMethod','Saltelli','ynom',ynom); 
sum(Tsa.Si)/sum(abs(Tsa.Si))
clf 
gsua_plot('Bar',Tsa,Tsa.STi)
%% 
gsua_eval(M',T7)


%% 
Range3= Range;
%% 

beta_m_estimado= T7.Nominal('beta_m');
%Se fijan por las estimaciones
Range(10,:)= beta_m_estimado;
 
%vars= [Hit Hi Me He Hr Hs Ms Mi];

Range(5,1)= Range(5,2)- (Range(5,2)-Range(5,1))*1.2;
Range(11,2)= Range(11,1) + (Range(11,2)-Range(11,1))*1.2;
Range(8,1)= 0;
Range(14,2)= Range(14,1) +(Range(14,2)-Range(14,1))*1.2;
Range(12,1)= 0;
Range(7,2) = Range(7,1) + (Range(7,2)-Range(7,1))*1.2;


%% Segunda Ronda de Estimacion 
Range4=Range;
%% 


[T,~]=gsua_dpmat(odes,vars,[0 80],'7m','output',1,'opt',opts,'Range',Range);

%% 
M=gsua_dmatrix(T,500);
[T7,J7,Y7]=gsua_sa(M,T,'parallel',false,'SensMethod','Saltelli');
save('Results7_sa2.mat','T7','J7','Y7','xdata')
%% 

opt=optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7,res]=gsua_pe(T,xdata,ydata2,'N',100,'opt',opt);
save('Results7_fix1.mat','T7','res','xdata','ydata2')
%% 
%load Results7_fix
y7=gsua_eval(T7.Estlsqc,T7,xdata,ydata2); %aca
savefig('familiacurvas.fig')

%% 

%Queremos ver los nuevos casos de contagios
ydata_a= [ydata2(1) diff(ydata2)];
y=y7;
for i=2:size(y7,2)
y(:,i)=y7(:,i)-y7(:,i-1);
end
clf
plot(xdata, ydata_a)
hold on
plot(xdata,y(1,:,1))
legend('Casos nuevos reales', 'estimados')
%% 
T7.Properties.CustomProperties.output=1;
T7.Nominal=T7.Estlsqc(:,1);
T7_2=gsua_ia(T7,T7.Estlsqc,false,true);
%% check
M=gsua_dmatrix(T7,100)
Ua = gsua_ua(M, T7, 'parallel', false, 'ynom',ydata2);

%% 
ynom= gsua_eval(T7.Nominal,T7);
M= gsua_dmatrix(T7,500);
Tsa= gsua_sa(M,T7,'parallel',false, 'SensMethod','Saltelli','ynom',ynom); 
sum(Tsa.Si)/sum(abs(Tsa.Si))
clf 
gsua_plot('Bar',Tsa,Tsa.STi)
%% 
gsua_eval(M',T7)

%% Tercera Ronda de estimación
Range5= Range;
%% 
gamma_h_estimado= T7.Nominal('gamma_h');
%Se fijan por las estimaciones
Range(11,:)= gamma_h_estimado;
%% 
mu_m_estimado= T7.Nominal('mu_m');
Range(14,:)= mu_m_estimado;
%% 

%vars=[Hit Hi Me He Hr Hs Ms Mi];;
vars= [Hit Hi Me He Ms Hr Hs Mi];

Range(5:7,:)= [Range5(7,:); Range5(5:6,:)] ;
Range(6,1)= Range(6,2)- (Range(6,2)-Range(6,1))*1.2;

%% 
[T,~]=gsua_dpmat(odes,vars,[0 80],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output=1;
M=gsua_dmatrix(T,500);
[T7,J7,Y7]=gsua_sa(M,T,'parallel',false,'SensMethod','Saltelli');
save('Results7_sa3.mat','T7','J7','Y7','xdata')
%% 
opt=optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7,res]=gsua_pe(T,xdata,ydata2,'N',100,'opt',opt);
save('Results7_fix2.mat','T7','res','xdata','ydata2')
%% 
th=sum(res<res(1)*1.01)
y7=gsua_eval(T7.Estlsqc(:,1:th),T7,xdata,ydata2); %aca
savefig('familiacurvas.fig')

%% 

%Queremos ver los nuevos casos de contagios
ydata_a= [ydata2(1) diff(ydata2)];
y=y7;
for i=2:size(y7,2)
y(:,i)=y7(:,i)-y7(:,i-1);
end
clf
plot(xdata, ydata_a)
hold on
plot(xdata,y(1,:,1))
legend('Casos nuevos reales', 'estimados')

%% 
T7.Properties.CustomProperties.output=1;
T7.Nominal=T7.Estlsqc(:,1);
T7_2=gsua_ia(T7,T7.Estlsqc,false,true);
%% check
M=gsua_dmatrix(T7,100)
Ua = gsua_ua(M, T7, 'parallel', false, 'ynom',ydata2);

%% 
ynom= gsua_eval(T7.Nominal,T7);
M= gsua_dmatrix(T7,500);
Tsa= gsua_sa(M,T7,'parallel',false, 'SensMethod','Saltelli','ynom',ynom); 
sum(Tsa.Si)/sum(abs(Tsa.Si))
clf 
gsua_plot('Bar',Tsa,Tsa.STi)
%% 
gsua_eval(M',T7)
%% 
Range6= Range;
%% 
beta_h_estimado= T7.Nominal('beta_h');
%Se fijan por las estimaciones
Range(9,:)= beta_h_estimado;
%% 

%vars= [Hit Hi Me He Ms Hr Hs Mi];
Range(6,1)= 0; Range(7,1)= 0;


%% 
[T,~]=gsua_dpmat(odes,vars,[0 80],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output=1;
M=gsua_dmatrix(T,500);
[T7,J7,Y7]=gsua_sa(M,T,'parallel',false,'SensMethod','Saltelli');
save('Results7_sa3.mat','T7','J7','Y7','xdata')
%% 
opt=optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7,res]=gsua_pe(T,xdata,ydata2,'N',100,'opt',opt);
save('Results7_fix2.mat','T7','res','xdata','ydata2')
%% 
th=sum(res<res(1)*1.01)
y7=gsua_eval(T7.Estlsqc(:,1:th),T7,xdata,ydata2); %aca
savefig('familiacurvas.fig')

%% 

%Queremos ver los nuevos casos de contagios
ydata_a= [ydata2(1) diff(ydata2)];
y=y7;
for i=2:size(y7,2)
y(:,i)=y7(:,i)-y7(:,i-1);
end
clf
plot(xdata, ydata_a)
hold on
plot(xdata,y(1,:,1))
legend('Casos nuevos reales', 'estimados')

%% 
T7.Properties.CustomProperties.output=1;
T7.Nominal=T7.Estlsqc(:,1);
T7_2=gsua_ia(T7,T7.Estlsqc,false,true);

%% 

M=gsua_dmatrix(T7,100)
Ua = gsua_ua(M, T7, 'parallel', false, 'ynom',ydata2);





