%Stephen Guarducci
%MATH 430-001
%UCID: 31221931
%Project
%Synchronization of Strongly Coupled Excitatory Neurons: Relating Network Behavior to Biophysics


clearvars;
close all;

%Ih Parameters
gks = 0;
VL=-65;
gnap = 0.5;
gh = 1.5;
gsyn = 0.013;
gL = 0.5;
I0 = -2.25;
Vhaks = 0;

% %Iks Parameters
% gks = 2.0;
% VL = -54;
% gnap = 0.21;
% gh = 0;
% %gsyn = 0.01; 
% gsyn = 0.015; %we see effects like in Ih 
% gL = 0.1;
% I0 = 1.791;
% Vhaks = -35;
 
%Parameters common to all figures
VNa=55;
VK = -90;
VH = -20;
gNa = 52; 
gK = 11;
C = 1.5;

Vsyn = -20;

if Vsyn > -20
    NT =0.001;
else
    NT = 0;
end

%Time
dt = 0.01;
t = 0:dt:10000;

%Square Wave
ti = 1000;
tf = 3000;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;
Iapp = I0*H;


V = zeros(1,length(t));
mna = zeros(1,length(t));
mnap = zeros(1,length(t));
mks = zeros(1,length(t));
mhf = zeros(1,length(t));
mhs = zeros(1,length(t));
msyn = zeros(1,length(t));
hna = zeros(1,length(t));
n = zeros(1,length(t));

V(1)=-60;
mna(1) = 0.1;
mnap(1) = 0.1;
mks(1) = 0.1;
mhf(1) = 0.1;
mhs(1) = 0.1;
msyn(1) = 0.1;
hna(1) = 0.1;
n(1) = 0.1;


alphamna = @(V) (-0.1*(V+23))./(exp(-0.1*(V+23))-1);
betamna = @(V) 4*exp(-(V+48)./18);

alphahna = @(V) 0.07*exp(-(V+37)/20);
betahna = @(V) 1./(exp(-0.1*(V+7))+1);

alphan = @(V) (-0.01*(V+27))./(exp(-0.1*(V+27))-1);
betan = @(V) 0.125.*exp(-(V+37)/80);

alphamnap = @(V) 1./(0.15*(1+exp(-(V+38)/6.5)));
betamnap = @(V) exp(-(V+38)/6.5)./(0.15*(1+exp(-(V+38)/6.5)));

alphasyn = 1100;
betasyn = 0.19;

%Activation and Inactivation Curves dependent on Voltage
mksinf = @(V) 1./(1+exp(-(V-Vhaks)/6.5));
mhfinf = @(V) 1./(1+exp((V+79.2)/9.78));
mhsinf = @(V) 1./(1+exp((V+71.3)/7.9));
%NOT IN THE PAPER???
mnainf = @(V) alphamna(V)./(alphamna(V)+betamna(V));
mnapinf = @(V) alphamnap(V)./(alphamnap(V)+betamnap(V));
ninf = @(V) alphan(V)./(alphan(V)+betan(V));
hnainf = @(V) alphahna(V)./(alphahna(V)+betahna(V));


%Time Scales dependent on Voltage
Tmks = @(V) 90;
Tmhf = @(V) 0.51./(exp((V-1.7)/10) + exp(-(V+340)/52)+1);
Tmhs = @(V) 5.6./(exp((V-1.7)/14) + exp(-(V+260)/43)+1);
%NOT IN PAPER???
Tmsyn = @(V) 1./(alphasyn+betasyn);
Tn = @(V) 1./(alphan(V)+betan(V));
Thna = @(V) 1./(alphahna(V)+betahna(V));
Tmna = @(V) 1./(alphamna(V)+betamna(V));
Tmnap = @(V) 1./(alphamnap(V)+betamnap(V));

dmksdt = @(V,mks) (mksinf(V)-mks)/Tmks(V);
dmhfdt = @(V,mhf) (mhfinf(V)-mhf)/Tmhf(V);
dmhsdt = @(V,mhs) (mhsinf(V)-mhs)/Tmhs(V);
dmsyndt = @(V,msyn) (alphasyn*NT*(1-msyn)-(betasyn*msyn));
%NOT IN PAPER???
dmnadt = @(V,mna) (alphamna(V)*(1-mna))-(betamna(V)*mna);
dmnapdt = @(V,mnap) (alphamnap(V)*(1-mnap)) -(betamnap(V)*mnap);
dhnadt = @(V,hna) (alphahna(V)*(1-hna))-(betahna(V)*hna);
dndt = @(V,n) (alphan(V)*(1-n))-(betan(V)*n);


for i=1:length(t)-1

  fv1 = (Iapp(i)-((gNa)*mna(i)^3*hna(i)+gnap*mnap(i))*(V(i)-VNa)-(gK*n(i)^(4)+gks*mks(i))*(V(i)-VK)-(gh*(0.65*mhf(i)+0.35*mhs(i))*(V(i)-VH))-(gL*(V(i)-VL))-gsyn*msyn(i)*(V(i)-Vsyn))/C;
  fmna1 = dmnadt(V(i),mna(i));
  fmks1 = dmksdt(V(i),mks(i));
  fmhf1 = dmhfdt(V(i),mhf(i));
  fmhs1 = dmhsdt(V(i),mhs(i));
  fmsyn1 = dmsyndt(V(i), msyn(i));
  fn1 = dndt(V(i),n(i));
  fhna1 = dhnadt(V(i), hna(i));
  fmnap1 = dmnapdt(V(i), mnap(i));
  
  iv = V(i) + fv1*dt;
  imna = mna(i) + fmna1*dt;
  imks = mks(i) + fmks1*dt;
  imhf = mhf(i) + fmhf1*dt;
  imhs = mhs(i) + fmhs1*dt;
  imsyn = msyn(i) + fmsyn1*dt;
  in = n(i) + fn1*dt;
  ihna = hna(i) + fhna1*dt;
  imnap = mnap(i) + fmnap1*dt;
  
  fv2 = (Iapp(i+1)-((gNa)*imna^3*ihna+gnap*imnap)*(iv-VNa)-(gK*in^(4)+gks*imks)*(iv-VK)-(gh*(0.65*imhf+0.35*imhs)*(iv-VH))-(gL*(iv-VL))-gsyn*imsyn*(iv-Vsyn))/C;
  fmna2 = dmnadt(iv,imna);
  fmks2 = dmksdt(iv,imks);
  fmhf2 = dmhfdt(iv, imhf);
  fmhs2 = dmhsdt(iv, imhs);
  fmsyn2 = dmsyndt(iv, imsyn);
  fn2 = dndt(iv,in);
  fhna2 = dhnadt(iv,ihna);
  fmnap2 = dmnapdt(iv, imnap);
    
  V(i+1) = V(i) + (fv1+fv2)*dt/2;
  mna(i+1) = mna(i) + (fmna1+fmna2)*dt/2;
  mks(i+1) = mks(i) + (fmks1+fmks2)*dt/2;
  mhf(i+1) = mhf(i) + (fmhf1+fmhf2)*dt/2;
  mhs(i+1) = mhs(i) + (fmhs1+fmhs2)*dt/2;
  msyn(i+1) = msyn(i) + (fmsyn1+fmsyn2)*dt/2;
  n(i+1) = n(i) + (fn1+fn2)*dt/2;
  hna(i+1) = hna(i) + (fhna1+fhna2)*dt/2;
  mnap(i+1) = mnap(i) + (fmnap1 + fmnap2)*dt/2;
  
end  
 
%Voltage vs. Time
figure
plot(t,V,'b','linewidth',2);
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('gsyn = 0.013');
set(gca,'fontsize',20);


%Current vs. Time
figure
plot(t,Iapp,'b','linewidth',2);
set(gca,'fontsize',20);
ylabel('Current(I)');
xlabel('Time(ms)')
title('Current vs Time');
 

V = -100:0.1:100;
%Voltage vs Activation and Inactivation Curves
figure
hold on
plot(V,mksinf(V),'g','linewidth',2);
plot(V,mhfinf(V),'b','linewidth',2);
plot(V,mhsinf(V),'k','linewidth',2);
plot(V,mnainf(V),'r','linewidth',2);
plot(V,mnapinf(V),'y','linewidth',2);
plot(V,ninf(V),'c','linewidth',2);
plot(V,hnainf(V),'m','linewidth',2);
set(gca,'fontsize',20);
xlabel('Voltage');
legend('mks', 'mhf', 'mhs', 'mna', 'mnap', 'n','hna');


%Voltage vs Time Scales
figure
hold on
plot([-100 100],[90 90],'g','linewidth',2);
plot(V,Tmhf(V),'b','linewidth',2);
plot(V,Tmhs(V),'k','linewidth',2);
plot([-100 100],[0.0009089339114 0.0009089339114],'r','linewidth',2);
plot(V,Tmna(V),'y','linewidth',2);
plot(V,Tmnap(V),'c','linewidth',2);
plot(V,Tn(V),'m','linewidth',2);
plot(V,Thna(V),'g','linewidth',2);
set(gca,'fontsize',20);
xlabel('Voltage');
legend('Tmks', 'Tmhf', 'Tmhs','Tmsyn', 'Tmna','Tmnap', 'Tn','Thna');
