%% References
%B.Gustavsen and A. Semlyen, "Rational approximation of frequency domain
%responses by Vector Fitting", IEEE Trans. Power Delivery, vol. 14, no. 3,
%pp. 1052 - 1061, July 1999. 

%B. Gustavsen, "Improving the pole relocating properties of vector fitting",
%IEEE Trans.Power Delivery, vol 21, no 3, pp 1587-1592, July 2006

%D. Deschrijver, M. Mrozowski, T. Dhaene and D. de Zutter: 
%"Macromodeling of Multiport Systems Using a Fast Implementation of the Vector Fitting Method",
%IEEE Microwave and Wireless Components Letters, vol. 18, no 6, pp 383-385, June 2008

%% Download the Matrix Fitting Toolbox
% https://www.sintef.no/projectweb/vectorfitting/downloads/

%%
clear
clc
close all

fileName = 'freq_dep_lr.dat';

%Acoustic constants
rho = 1.2;
c = 343;

%Material properties
sigma = 50000; % Flow resistivity
dmat = 0.2; % Material thickness
d0 = 0; % Backing air cavity thickness


%Frequency range
f_range = [5,2000];

%Amount of poles
Npoles = 14;

if(mod(Npoles,2)==1)
    print('Npoles must be even');
    return;
end

% Material model
f = f_range(1):1:f_range(2); % Range where boundary conditions are defined
Zc = rho*c*(1+0.07*(f./sigma).^(-0.632) - 1i*0.107*(f./sigma).^(-0.632));
kt = 2*pi*f/c.*(1+0.109*(f./sigma).^(-0.618) - 1i*0.16*(f./sigma).^(-0.618));
k0 = 2*pi*f/c;


% Determine the material surface admittance (frequency response) 
thetai = 0;
kx = sqrt(kt.^2 - k0.^2*(sin(thetai))^2);
if d0 == 0
    Z = -1i*Zc.*(kt./kx).*cot(kx*dmat);
else
    Zd = -1i*rho*c/cos(thetai)*cot(k0*cos(thetai)*d0);
    Z = Zc.*kt./kx.*((-1i*Zd.*cot(kx*dmat) + Zc.*kt./kx)./(Zd - 1i*Zc.*(kt./kx).*cot(kx*dmat)));
end

Y = 1./Z; % Function to be fitted
     


% Fit impedance to a rational function
w = 2*pi*f;
ss = 1i*w;
Ns = length(ss);
weight = ones(1,Ns);

bigY(1,1,:)=Y; %Make 3D array
%================================================
%=           POLE-RESIDUE FITTING               =
%================================================ 
opts.N=Npoles ;%           %Order of approximation. 
opts.poletype='logcmplx'; %Mix of linearly spaced and logarithmically spaced poles
opts.weightparam=2; %5 --> weighting with inverse magnitude norm
opts.Niter1=4;    %Number of iterations for fitting sum of elements (fast!) 
opts.Niter2=4;    %Number of iterations for matrix fitting 
opts.asymp=1;      %Fitting includes D   
opts.logx=0;       %=0 --> Plotting is done using linear abscissa axis 
poles=[];      
[SER,rmserr,bigYfit,opts2]=VFdriver(bigY,ss,poles,opts); %Creating state-space model and pole-residue model 


f2=logspace(-3,5,101); s2=2*pi*i*f2; Ns2=length(s2);
A=diag(SER.A); B=SER.B; C=SER.C; D=SER.D; E=SER.E; I=ones(Npoles,1);
for k=1:Ns2
  Yfit2(k)=D +s2(k)*E+ C*diag((s2(k)*I-A).^(-1))*B; 
  G(k)=real(Yfit2(k));
end    
figure(1); hold on, semilogx(f2,abs(Yfit2),'c-.'); hold off
figure(2001); semilogx(f2,G);

%================================================
%=           Passivity Enforcement              =
%================================================ 
opts.Niter_in=2;
opts.parametertype='Y';
opts.plot.s_pass=2*pi*i*logspace(-3,1,1001).'; 
opts.plot.ylim=[-2e-4 2e-4];
opts.outputlevel=1;

[SER2,bigYfit_passive,opts2]=RPdriver(SER,ss,opts);

SER = SER2;


A = [];
B = [];
C = [];
lambda = [];
al = [];
be = [];
flag = 0;
for i=1:length(SER.C)
    if flag
        flag = 0;
    elseif isreal(SER.C(i))
        A = [A,SER.C(i)];
        lambda = [lambda,-full(SER.A(i,i))];
    else
        flag = 1;
        B = [B,real(SER.C(i))];
        C = [C,imag(SER.C(i))];
        al = [al,-full(real(SER.A(i,i)))]; 
        be = [be,-full(imag(SER.A(i,i)))];
    end
end
Yinf = SER.D;


fileID = fopen(fileName,'w');
fprintf(fileID,'%d %d %d\n',[Npoles,length(A),length(B)]);
if(length(A))
    fprintf(fileID,'%.12f\n',A);
end
if(length(B))
    fprintf(fileID,'%.12f\n',B);
    fprintf(fileID,'%.12f\n',C);
end
if(length(lambda))
    fprintf(fileID,'%.12f\n',lambda);
end
if(length(al))
    fprintf(fileID,'%.12f\n',al);
    fprintf(fileID,'%.12f\n',be);
end
fprintf(fileID,'%.12f\n',Yinf);
fprintf(fileID,'-----\n');
fprintf(fileID,'sigma = %.12f\n',sigma);
fprintf(fileID,'dmat = %.12f\n',dmat);
fprintf(fileID,'freqRange = [%d,%d]\n',f_range);
fclose(fileID);