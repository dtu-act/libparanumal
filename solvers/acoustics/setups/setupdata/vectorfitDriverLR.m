clear
clc
close all

%Name of output file
fileName = 'LRData.dat';

%Acoustic constants
rho = 1.2;
c = 343;

%Material properties
sigma = 9323.991;
dmat = 0.1;


%Frequency range
f_range = [50,2000];

%Amount of poles
Npoles = 6;

%DONT CHANGE ANYTHING FROM HERE
if(mod(Npoles,2)==1)
    display('Npoles must be even');
    return;
end

opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=1;      %Include just D, not E in fitting
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=0;   %DO skip identification of residues (C,D,E)
opts.cmplx_ss=1;   %Create complex state space model
opts.spy1=0;       %No plotting for first stage of vector fitting
opts.spy2=0;       %Create magnitude plot for fitting of f(s)
opts.logx=1;       %Use linear abscissa axis
opts.logy=1;       %Use logarithmic ordinate axis
opts.errplot=0;    %Include deviation in magnitude plot
opts.phaseplot=0;  %Do NOT produce plot of phase angle
opts.legend=1;     %Include legends in plots

f = f_range(1):1:f_range(2); % Range where boundary conditions are defined
Zc = rho*c*(1+0.07*(f./sigma).^(-0.632) - 1i*0.107*(f./sigma).^(-0.632));
kt = 2*pi*f/c.*(1+0.109*(f./sigma).^(-0.618) - 1i*0.16*(f./sigma).^(-0.618));

Z = -1i*Zc.*cot(kt*dmat);

Y = 1./Z;

% Fit impedance to a rational function
w = 2*pi*f;
ss = 1i*w;
Ns = length(ss);
weight = ones(1,Ns);

%Complex starting poles :
bet=linspace(w(1),w(Ns),Npoles/2);
poles=[];
for n=1:length(bet)
  alf=-bet(n)*1e-2;
  poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end

Niter=10;
for iter=1:Niter
  if iter==Niter, opts.skip_res=0; end
  disp(['   Iter ' num2str(iter)])
  [SER,poles,rmserr,fit]=vectfit3(Y,ss,poles,weight,opts);
  rms(iter,1)=rmserr;
end

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
