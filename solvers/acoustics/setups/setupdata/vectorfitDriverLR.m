clear
clc
close all

%Name of output file
fileName = 'LRData6.dat';


%Acoustic constants
rho = 1.2;
c = 343;

%Material properties
sigma = 10900;
dmat = 0.04;

%Frequency range
f_range = [50,2000];

%Amount of poles
Npoles = 6;

%DONT CHANGE ANYTHING FROM HERE
if(mod(Npoles,2)==1)
    display('Npoles must be even');
    return;
end


f = f_range(1):1:f_range(2); % Range where boundary conditions are defined
Zc = rho*c*(1+0.07*(f./sigma).^(-0.632) - 1i*0.107*(f./sigma).^(-0.632));
kt = 2*pi*f/c.*(1+0.109*(f./sigma).^(-0.618) - 1i*0.16*(f./sigma).^(-0.618));

Z = -1i*Zc.*cot(kt*dmat);

Y = 1./Z;

% Fit impedance to a rational function
w = 2*pi*f;
ss = 1i*w;
opts.cmplx_ss = 1;
opts.spy2 = 0;
Ns = length(ss);
weight = ones(1,Ns);

%Complex starting poles :
bet=linspace(w(1),w(Ns),Npoles/2);
poles=[];
for n=1:length(bet)
  alf=-bet(n)*1e-2;
  poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end

Niter=3;
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

