clear
clc
close all

%Name of output file
fileName = 'ERDATA14.dat';

%Acoustic constants
rho = 1.2;
c = 343;

%Material properties
%sigma = 47700;
%dmat = 0.05;
%d0 = 0.15;
sigma = 9.323991e3; %2
%sigma = 1.062646e4;%1


dmat = 0.1;
d0 = 0;

%Frequency range
f_range = [50,2000];

%Amount of poles
Npoles = 14;

%DONT CHANGE ANYTHING FROM HERE
if(mod(Npoles,2)==1)
    display('Npoles must be even')
    return;
end

f = f_range(1):1:f_range(2); % Range where boundary conditions are defined
Zc = rho*c*(1+0.07*(f./sigma).^(-0.632) - 1i*0.107*(f./sigma).^(-0.632));
kt = 2*pi*f/c.*(1+0.109*(f./sigma).^(-0.618) - 1i*0.16*(f./sigma).^(-0.618));
k0 = 2*pi*f/c;

thvec = 0:90;


Y = zeros(length(thvec),length(f));
Z = zeros(length(thvec),length(f));


for jj = 1:length(thvec)
    % Determine the material surface admittance (frequency response) for
    % the different angles considered
    thetai = thvec(jj)*pi/180;
    kx = sqrt(kt.^2 - k0.^2*(sin(thetai))^2);
    if d0 == 0
        Z(jj,:) = -1i*Zc.*(kt./kx).*cot(kx*dmat);
    else
        Zd = -1i*rho*c/cos(thetai)*cot(k0*cos(thetai)*d0);
        Z(jj,:) = Zc.*kt./kx.*((-1i*Zd.*cot(kx*dmat) + Zc.*kt./kx)./(Zd - 1i*Zc.*(kt./kx).*cot(kx*dmat)));
    end
    
    Y(jj,:) = 1./Z(jj,:);
    
    
end

% Fit the admittance to a rational function with vector fitting
w = 2*pi*f;
ss = 1i*w;
Ns = length(ss);

%Complex starting poles :
bet=linspace(w(1),w(Ns),Npoles/2);
poles=[];
for n=1:length(bet)
    alf=-bet(n)*1e-2;
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ];
end
% Curve fitting parameters
weight=1./abs(Y); %Inverse magnitude wighting
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

%=============================================
disp('vector fitting...')
Niter=10;
for iter=1:Niter
    if iter==Niter, opts.skip_res=0; end
    disp(['   Iter ' num2str(iter)])
    [SER,poles,rmserr,Yfit]=vectfit3(Y,ss,poles,weight,opts);
    rms(iter,1)=rmserr
    
end

% Get the parameters out and in the right format and construct the struct
for jj = 1:length(thvec)
    PPmax = Npoles; % Max number of possible real poles
    Smax = Npoles/2; % Max number of possible complex pole pairs
    
    % Count number of real poles
    PP = 0;
    for i = 1:Npoles
        if isreal(SER.C(jj,i))
            PP = PP+1;
        end
    end
    
    S2 = Npoles-PP; % number of complex poles (nb number of complex pole pairs S = S2/2)
    S = S2/2;
    
    %fprintf('Angle: %d deg, # real poles = %d, # complex pole pairs = %d\n********\n',thvec(jj),PP,S);
    
    
    % I am assuming that the vectfit3 returns first the real poles and the the
    % complex poles...
    A = zeros(1,PP);
    lambda = zeros(1,PP);
    % Retrieve parameters for real poles
    for i = 1:PP
        A(i) = SER.C(jj,i);
        lambda(i) = -full(SER.A(i,i));
    end
    
    % Retrieve parameters for complex poles
    cntr = 1;
    B = zeros(1,S);
    C = zeros(1,S);
    al = zeros(1,S);
    be = zeros(1,S);
    for i = 1:S2
        if mod(i,2) == 1 % Only look at every other poles, since the come in pairs
            B(cntr) = real(SER.C(jj,i+PP));
            C(cntr) = imag(SER.C(jj,i+PP));
            al(cntr) = -full(real(SER.A(i+PP,i+PP)));
            be(cntr) = -full(imag(SER.A(i+PP,i+PP)));
            cntr = cntr+1;
        end
    end
    Yinf = SER.D(jj);
    
    % Set up structs for the different angles
    BCdata(jj) = struct('A',A,'B',B,'C',C,'Yinf',Yinf);
end



fileID = fopen(fileName,'w');
fprintf(fileID,'%d %d %d\n',[Npoles,PP,S]);
if(length(A))
    for i = 1:length(thvec)
        fprintf(fileID,'%.12f\n',BCdata(i).A);
    end
end
if(length(B))
    for i = 1:length(thvec)
        fprintf(fileID,'%.12f\n',BCdata(i).B);
    end
    for i = 1:length(thvec)
        fprintf(fileID,'%.12f\n',BCdata(i).C);
    end
end
if(length(lambda))
    fprintf(fileID,'%.12f\n',lambda);
end
if(length(al))
    fprintf(fileID,'%.12f\n',al);
    fprintf(fileID,'%.12f\n',be);
end
for i = 1:length(thvec)
    fprintf(fileID,'%.12f\n',BCdata(i).Yinf);
end
fprintf(fileID,'-----\n');
fprintf(fileID,'sigma = %.12f\n',sigma);
fprintf(fileID,'dmat = %.12f\n',dmat);
fprintf(fileID,'d0 = %.12f\n',d0);
fprintf(fileID,'freqRange = [%d,%d]\n',f_range);
fclose(fileID);

%length(BCdata(1).A)*91+length(BCdata(1).B)*91+length(BCdata(1).C)*91+length(lambda)+length(al)+length(be)+1+91