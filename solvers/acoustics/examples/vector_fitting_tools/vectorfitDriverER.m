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

% Initialize logging file
% Note: one should check the logging file (e.g. ctrl+F "NOT SUCCESSFUL") to
% make sure that the passivity was successfully enforced for all angles


%Name of output file
fileName = 'ERData.dat';

%Acoustic constants
rho = 1.2;
c = 343;

%Material properties
sigma = 3000;
dmat = 0.05;
d0 = 0;

%Frequency range
f_range = [15,900];

%Amount of poles
Npoles = 20;

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
weight=1./abs(Y); %Inverse magnitude wighting

% Complex log spaced starting poles
bet=logspace(log10(w(1)),log10(w(Ns)),Npoles/2);
poles=[];
nu = 1e-3;
for n=1:length(bet)
alf=-nu*bet(n);
poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end


opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=1;      %Include just D, not E in fitting
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=1;   %DO skip identification of residues (C,D,E)
opts.cmplx_ss=1;   %Create complex state space model
opts.spy1=0;       %No plotting for first stage of vector fitting
opts.spy2=0;       %Create magnitude plot for fitting of f(s)
opts.logx=1;       %Use linear abscissa axis
opts.logy=1;       %Use logarithmic ordinate axis
opts.errplot=0;    %Include deviation in magnitude plot
opts.phaseplot=0;  %Do NOT produce plot of phase angle
opts.legend=1;     %Include legends in plots

Niter=4;
for iter=1:Niter
  if iter==Niter, opts.skip_res=0; end
  disp(['   Iter ' num2str(iter)])
  [SER,poles,rmserr,fit]=vectfit3(Y,ss,poles,weight,opts);
  rms(iter,1)=rmserr;
end

figure()
plot(f,real(Y(1,:)))
hold on
plot(f,real(fit(1,:)),'--','linewidth',2)
ylabel('Real(Y)')
legend('Ref.','Fit')

figure()
plot(f,imag(Y(1,:)))
hold on
plot(f,imag(fit(1,:)),'--','linewidth',2)
ylabel('imag(Y)')
legend('Ref.','Fit')

[nierrR,meanerrR,maxerrR,stderrR,nierrI,meanerrI,maxerrI,stderrI] = computeRelativeMappingError(thvec,Y,fit);

fprintf('Material fitted for the frequency range %d to %d Hz\n',f(1),f(end));
fprintf('Rational function error:\n');
fprintf('Real part: Normal incidence error: %1.6f%%, Mean error: %1.6f%%, Max error: %1.6f%%, Std. dev: %1.6f%%\n', nierrR,meanerrR,maxerrR,stderrR);
fprintf('Imaginary part: Normal incidence error: %1.6f%%, Mean error: %1.6f%%, Max error: %1.6f%%, Std. dev: %1.6f%%\n', nierrI,meanerrI,maxerrI,stderrI);




% Do passivity test one angle at a time

SERP.A = SER.A;
SERP.B = SER.B;
SERP.C = zeros(size(SER.C));
SERP.D = zeros(size(SER.D));
SERP.E = zeros(size(SER.E));
fitpassive = zeros(size(fit));



for jj = 1:length(thvec)
    fprintf('Testing angle: %d\n',thvec(jj));
    SERTMP.A = SER.A;
    SERTMP.B = SER.B;
    SERTMP.C = SER.C(jj,:);
    SERTMP.D = SER.D(jj);
    SERTMP.E = SER.E(jj); % Should be zero
    
    [R,a]=ss2pr(SERTMP.A,SERTMP.B,SERTMP.C);
    SERTMP.R=R;
    SERTMP.poles=a;
    
    opts.N=Npoles;%           %Order of approximation. 
    opts.poletype='logcmplx'; %Mix of linearly spaced and logarithmically spaced poles
    opts.weightparam=2; %5 --> weighting with inverse magnitude norm
    opts.Niter1=4;    %Number of iterations for fitting sum of elements (fast!) 
    opts.Niter2=4;    %Number of iterations for matrix fitting 
    opts.asymp=1;      %Fitting includes D   
    opts.logx=0;       %=0 --> Plotting is done using linear abscissa axis 

    opts.Niter_in=2;
    opts.Niter_out=200;
    opts.parametertype='Y';
    opts.plot = [];
    opts.outputlevel = 1;
    [SER2,fitpassive(jj,:),opts2]=RPdriver(SERTMP,ss,opts);
    SERP.C(jj,:) = SER2.C;
    SERP.D(jj) = SER2.D;
    SERP.E(jj) = SER2.E;
    
    %pause;
    
    
end


% Then just work with SERP...
% Need to do some verification on this stuff - e.g. comparing against the
% LR case
% And need to store the passive Yfit and analyze the accuracy...
% And then run some simulations and check that I get similar results when
% using passive vs non-checked (but stable) ER conditions





% Compute the relative error in the rational function approximation
[nierrR,meanerrR,maxerrR,stderrR,nierrI,meanerrI,maxerrI,stderrI,err_real,err_imag] = computeRelativeMappingError(thvec,Y,fitpassive);

fprintf('Material fitted for the frequency range %d to %d Hz\n',f(1),f(end));
fprintf('Rational function error:\n');
fprintf('Real part: Normal incidence error: %1.6f%%, Mean error: %1.6f%%, Max error: %1.6f%%, Std. dev: %1.6f%%\n', nierrR,meanerrR,maxerrR,stderrR);
fprintf('Imaginary part: Normal incidence error: %1.6f%%, Mean error: %1.6f%%, Max error: %1.6f%%, Std. dev: %1.6f%%\n', nierrI,meanerrI,maxerrI,stderrI);


idx = 91;
figure()
plot(f,real(Y(idx,:)),'k--','linewidth',2)
hold on
plot(f,real(fit(idx,:)))
plot(f,real(fitpassive(idx,:)))
legend('reference','before','after')


figure()
plot(f,imag(Y(idx,:)),'k--','linewidth',2)
hold on
plot(f,imag(fit(idx,:)))
plot(f,imag(fitpassive(idx,:)))
legend('reference','before','after')

figure()
plot(err_real)

figure()
plot(err_imag)



% Get the parameters out and in the right format and construct the struct
for jj = 1:length(thvec)
    PPmax = Npoles; % Max number of possible real poles
    Smax = Npoles/2; % Max number of possible complex pole pairs
    
    % Count number of real poles
    PP = 0;
    for i = 1:Npoles
        if isreal(SERP.C(jj,i))
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
        A(i) = SERP.C(jj,i);
        lambda(i) = -full(SERP.A(i,i));
    end
    
    % Retrieve parameters for complex poles
    cntr = 1;
    B = zeros(1,S);
    C = zeros(1,S);
    al = zeros(1,S);
    be = zeros(1,S);
    for i = 1:S2
        if mod(i,2) == 1 % Only look at every other poles, since the come in pairs
            B(cntr) = real(SERP.C(jj,i+PP));
            C(cntr) = imag(SERP.C(jj,i+PP));
            al(cntr) = -full(real(SERP.A(i+PP,i+PP)));
            be(cntr) = -full(imag(SERP.A(i+PP,i+PP)));
            cntr = cntr+1;
        end
    end
    Yinf = SERP.D(jj);
    
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
 No newline at end of file
