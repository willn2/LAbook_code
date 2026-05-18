% MatlabCodeLABook  3.5.2026 mfi
% A listing of Matlab scripts used in the LABook
%
% Solution to Exercise E8.6 generating results in Figure 8.3 
%
close all;t=0:32;                                   % time in seconds
% generate two impulse responses, a,b and implement as circulant matrices
a=0.2*[[7 -2.5 -1 -0.5 -0.25 -0.1] zeros(1,22) [-0.1 -0.25 -0.5 -1 -2.5]];
b=0.7*[1 0.5 0.1 zeros(1,27) 0.1 0.5 1.0];
N=length(b);icnt=0;
A(1,:)=a; B(1,:)=b;                                 % transformation matrices...
for j=3:2:N                                         % ...in circulant form
     icnt=(j+1)/2;                                  % downsample output by 2
     A(icnt,1:j-1)=a(N-j+2:N);
     A(icnt,j:end)=a(1:N-j+1);
     B(icnt,1:j-1)=b(N-j+2:N);
     B(icnt,j:end)=b(1:N-j+1);
end
x=zeros(size(t));x(2:2:10)=1;x(16:20)=1;x(26:30)=1; % input function
rng(0);n=0.2*randn(1,17);                           % measurement noise
figure(1)
subplot(3,2,1),stem(t,A(9,:),'k','linewidth',1),ylim([-1 2])
xlabel('t'),ylabel('A(9,:)'),set(gca,'FontSize',14)
subplot(3,2,2),stem(t,B(9,:),'k','linewidth',1),ylim([-1 2])
xlabel('t'),ylabel('B(9,:)'),set(gca,'FontSize',14)
subplot(3,2,3),imagesc(A);colormap gray,axis equal,title('A')
subplot(3,2,4),imagesc(B);colormap gray,axis equal,title('B')
subplot(3,2,5),plot(t,x,'k','linewidth',2),ylim([-1 2])
subplot(3,2,6),plot(t,x,'k','linewidth',2),ylim([-1 2])
%
u=0:1/33:1/2;                     % frequency axis from time axis param's
FA=abs(fft(a));FB=abs(fft(b));    % impulse response spectra
figure,plot(u,FA(1:icnt),'r'),hold on,plot(u,FB(1:icnt),'b'),hold off
xlabel('u'),ylabel('abs(fft)')
y1A=A*x';y1B=B*x';y2A=y1A+n';y2B=y1B+n'; % measurements w & w/o noise
figure,subplot(2,2,1),plot(1:icnt,y1A),xlabel('t'),ylabel('y1A')
subplot(2,2,2),plot(1:icnt,y1B),xlabel('t'),ylabel('y1B')
subplot(2,2,3),plot(1:icnt,y2A),xlabel('t'),ylabel('y2A')
subplot(2,2,4),plot(1:icnt,y2B),xlabel('t'),ylabel('y2B')
figure(1)
subplot(3,2,5),hold on,plot(t,pinv(A)*y1A,'r','LineWidth',1) % est x
plot(t,pinv(B)*y1B,'b','LineWidth',1)               % estimate x from y
subplot(3,2,6),hold on,plot(t,pinv(A)*y2A,'r','LineWidth',1)
plot(t,pinv(B)*y2B,'b','LineWidth',1)
%
%
%%
% Solution to Exercise E8.7 resulting in Figure 8.4.
%
% b is beta, the coefficient vector used to generate yp in y1 and y2.
% b1 and b2 are estimates of b made from analyzing noisy data.
%
t1=0:0.1:20;b=[1000 -20 -8 0.5];                  % 20s time trace, b=beta
R1=zeros(1,6);R2=R1;                              % initialization
n1=200*randn(size(t1));                           % noise 1 (noise variance = 200^2)
yp = b(1) + b(2)*t1 + b(3)*t1.^2 + b(4)*t1.^3;    % signal 1
SNR=10*log10(sum(yp.^2)/t1(end)/200^2);           % SNR in dB
y1=yp + n1;                                       % measurement 1
subplot(1,3,1),plot(t1,y1,'k.','linewidth',2),hold on
set(gca,'FontSize',14),xlabel('t'),ylabel('Signal Amplitude')
t2=0:0.05:10;                                     % 10s time trace
n2=200*randn(size(t2));                           % noise 2 (same distribution)
y2 = b(1) + b(2)*t2 + b(3)*t2.^2 + b(4)*t2.^3 + n2;  % measurements 2
subplot(1,3,2),plot(t2,y2,'k.','linewidth',2),hold on
set(gca,'FontSize',14),xlabel('t'),ylabel('Signal Amplitude')
%
T1=[ones(length(t1),1) t1' t1.^2' t1.^3' t1.^4' t1.^5']; % time matrix 1
T2=[ones(length(t2),1) t2' t2.^2' t2.^3' t2.^4' t2.^5']; % time matrix 2
for j=1:6       % fit data with polynomial models of order 0-5
    b1=pinv(T1(:,1:j))*y1';             % estimate beta 1
    if j==4; bb1=b1; end                % save coeff from N=3 model
    z1=T1(:,1:j)*b1;
    if j==2||j==4||j==6                 % plot only these three models
        subplot(1,3,1), plot(t1,z1)
    end
    R1(j)=rms(z1-y1');                  % RMS 1 values
    b2=pinv(T2(:,1:j))*y2';             % estimate beta 2
    if j==4; bb2=b2;end
        z2=T2(:,1:j)*b2;                % plot three models
    if j==2||j==4||j==6
        subplot(1,3,2), plot(t2,z2)
    end
    R2(j)=rms(z2-y2');                  % RMS 2 values
end
subplot(1,3,3),plot(0:5,R1,'ko-','linewidth',2),hold on,
plot(0:5,R2,'kx-','linewidth',2),hold off
set(gca,'FontSize',14),xlabel('polynomial order'),ylabel('RMSE')
%
%
%%
% Solution to Exercise E9.1.
%
V=[1 2 1;1 2 2;0.5 2 3];        % original vectors
plot3([0 V(1,1)],[0 V(2,1)],[0 V(3,1)],'k','linewidth',3),hold on
plot3([0 V(1,2)],[0 V(2,2)],[0 V(3,2)],'r','linewidth',3)
plot3([0 V(1,3)],[0 V(2,3)],[0 V(3,3)],'b','linewidth',3)
xlabel('x'),ylabel('y'),zlabel('z'),set(gca,'FontSize',18)
xlim([0 3]),ylim([0 3]),zlim([0 3])
text(V(1,1)+0.05,V(2,1)+0.05,V(3,1)+0.05,'v_1','fontsize',18)
text(V(1,2)+0.05,V(2,2)+0.05,V(3,2)+0.05,'v_2','fontsize',18)
text(V(1,3)+0.05,V(2,3)+0.05,V(3,3)+0.05,'v_3','fontsize',18)
grid
figure
[U,~]=MGSO(V);                  % apply the Gram-Schmidt algorithm
plot3([0 U(1,1)],[0 U(2,1)],[0 U(3,1)],'k','linewidth',3),hold on
plot3([0 U(1,2)],[0 U(2,2)],[0 U(3,2)],'r','linewidth',3)
plot3([0 U(1,3)],[0 U(2,3)],[0 U(3,3)],'b','linewidth',3)
xlabel('x'),ylabel('y'),zlabel('z'),set(gca,'FontSize',18)
xlim([-1 1]),ylim([-1 1]),zlim([-1 1])
text(U(1,1)+0.05,U(2,1)+0.05,U(3,1)+0.05,'u_1','fontsize',18)
text(U(1,2)+0.05,U(2,2)+0.05,U(3,2)+0.05,'u_2','fontsize',18)
text(U(1,3)+0.05,U(2,3)+0.05,U(3,3)+0.05,'u_3','fontsize',18)
grid
%
%
%
% function [Q,R] = MGSO(A)
% %
% % The Modified Gram-Schmitt Orthonormalization via
% % See http://web.mit.edu/18.06/www/Essays/gramschmidtmat.pdf
% % A is an input matrix of column vectors to be analyzed
% % Q is the orthonormalized output matrix
% %
% [M,N]=size(A);                  % Preset the array spaces
% Q=zeros(M,N);R=zeros(N);        % same
% for j=1:N
%     v=A(:,j);                   % v is the j column of A
%     for i=1:j-1
%         R(i,j)=Q(:,i)'*A(:,j);  % modify to improve accuracy
%         v=v-R(i,j)*Q(:,i);      % subtract projection (q_i^t*v)*q_i
%     end                         % v is now perp to all q_1,...,q_(j-1)
%     R(j,j)=norm(v);
%     Q(:,j)=v/R(j,j);            % normalize v to be unit vector q_j
% end                             % Q is orthonormal version of A
% end
% %
%
%%
% Solution to Example 10.5 generating results in Figure 10.3.
%
n=0:31;h=zeros(1,32);h(1:9)=1-n(1:9)/8;h(25:32)=n(1:8)/8; % h = impulse response
H=toeplitz([h(1) fliplr(h(2:end))],h);N=length(n);      % H = system matrix
u=0:1/N:1/2-1/(2*N);                                    % frequency axis
imagesc(n,n,H),colormap gray                            % view circulant system matrix
for k=1:N                                               % generate Fourier basis from N=32 and T=1
    for j=1:N
        Q(k,j)=exp(1i*2*pi*(j-1)*(k-1)/N);
    end
end
F=Q'*H*Q/N;HHH=abs(diag(F));                            % eigen-decompose H with Q
Fhp=abs(fft(h));Fh=Fhp(1:N/2);                          % fft comparison
figure,subplot(2,1,1),stem(n,h,'k','linewidth',1.5)     % view impulse response
set(gca,'FontSize',14),xlabel('n'),ylabel('h[n]')
subplot(2,1,2),plot(u,HHH(1:N/2),'k',LineWidth=1.5)     % view eigen-decomposed H result
hold on,plot(u,Fh,'o','linewidth',3)                    % view fft result
set(gca,'FontSize',14),xlabel('k'),ylabel('DFT')
%
for j=1:16; x(j)=(sin(8*pi*u(j)))^2/(8*(sin(pi*u(j)))^2);end
plot(u,x,'gx','linewidth',1.5)                          % view analytic result
legend('eigenanalysis','fft','analytic')
%
%
%%
% Solution to Exercise E10.6 with results plotted in Figure 10.5a.
%
dt=0.1;t=dt:dt:6;               % initialize
% A=0.1*(fix(16*rand(4))-8);    % for exploration only; next line has A
A=[-0.3 0.7 -0.4 -0.2;-0.3 0.6 -0.2 0.7;-0.1 -0.8 0 -0.4;-0.3 0.3 0.7 0.3];
y0=[1 0 0 0]';                  % initial conditions
[U,D] = eig(A);c=U\y0;          % compute eigensystem
y=c(1)*U(:,1)*exp(D(1,1)*t) + c(2)*U(:,2)*exp(D(2,2)*t) ...
    + c(3)*U(:,3)*exp(D(3,3)*t) + c(4)*U(:,4)*exp(D(4,4)*t); % compute solution
figure,hold on
plot(t,real(y(1,:)),'bo-','linewidth',2),plot(t,real(y(2,:)),'ro-','linewidth',2)
plot(t,real(y(3,:)),'ko-','linewidth',2),plot(t,real(y(4,:)),'mo-','linewidth',2)
set(gca,'FontSize',14),xlabel('t')
legend('y_1','y_2','y_3','y_4','location','northwest')
%
%
%%
% Solution to Exercise E10.7 with results plotted in Figure 10.5b and 10c.
%
dx=0.1;x=-10:dx:1;N=length(x);a=[0,0.1,0.2,0.5];  % initialize
y=zeros(2,N);cx=zeros(1,N);                       % initialize
y0=[0;0.5];                                       % initial conditions
for k=1:4
    for j=1:N
        A=[0 1;0.5*x(j) a(k)];                    % matrix A(t;a)
        cx(j)=cond(A);                            % find condition of A
        [U,D] = eig(A);c=U\y0;          % compute eigensystem, constants
        if j==25;DD(1,:)=[D(1,1) D(2,2)];end      % record some eigenvalues
        if j==50;DD(2,:)=[D(1,1) D(2,2)];end
        if j==75;DD(3,:)=[D(1,1) D(2,2)];end
        if j==100;DD(4,:)=[D(1,1) D(2,2)];end
        y1(k,j)=c(1)*U(1,1)*exp(D(1,1)*x(j)) ...
            + c(2)*U(1,2)*exp(D(2,2)*x(j));       % compute y(t;a) (4 plots)
    end
end
plot(x,real(y1),'linewidth',1.5)                  % y(x) for 4 values of a
set(gca,'FontSize',14),xlabel('x'),ylabel('y(x;a)')
legend('a=0.0','a=0.1','a=0.2','a=0.5','location','northwest')
figure,semilogy(x,cx,'k','linewidth',2)           % plot condition(A)
set(gca,'FontSize',18),xlabel('x'),ylabel('condition(A)'),ylim([1 100])
%
%
%%
% Solution to Exercise E10.8 with results displayed in Figure 10.6.
%
sx=1;sy=2;r=1/sqrt(2);          % sx,sy are std in x and y.  r=rho
K=[sx^2 r*sy*sx;r*sx*sy sy^2];  % compute covariance matrix K from parts
[U,D]=eig(K);                       % find eigensystem of K
[~,ind] = sort(diag(D),'descend');  % reorder in descending order
Dp = D(ind,ind);                    % reordered eigenvalues
Up = U(:,ind);                      % reordered eigenvectors
DD1=sqrt(Dp(1,1));DD2=sqrt(Dp(2,2)); % compute the singular values
plot([-DD1*Up(1,1) DD1*Up(1,1)], ...
    [-DD1*Up(2,1) DD1*Up(2,1)],'r','linewidth',2)
hold on                             % plot eigenvectors
plot([-DD2*Up(1,2) DD2*Up(1,2)], ...
    [-DD2*Up(2,2) DD2*Up(2,2)],'b','linewidth',2)
xlabel('x'),ylabel('y'),set(gca,'FontSize',16)
xlim([-4 4]),ylim([-4 4]),axis square
%
% Compute y for Mahalanobis distance = 1 to visualize K
%
x=-4:0.1:5;N=length(x);y1=zeros(1,N);y2=y1;   % find an arbitrary x-axis range
for j=1:N                           % compute y using quadratic equation
    b=2*r*x(j)*sy/sx;c=(sy*x(j)/sx)^2-sy^2*sqrt(1-r^2);
    y1(j)=b/2+real(sqrt(b^2-4*c)/2);
    y2(j)=b/2-real(sqrt(b^2-4*c)/2);
end
plot(x,y1,x,y2)
m=zeros(70,2);X=mvnrnd(m,K);hold on;plot(X(:,1),X(:,2),'ko')
%
%
%%
% Solution to Exercise E11.2 with results plotted in Figure 11.1.
% First, simulate the measurement data
%
M=500;mx=40;my=60;                              % initialize parameters
m=ones(M,2);m(:,1)=m(:,1)*mx;m(:,2)=m(:,2)*my;  % initialize parameters
sx=10;sy=20;r=0.9;                              % initialize parameters
K=[sx^2 r*sy*sx;r*sx*sy sy^2];                  % population cov matrix
X=mvnrnd(m,K,M);                                % draw 500 samples of x,y
Xb=mean(X);Kx=cov(X);                           % sample mean, cov matrix
%
%  Method (1) Use Eigenanalysis to compute PCs from sample statistics
%
[U,D]=eig(Kx);                                  % eigensystem of Kx
[d,ind] = sort(diag(D),'descend');              % place in descending order
Dp = D(ind,ind);DD1=sqrt(Dp(1,1));DD2=sqrt(Dp(2,2));
Up = U(:,ind);          % arrange eigenvectors for plotting the PCs
plot(X(:,1),X(:,2),'ko'),axis square,hold on    % plot the M data pairs as points
plot([-50*Up(1,1) 50*Up(1,1)]+mx,[-50*Up(2,1) 50*Up(2,1)]+my,'r','LineWidth',2)
plot([-10*Up(1,2) 10*Up(1,2)]+mx,[-10*Up(2,2) 10*Up(2,2)]+my,'b','LineWidth',2)
xlim([0 100]),ylim([0 100]),set(gca,'FontSize',14)
xlabel('x'),ylabel('y')
% Y=X*Up;                                       % transform data (not used)
% figure,plot(Y(:,1),Y(:,2),'kx'),axis square
%
% Method (2) Use SVD to compute PCs
%
[Vs,Ds,Us]=svd(X-Xb);                           % To use SVD correctly, subtract mean...
figure,plot(X(:,1),X(:,2),'ko'),axis square,hold on %...but add the mean back in & display
plot([-50*Us(1,1) 50*Us(1,1)]+mx,[-50*Us(2,1) 50*Us(2,1)]+my,'r','LineWidth',2)
plot([-10*Us(1,2) 10*Us(1,2)]+mx,[-10*Us(2,2) 10*Us(2,2)]+my,'b','LineWidth',2)
xlim([0 100]),ylim([0 100]),set(gca,'FontSize',14)
xlabel('x'),ylabel('y')
%
%  Method (3) Use Matlab's PCA function to compute PCs
%
C=pca(X);                                       % use Matlab's pca...
figure, plot(X(:,1),X(:,2),'ko'),axis square    % ...to find same result
hold on,plot([-50*C(1,1) 50*C(1,1)]+mx,[-50*C(2,1) 50*C(2,1)]+my,'r','LineWidth',2)
plot([-10*C(1,2) 10*C(1,2)]+mx,[-10*C(2,2) 10*C(2,2)]+my,'b','LineWidth',2)
xlim([0 100]),ylim([0 100]),set(gca,'FontSize',14)
xlabel('x'),ylabel('y')
%
%
%%
% Solution to Exercise E11.3, Part (a). Results are plotted in Figure 11.2.
% Compute H1, four examples of sample output Y, and the MTF curve for an LTI system.
%
M=128;tt=1:M;h=zeros(size(tt));             % tt is time axis, initialize h
h(1:M/2+1)=exp(-(tt(1:M/2+1)-1).^2/(2*5^2))/(5*sqrt(2*pi)); % half of h 
h(M:-1:M/2+1)=h(2:M/2+1);                   % other half of h w/o 1st value
H1=toeplitz(h);                             % circulant system matrix
imagesc(H1),colormap gray, axis square, title('Time-invariant H')
%
T0=[8 16 32 48];                            % input sinusoidal periods: u0=1/T0
for j=1:4                                   % input four sinusoids
    X1(j,:)=sin(2*pi*tt/T0(j));             % input matrix X1 is 4x128
end
Y1=H1*X1';                                  % output matrix Y1 is 128x4
figure
for j=1:4                                   % plot I/O pairs in 2x2 plots
    subplot(2,2,j),plot(tt,X1(j,:),'k--','linewidth',1),hold on % input
    plot(1:M,Y1(:,j),'k','linewidth',2)     % output 
    xlabel('t (samples)'),ylabel('x(t) or y(t)'),set(gca,'fontsize',12)
    title(['T_0 = ',(num2str(T0(j)))])      % label period of the sine wave
end
u0=[0 1./T0(4:-1:1)];mtf=[(max(Y1)) 1];     % measured mtf points versus frequency
figure,plot(u0,mtf(5:-1:1),'ko','linewidth',2),hold on % flip so freq axis increases
uu=0:0.01:0.15;MTF=exp(-2*pi^2*5^2*uu.^2);  % predicted MTF to compare with mtf
plot(uu,MTF,'k','LineWidth',2)
xlabel('u_0 (Hz)'),ylabel('MTF'),set(gca,'fontsize',14)
%
%
%%
% Solution to Exercise E11.3, Part (b). Results are plotted in Figure 11.3.
% Compute H, four examples of output for LTV measurements, and singular spectra.
% Note:  You need H1 from the previous script to compare results in last plot.
%
% close all
M=128;N=200;t=1:N;s1=7;ds=2;    % H is MxN, s1,ds define sigma range
for j=1:M/2                     % generate time-dependent sigma for h
    s(j)=s1-ds*j/(M/2);         % first half (s1 is start, ds is change)
    s(M-j+1)=s(j);              % second half
end
H=zeros(M,N);                   % initialize H (system matrix)
for j=1:M                       % time-varying Gaussian H: MxN
    H(j,:)=exp(-(t-((j-1)+M/4)).^2/(2*s(j)^2))/(s(j)*sqrt(2*pi));
end
figure,imagesc(H),colormap gray,axis equal,title('Time-varying H')
%
T0=[8 16 32 48];                % input sinusoidal periods: u0=1/T0
for j=1:4                       % input four sinusoids
    X(j,:)=sin(2*pi*t/T0(j));   % input matrix X is Nx4 
end
Y=H*X';                         % output matrix Y is Mx4
figure
for j=1:4                       % plot I/O pairs in 2x2 plots
    subplot(2,2,j),plot(1:M,X(j,1:M),'k--','linewidth',1),hold on % input
    plot(1:M,Y(:,j),'k','linewidth',2)                            % output
    xlabel('t (samples)'),ylabel('x(t) or y(t)'),set(gca,'fontsize',12)
    title(['T_0 = ',(num2str(T0(j)))])
end
%
[~,D,~]=svd(H);DD=diag(D);      % SVD of LSV system matrix (singular values only)
figure,plot(1:50,DD(1:50),'r','linewidth',2),hold on
[~,D1,~]=svd(H1);DD1=diag(D1);  % SVD of LTI system matrix
plot([1 2:2:51],[1;DD1(2:2:51)],'k','linewidth',2),hold on      % discard every other point
xlabel('singular value'),ylabel('singular spectrum'),set(gca,'fontsize',14)
legend('LTV','LTI')
%
%
%%
% Solution to Example 13.1 with results displayed in Figure 13.3b.
%
S=[1 -1 -1 0 0 0;0 1 0 1 -1 0;0 0 1 -1 0 0;0 0 0 0 1 -1]; % stoichiometric matrix
T=0.01;T0=20;t=0:T:T0;N=length(t);                        % time base where T<<T0
X=zeros(4,N);               % rows are conc's and columns are time pts
x0=[0 0 0 0]';X(:,1)=x0;    % concentration of 4 metabolites, initialized
k=[1 1 1 1 1 1];            % rate constants: forward fluxes
kn=[0 0 1 1 0 0];           % rate constants: reverse fluxes
for j=2:N                   % j is time index
    f=[k(1);k(2)*X(1,j-1);k(3)*X(1,j-1)-kn(3)*X(3,j-1);...  % flux vector 
        k(4)*X(3,j-1)-kn(4)*X(2,j-1);k(5)*X(2,j-1);k(6)*X(4,j-1)];
    dx=S*f*T;               % incremental change in molecular concentration
    X(:,j)=X(:,j-1)+dx;     % concentration at time t(j) 
end
plot(t,X,'LineWidth',1.5) % plot the four concentration time series
xlabel('Time'),ylabel('Concentration'),set(gca,'fontsize',14) 
legend('x_1','x_2','x_3','x_4')
%
%
%%
% Solution to Exercise E13.2 with results plotted in Figure 13.6.
%
S=[1 -1 -1 0 0 0;0 1 0 1 -1 0;0 0 1 -1 0 0;0 0 0 0 1 -1]; % stoichiometric matrix
T=0.01;T0=40;t=0:T:T0;N=length(t);Np=(N-1)/2;    % time base where T<<T0
X=zeros(4,N);               % rows are conc's and columns are time pts
x0=[0 0 0 0]';X(:,1)=x0;    % concentration of 4 metabolites, initialized
k=[1 1 1 1 1 1];            % rate constants: forward fluxes
kn=[0 0 1 1 0 0];           % rate constants: reverse fluxes
for j=2:Np                  % j is time index
    f=[k(1);k(2)*X(1,j-1);k(3)*X(1,j-1)-kn(3)*X(3,j-1);...  % flux vector 
        k(4)*X(3,j-1)-kn(4)*X(2,j-1);k(5)*X(2,j-1);k(6)*X(4,j-1)];
    dx=S*f*T;               % incremental change in molecular concentration
    X(:,j)=X(:,j-1)+dx;     % concentration at time t(j) 
end
X(:,Np+1)=X(:,Np);          % initialize conc for 20<t<40
X(1,Np+1)=1.2*X(1,Np);X(2,Np+1)=0.8*X(2,Np);
for j=Np+2:N                % j is time index
    f=[k(1);k(2)*X(1,j-1);k(3)*X(1,j-1)-kn(3)*X(3,j-1);...  % flux vector 
        k(4)*X(3,j-1)-kn(4)*X(2,j-1);k(5)*X(2,j-1);k(6)*X(4,j-1)];
    dx=S*f*T;               % incremental change in molecular concentration
    X(:,j)=X(:,j-1)+dx;     % concentration at time t(j) 
end
hold on,plot(t,X,'LineWidth',1.5) % plot the four concentration time series
xlabel('Time'),ylabel('Concentration'),set(gca,'fontsize',14) 
legend('x_1','x_2','x_3','x_4')
%
%
%%
% Solution to Example 14.1. Results plotted in Figure 14.1. M = M1+M2 = 200.  m1x is population mean along
% x axis for class 1, s2y^2 is population variance along the y axis for class 2, etc.
%
M1=100;M2=M1;M=M1+M2;m1x=40;m1y=60;m2x=60;m2y=40;       % initialize parameters
m1=ones(M1,2);m1(:,1)=m1(:,1)*m1x;m1(:,2)=m1(:,2)*m1y;  % form vectors of mean values
m2=ones(M2,2);m2(:,1)=m2(:,1)*m2x;m2(:,2)=m2(:,2)*m2y;  % r1,r2: corr coeff btw x&y, class 1,2
s1x=10;s1y=20;r1=0.9;s2x=20;s2y=10;r2=-0.6;             % other covariance matrix parameters
K1=[s1x^2 r1*s1y*s1x;r1*s1x*s1y s1y^2];                 % population cov matrix
K2=[s2x^2 r2*s2y*s2x;r2*s2x*s2y s2y^2];                 % population cov matrix
X1=mvnrnd(m1,K1,M1);X2=mvnrnd(m2,K2,M2);                % simulate two classes of data, X1 and X2
X1b=mean(X1);K1x=cov(X1);X2b=mean(X2);K2x=cov(X2);      % sample means and cov matrices
plot(X1(:,1),X1(:,2),'ko'),axis square,hold on,plot(m1x,m1y,'kx','linewidth',3)
plot(X2(:,1),X2(:,2),'ro'),plot(m2x,m2y,'rx','linewidth',3)
xlim([-50 150]),ylim([-50 150]),set(gca,'FontSize',18)
xlabel('x'),ylabel('y')
a=2*(X2b-X1b)/((K1x+K2x)/2);            % linear discriminant parameter a
b=X2b/((K1x+K2x)/2)*X2b'-X1b/((K1x+K2x)/2)*X1b'; % lin disc parameter b
x=-50:150;                              % x axis range 
y=(b-a(1)*x)/a(2);                      % generate linear discriminant line 
plot(x,y,'k--',LineWidth=2)
e1=0;e2=e1;                             % initialize classification error counters
for j=1:M1                              % if M1 is not = M2 adjust M1 in for loop
    yy1=(b-a(1)*X1(j,1))/a(2);          % find LD function for each value of x
    yy2=(b-a(1)*X2(j,1))/a(2);
    if X1(j,2)<yy1;e1=e1+1;end          % index error count if below line
    if X2(j,2)>yy2;e2=e2+1;end          % index error count if above line
end
E=(e1+e2)/M;                            % net fractional error
title(['e1 = ',num2str(e1),', e2 = ',num2str(e2),' E = ',num2str(E)])
%
%
%%
% Solution to Exercise E14.1. M is the sample size per class (M_i in the text).
% Class separability is from the difference in means along the x axis and the covariance matrix.
%
M=100;                                              % sample size per class
m1x=43;m1y=50;m1z=50;                               % means for 3 variables of class 1
m2x=57;m2y=50;m2z=50;                               % means for 3 variables of class 2
m1=ones(M,3);m1(:,1)=m1(:,1)*m1x;                   % mean vectors for class 1
m1(:,2)=m1(:,2)*m1y;m1(:,3)=m1(:,3)*m1z;
m2=ones(M,3);m2(:,1)=m2(:,1)*m2x;                   % mean vectors for class 2
m2(:,2)=m2(:,2)*m2y;m2(:,3)=m2(:,3)*m2z;
s1x=10;s1y=10;s1z=7;r1xy= 0.4;r1xz= 0.1;r1yz=-0.1;  % cov matrix 1 parameters
s2x=20;s2y=10;s2z=10;r2xy=-0.6;r2xz=-0.1;r2yz= 0.8; % cov matrix 2 parameters
K1=[s1x^2 r1xy*s1x*s1y r1xz*s1x*s1z;...             % population cov matrix 1
    r1xy*s1x*s1y s1y^2 r1yz*s1y*s1z;...
    r1xz*s1x*s1z r1yz*s1y*s1z s1z^2];
% K2=[s2x^2 r2xy*s2x*s2y r2xz*s2x*s2z;...           % population cov matrix 2 (unused)
%     r2xy*s2x*s2y s2y^2 r2yz*s2y*s2z;...
%     r2xz*s2x*s2z r2yz*s2y*s2z s2z^2];
X1=mvnrnd(m1,K1,M);X2=mvnrnd(m2,K1,M);              % simulate the data X1, X2
X1b=mean(X1);K1x=cov(X1);X2b=mean(X2);K2x=cov(X2);  % sample statistics estimates
a=2*(X2b-X1b)/((K1x+K2x)/2);                        % lin discrim function parameter
b=X2b/((K1x+K2x)/2)*X2b'-X1b/((K1x+K2x)/2)*X1b';    % lin discrim function parameter
xp=[0 100 100 0];yp=[0 0 100 100];                  % x,y base plane for LDF
zp(1)=(b-a(2)*yp(1)-a(1)*xp(1))/a(3);               % 4 z-axis corners of LDF
zp(2)=(b-a(2)*yp(2)-a(1)*xp(2))/a(3);
zp(3)=(b-a(2)*yp(3)-a(1)*xp(3))/a(3);
zp(4)=(b-a(2)*yp(4)-a(1)*xp(4))/a(3);
plot3(X1(:,1),X1(:,2),X1(:,3),'k.'),axis square,hold on
plot3(X2(:,1),X2(:,2),X2(:,3),'r.'),plot3([m1x m2x],[m1y m2y],[m1z m2z],'b-o','linewidth',4)
patch(xp,yp,zp,'y')                                 % plot LDF as yellow plane
set(gca,'FontSize',18),xlabel('x'),ylabel('y'),zlabel('z')
xlim([0 100]),zlim([0 100]),zlim([0 100])
e1=0;e2=e1;                                         % initialize error counters
for j=1:M                                           % detect and count errors
    zz1=(b-a(2)*X1(j,2)-a(1)*X1(j,1))/a(3);         % find LDF for x,y values of X1
    zz2=(b-a(2)*X2(j,2)-a(1)*X2(j,1))/a(3);         % find LDF for x,y values of X2
    if X1(j,3)<zz1;e1=e1+1;end          % index error count if below line
    if X2(j,3)>zz2;e2=e2+1;end          % index error count if above line
end
E=(e1+e2)/(2*M);                        % net fractional classification error
title(['e1 = ',num2str(e1),', e2 = ',num2str(e2),' E = ',num2str(E)])
%
% Compute principal components via SVD
%
X=[X1;X2];Xb=mean(X);                   % combine data from both classes
[~,~,U]=svd(X-Xb);                      % apply SVD after subtracting overall sample mean
figure,subplot(2,2,2),plot3(X1(:,1),X1(:,2),X1(:,3),'k.'),axis square,hold on
plot3(X2(:,1),X2(:,2),X2(:,3),'r.')     % plot data again, then the 3 PCs
plot3([-50*U(1,1) 50*U(1,1)]+Xb(1),[-50*U(2,1) ...
50*U(2,1)]+Xb(2),[-50*U(3,1) 50*U(3,1)]+Xb(3),'k','LineWidth',4)
plot3([-50*U(1,2) 50*U(1,2)]+Xb(1),[-50*U(2,2) ...
50*U(2,2)]+Xb(2),[-50*U(3,2) 50*U(3,2)]+Xb(3),'b','LineWidth',4)
plot3([-50*U(1,3) 50*U(1,3)]+Xb(1),[-50*U(2,3) ...
50*U(2,3)]+Xb(2),[-50*U(3,3) 50*U(3,3)]+Xb(3),'m','LineWidth',4)
set(gca,'FontSize',18), xlabel('x'),ylabel('y'),zlabel('z')
xlim([0 100]),zlim([0 100]),zlim([0 100])
%
% Select largest 2 of 3 singular vectors, compute LDF and count classification errors
%
Y=X*U;Y1=Y(1:M,1:2);Y2=Y(M+1:2*M,1:2);  % principal components
xp=round(min([Y1(:,1);Y2(:,1)])):round(max([Y1(:,1);Y2(:,1)])); % x-axis for LDF
subplot(2,2,3),plot(Y1(:,1),Y1(:,2),'k.'),axis square,hold on
plot(Y2(:,1),Y2(:,2),'r.')
set(gca,'FontSize',18), xlabel('x'),ylabel('y')
Y1b=mean(Y1);K1y=cov(Y1);Y2b=mean(Y2);K2y=cov(Y2);  % estimate sample statistics
a=2*(Y2b-Y1b)/((K1y+K2y)/2);                        % lin discrim function parameter
b=Y2b/((K1y+K2y)/2)*Y2b'-Y1b/((K1y+K2y)/2)*Y1b';    % lin discrim function parameter
yp=(b-a(1)*xp)/a(2);                                % LDF
plot(xp,yp,'k','linewidth',2)
e1=0;e2=e1;                                         % initialize error counters
for j=1:M
    yy1=(b-a(1)*Y1(j,1))/a(2);          % find LDF for x,y values of X1
    yy2=(b-a(1)*Y2(j,1))/a(2);          % find LDF for x,y values of X2
    if Y1(j,2)<yy1;e1=e1+1;end          % index error count if below line
    if Y2(j,2)>yy2;e2=e2+1;end          % index error count if above line
end
E=(e1+e2)/(2*M);                        % net fractional classification error
title(['e1 = ',num2str(e1),', e2 = ',num2str(e2),' E = ',num2str(E)])
%
% Select largest 1 of 3 variables, compute LDFs and count classification errors
%
subplot(2,2,4),histogram(Y1(:,1)),hold on,histogram(Y2(:,1))
set(gca,'FontSize',18), xlabel('u_1'),ylabel('occurence')
a=2*(Y2b(1)-Y1b(1))/((K1y(1,1)+K2y(1,1))/2);                  % lin discrim function parameter
b=Y2b(1)^2/((K1y(1,1)+K2y(1,1))/2) - Y1b(1)^2/((K1y(1,1)+K2y(1,1))/2);
xpp=b/a;                                % LDF for 1D data
e1=0;e2=e1;                             % initialize error counters
for j=1:M
    if Y1(j,1)>xpp;e1=e1+1;end          % index error count if above threshold
    if Y2(j,1)<xpp;e2=e2+1;end          % index error count if below threshold
end
E=(e1+e2)/(2*M);                        % net fractional classification error
title(['threshold',num2str(xpp),' e1 = ',num2str(e1),' e2 = ',num2str(e2),' E = ',num2str(E)])
%
%
%%
% Solution to Example 16.19. Results are plotted in Figure 16.13.
% (a) program DFT matrices Qd (see Eq A.10) and Qq (see Eq 16.60) to compute a frequency spectra
%
n=5;N=2^n;q=exp(1i*2*pi/N);     % n=# qubits, N=2^n=# basis states; Fourier kernel
J=zeros(N);Qd=J;Qq=J;           % initialize
for j=1:N                       % form index matrix, 1st row & column are all zeros
    J(j+1,:)=0:j:j*(N-1);       % QFT matrix index
    for k=1:N
        Qd(j,k)=q^(-(j-1)*(k-1));   % classical DFT matrix (negative exponent)
        Qq(j,k)=q^(mod(J(j,k),N));  % quantum DFT matrix (positive exponent)
    end
end
% Form time and frequency axes
%
m=1:N;du=1/N;U=1;u=0:du:U-du;                       % m=time axis, u=frequency axis
g=(sin(2*pi*(m-1)/(N/2))+2*sin(2*pi*(m-1)/(N/4)));  % input = sum of two sine waves
qsum=sum(g.^2);g=g/sqrt(qsum);                      % normalize input energy
subplot(4,1,1),plot(m,g,'LineWidth',1.5),set(gca,'fontsize',18)
xlabel('$t_m$ (s)','interpreter','latex'),ylabel('$g(t_m)$','interpreter','latex')
G=Qd*g'/sqrt(N);                                    % DFT result via Qd
subplot(4,1,2),stem(u,abs(G),'LineWidth',1.5),xlim([0,0.5]),set(gca,'fontsize',18)
xlabel('$u_k$ (Hz)','interpreter','latex'),ylabel('$\tilde{g}[k]$','interpreter','latex')
gt=Qq*g'/sqrt(N);                                   % DFT result using Qq
subplot(4,1,3),stem(u,abs(gt),'LineWidth',1.5),xlim([0,0.5]),set(gca,'fontsize',18)
xlabel('$u_k$ (Hz)','interpreter','latex'),ylabel('$\tilde{g}[k]$','interpreter','latex')
% figure,subplot(2,1,1),imagesc(real(Qd)),subplot(2,1,2),imagesc(real(Qq)) % to image the matrices
%
% (b) Program the QFT circuit for n=5 qubits.
%
ket0=[1;0];ket1=[0;1];x=zeros(1,n);Gq=0;  % 1-qubit basis, x=decimal array of binary value numbers 
for k=N:2*N-1                   % to control binary digit length, begin with N and go to 2N-1
    xx=dec2bin(k);              % convert decimal index to binary
    for i=1:n                   % separate binary digits from xx into decimal array x
        if xx(i+1)=='1'         % e.g., k'=34 -> xx='100010' -> x=0,0,0,1,0 as a string...
            x(i)=1;             %...of decimal number that indicate k=34-32=2 
        else
            x(i)=0;
        end 
    end
% Form the n=5 qubit states. In Fig 16.12: X1=\tilde{x}_0, X2=\tilde{x}=1, etc.
    X1=ket0+ket1*exp(1i*2*pi*x(5)/2);                                           % LSB
    X2=ket0+ket1*exp(1i*2*pi*(x(4)/2+x(5)/2^2));
    X3=ket0+ket1*exp(1i*2*pi*(x(3)/2+x(4)/2^2+x(5)/2^3));
    X4=ket0+ket1*exp(1i*2*pi*(x(2)/2+x(3)/2^2+x(4)/2^3+x(5)/2^4));
    X5=ket0+ket1*exp(1i*2*pi*(x(1)/2+x(2)/2^2+x(3)/2^3+x(4)/2^4+x(5)/2^5));     % MSB
    Gq=Gq+g(k-N+1)*kron(X1,kron(X2,kron(X3,kron(X4,X5))))/2^(n/2); % weight & sum tensor products
end
subplot(4,1,4),stem(u,abs(Gq),'LineWidth',1.5),xlim([0,0.5]),set(gca,'fontsize',18)% (freq spect)^2
xlabel('$u_k$ (Hz)','interpreter','latex'),ylabel('$\tilde{g}_x$','interpreter','latex')
%