clear all
close all
clc

q=1.6e-19;
eps0=8.85e-12;
k=1.38e-23;
T=300;
eps=11;
ni=2e16;
Mn=0.14;
Mp=0.06;
Nd=1e17*1e6;
Na=1e17*1e6;
W=[300 300];
W=W.*1e-9;
NL=[Nd -Na];
Delta=1e-5;
Eg=1,1;
U=[0   0.2 0.3 0.4 0.5 0.55 0.58 0.6 0.62 0.64 0.66 0.68 0.69 0.7];
K=[100 100 100 100 100 100  100  100  120 130  150  170  200 300];
Sx=20;
X(1)=0;
for i=1:length(W)
    X=[X max(X)+W(i)/Sx:W(i)/Sx:max(X)+W(i)];
end
I=length(X);
for i=1:I-1
    Dx(i)=X(i+1)-X(i);
end
N(1)=NL(1);
for i=1:length(NL)
    N=[N ones(1,Sx).*NL(i)];
end
MN=ones(1,I).*Mn;
MP=ones(1,I).*Mp;
ft=k*T/q;
L0=sqrt(eps*eps0*ft/q/ni);
t0=1e-6;
m0=L0^2/t0/ft;
x=X./L0;
dx=Dx./L0;
N=N./ni;
mn=MN./m0;
mp=MP./m0;
u=U./ft;
fip(1)=log(N(1)/2+sqrt((N(1)/2)^2+1));
fip(I)=-log(-N(I)/2+sqrt((N(I)/2)^2+1));
for i=2:I-1
    fip(i)=fip(1)+(fip(I)-fip(1))/(x(I)-x(1))*(x(i)-x(1));
end
fip='fip';
nn='graf_';
for j=1:length(u)
    if j>1
    for i=1:I
        D(i)=(u(j)-u(j-1))*(1-(fi(I)-fi(i))/(fi(I)-fi(1)));
        fip(i)=fi(i)+D(i);
    end
    %plot(X.*1e9,fip.*ft)
    end
    if j==1
            fn=ones(I,1);
            fp=ones(I,1);
    else
        A=zeros(I,I);
        B=zeros(I,1);
        A(1,1)=1;
        A(I,I)=1;
        B(1)=1;
        B(I)=exp(-u(j));
        for i=2:I-1
            A(i,i)=2/(dx(i)+dx(i-1))*((mn(i)*exp(fip(i))*(-1)/dx(i))-(mn(i-1)*exp(fip(i-1))*(1/dx(i-1))));
            A(i,i-1)=2/(dx(i)+dx(i-1))*(mn(i-1)*exp(fip(i-1))/dx(i-1));
            A(i,i+1)=2/(dx(i)+dx(i-1))*(mn(i)*exp(fip(i))/dx(i));
        end
        fn=A^(-1)*B

        A=zeros(I,I);
        B=zeros(I,1);
        A(1,1)=1;
        A(I,I)=1;
        B(1)=1;
        B(I)=exp(u(j));
        for i=2:I-1
            A(i,i)=2/(dx(i)+dx(i-1))*((mn(i)*exp(-fip(i))*(-1)/dx(i))-(mn(i-1)*exp(-fip(i-1))*(1/dx(i-1))));
            A(i,i-1)=2/(dx(i)+dx(i-1))*(mn(i-1)*exp(-fip(i-1))/dx(i-1));
            A(i,i+1)=2/(dx(i)+dx(i-1))*(mn(i)*exp(-fip(i))/dx(i));
        end
        fp=A^(-1)*B
    end
    err=1;
    while err>Delta
        A=zeros(I,I);
        B=zeros(I,1);
        A(1,1)=1;
        A(I,I)=1;
        B(1)=log(N(1)/2+sqrt((N(1)/2)^2+1));
        B(I)=-log(-N(I)/2+sqrt((N(I)/2)^2+1))+u(j);

        for i=2:I-1
            A(i,i)=2/(dx(i)+dx(i-1))*(-1/dx(i)-1/dx(i-1));
            A(i,i-1)=2/(dx(i)+dx(i-1))/dx(i-1);
            A(i,i+1)=2/(dx(i)+dx(i-1))/dx(i);
            B(i)=fn(i)*exp(fip(i))-fp(i)*exp(-fip(i))-N(i);
        end

        fi=A^(-1)*B;
        clc
        err=max(abs(fi-fip))/max(abs(fip))
        V=U(j)
        pause(1e-6)
        fip=fip+(fi-fip)/K(j);
        subplot(1,3,1)
        plot(X.*1e9,-fi.*ft,'--',X.*1e9,-fi.*ft+Eg/2,X.*1e9,-fi.*ft-Eg/2,'LineWidth',1.5)
        xlabel('X, nm','FontSize',14)
        ylabel('Fi, eV','FontSize',14)
        xlim([min(X.*1e9) max(X.*1e9)])
        grid on
        subplot(1,3,2)
        plot(X.*1e9,fn.*exp(fi).*(ni*1e-6),'--',X.*1e9,fp.*exp(-fi).*(ni*1e-6),'LineWidth',1.5)
        xlabel('X, nm','FontSize',14)
        ylabel('n, p, cm^-^3','FontSize',14)
        xlim([min(X.*1e9) max(X.*1e9)])
        grid on
        subplot(1,3,3)
        semilogy(X.*1e9,fn.*exp(fi).*(ni*1e-6),'--',X.*1e9,fp.*exp(-fi).*(ni*1e-6),'LineWidth',1.5)
        xlabel('X, nm','FontSize',14)
        ylabel('n, p, cm^-^3','FontSize',14)
        xlim([min(X.*1e9) max(X.*1e9)])
        grid on
        set(gcf,'Position',[3 378 1427 420])
    end
    nnn=[nn num2str(j)];
    print(gcf,'-dmeta',nnn)
    jn(j)=mn(I)*exp(fi(I))*(fn(I)-fn(I-1))/dx(I-1);
    jp(j)=-mp(I)*exp(-fi(I))*(fp(I)-fp(I-1))/dx(I-1);
    J(j)=(jn(j)+jp(j))*q*ni*L0/t0;
    pause
end
figure
plot(U,-J.*1e-4,'LineWidth',1.5)
xlabel('U, V','FontSize',14)
ylabel('J, A/cm^2','FontSize',14)
xlim([min(U) max(U)])
grid on
print(gcf,'-dmeta','vah')
