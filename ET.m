%********************************************************************************
%*********************************Solves for ET model****************************
%** Usage:***********************************************************************
% [net,Jet,Tet] =  ET(xvec,V,mu0,n0,nL,Tbc1,Tbc2,te0,s,ave,beta)*****************
%************ Written by Kausar Banoo / Purdue University************************
%********************************************************************************

function[nhd,Jhd,oldT]=ET(V,mu0,n0,nL,Tbc1,Tbc2,te0,ELECQ,vsat,beta)

global Lsd Lg_top Nx;

XI=linspace(-(Lsd+Lg_top/2)*1e9,(Lsd+Lg_top/2)*1e9,Nx);%in nanometers
xvec=XI;
% xvec in cm, V in V, mu0 in cm^2/V-s
% UNIT CONV.=========================
xvec=xvec*1e-7;		% nm --> cm
mu0=mu0*1e4;		% m^2/V-sec --> cm^2/V-sec
n0=n0*1e-6; nL=nL*1e-6;	% m^-3 --> cm^-3
% UNIT CONV.=========================

% Default: s = 0, ave = 3/2, Tbc1 = Tbc2 = 300
s=0; ave=1;
% beta = 2, te0 = 3e-13s

vt = vsat*1e2;	% cm/s
q = 1.6e-19;

fluxfac = (1/2)*n0*vt*q;

if(size(xvec,2)==1)
	n = size(xvec,1);
else
	n = size(xvec,2);
end

T0 = 300;
kb = 1.3807e-23; %J K^(-1)
Vth = kb*T0/q;
xi0 = vt/(2*mu0*Vth);
Ecrit = vt/mu0;
wfc = 5/2+s;

%For bulk pbc=1, for device pbc=0
%pbc = 1;
pbc = 0;
Tlower = 100;

maxiters = 40;

clear newEtvec;
clear Etvec;
clear newT;
clear oldT;
clear T;
clear Jhd;
clear nhd;
clear errhd;
clear nethd;
clear delT;

j=1;

clear Exvec;
clear midxvec;

for i=1:n-1,
    Exvec(i,1) = -(V(i+1) - V(i))/(xvec(i+1) - xvec(i))+1;
    edge = xvec(i);
    width = (xvec(i+1) - xvec(i));
    midxvec(i,1) = edge + width/2;
end
 
%First iteration
te = ones(n,1)*te0; %sec
Tbc1 = T0;
Tbc2 = T0; 
Etvec(:,j) = zeros(n-1,1);
T(:,j) = ones(n,1)*T0; %K


%First iteration with Et = 0
iters = 1;
A = sparse(zeros(2*n-2,2*n-2));
b = sparse(zeros(2*n-2,1));

i = 1;
L1 = (xvec(i+1) - xvec(i));
Ex1 = Exvec(i);

%mu1 = mu0*T0/T(i,j);

if(T(i,j) > T0)
%   mu1 = mu0/(1 + (3*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
    mu1 = mu0/(1 + (2*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
else
    mu1 = mu0;
end

A(2*i-1,2*i-1) = 1;
A(2*i-1,2*i) = 1/Ex1;

b(2*i-1) = q*mu1*n0;

for i=1:n-2,	

   L1 = (xvec(i+1) - xvec(i));
   Ex1 = Exvec(i);
   Ex2 = Exvec(i+1);

   A(2*i,2*i) = 1;
   A(2*i,2*i+2) = -1;
   A(2*i+1,2*i-1) = -exp(-Ex1*L1*q/(kb*T(i,j)));
   A(2*i+1,2*i) = -1/Ex1;
   A(2*i+1,2*i+1) = 1;
   A(2*i+1,2*i+2) = 1/Ex2;
      
end

i = n-1;
L1 = (xvec(i+1)-xvec(i));
Ex1 = Exvec(i);

A(2*i,2*i-1) = exp(-Ex1*L1*q/(kb*T(i,j)));
A(2*i,2*i) = 1/Ex1;

%mu1 = mu0*T0/T(i+1,j);     

if(T(i+1,j)>T0)
%     mu1 = mu0/(1 + (3*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
      mu1 = mu0/(1 + (2*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
else
      mu1 = mu0;
end

b(2*i) = q*mu1*nL;

solhd = A\b;

c1 = solhd(1);
c2 = solhd(2);

for i=2:n-1,
    c1 = [c1;solhd(2*i-1)];
    c2 = [c2;solhd(2*i)];
end

i = 1;
Ex = Exvec(i);

%mu = mu0*T0/T(i,j);

if(T(i,j) > T0)
%   mu = mu0/(1 + (3*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
    mu = mu0/(1 + (2*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
else
    mu = mu0;
end

Jhd(i,j) = c2(1);
nhd(i,j) = (c1(1) + c2(1)/Ex)/(q*mu);

for i=2:n,
    Ex = Exvec(i-1);
   
    %mu = mu0*T0/T(i,j);

    if(T(i,j) > T0)
%      mu = mu0/(1 + (3*kb*(T(i,j) - T0)*mu0/(2*q*te(i-1)*vt^2))^beta)^(1/beta);
       mu = mu0/(1 + (2*kb*(T(i,j) - T0)*mu0/(2*q*te(i-1)*vt^2))^beta)^(1/beta);
    else
       mu = mu0;
    end

    L = (xvec(i) - xvec(i-1));
    Jhd(i,j) = c2(i-1);
    nhd(i,j) = (exp(-Ex*L*q/(kb*T(i,j)))*c1(i-1) + c2(i-1)/Ex)/(q*mu);

end

%errhd(j,1) = max(Jhd(:,j)) - min(Jhd(:,j));

%Find the new temperature
if(pbc==1)
   Tbc1 = sum(T(:,j))/size(T(:,j),1);
   Tbc2 = Tbc1;
end

% Use the Band and Rose Newton method to calculate T
tmat = [0;0.5;1.0;2.0];
AEt = sparse(zeros(n-2,n-2));
bEt = sparse(zeros(n-2,1));
  
oldT(:,j) = T(:,j);

normdelT = 1000;

clear resiterT;
clear gofu;   
clear Vofk;
clear Vval;
iterT = 0;
t = 1;
tmax = 1;
while(normdelT > 1 & iterT < 20 & t>1e-3)
     
    iterT = iterT+1;
    for i=1:n-2,

        L1 = xvec(i+1) - xvec(i);
        L2 = xvec(i+2) - xvec(i+1);
  
        if(i>1)
           mu1 = 0;
           if(T(i+1,j) > T0)
%             mu1 = mu1 + 0.5*mu0/(1 + (3*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i+1)*vt^2))^beta)^(1/beta);
              mu1 = mu1 + 0.5*mu0/(1 + (2*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i+1)*vt^2))^beta)^(1/beta);
           else
              mu1 = mu1 + 0.5*mu0;
           end

           if(T(i,j) > T0)
%            mu1 = mu1 + 0.5*mu0/(1 + (3*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
             mu1 = mu1 + 0.5*mu0/(1 + (2*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
           else
             mu1 = mu1 + 0.5*mu0;
           end 
        end

        mu2 = 0;
        if(T(i+1,j) > T0)
%          mu2 = mu2 + 0.5*mu0/(1 + (3*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i+1)*vt^2))^beta)^(1/beta);
           mu2 = mu2 + 0.5*mu0/(1 + (2*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i+1)*vt^2))^beta)^(1/beta);
        else
           mu2 = mu2 + 0.5*mu0;
        end

        if(T(i+2,j) > T0)
%          mu2 = mu2 + 0.5*mu0/(1 + (3*kb*(T(i+2,j) - T0)*mu0/(2*q*te(i+2)*vt^2))^beta)^(1/beta);
           mu2 = mu2 + 0.5*mu0/(1 + (2*kb*(T(i+2,j) - T0)*mu0/(2*q*te(i+2)*vt^2))^beta)^(1/beta);
        else
           mu2 = mu2 + 0.5*mu0;
        end 
   
        n1 = 0.5*(nhd(i+1,j) + nhd(i,j));

        n2 = 0.5*(nhd(i+1,j) + nhd(i+2,j));


        coeffT32(i) = wfc*(kb^2/q)*ELECQ*mu2*n2/(2*L2);

        coeffT31(i) = wfc*(kb/q)*Jhd(i+1,j)/2;

        coeffT22(i) = -wfc*(kb^2/q)*ELECQ*(mu2*n2/(2*L2) + mu1*n1/(2*L1));

        coeffT21(i) = (ave+1)*(kb/q)*Jhd(i+1,j)*(1/2 - 1/2) - ave*(kb/te(i))*nhd(i+1,j)*(L1+L2)/2;

        coeffT12(i) = wfc*(kb^2/q)*ELECQ*mu1*n1/(2*L1);

        coeffT11(i) = -wfc*(kb/q)*Jhd(i+1,j)/2;

        const(i) = (Exvec(i+1)+Exvec(i))*Jhd(i+1,j)*(L1+L2)/4 + ave*(kb/te(i+1))*nhd(i+1,j)*T0*(L1+L2)/2;

        AEt(i,i) = 2*coeffT22(i)*oldT(i+1,j)+ coeffT21(i);
        if(i<n-2)
          AEt(i,i+1) = 2*coeffT32(i)*oldT(i+2,j) + coeffT31(i);
        end
        if(i>1)
           AEt(i,i-1) = 2*coeffT12(i)*oldT(i,j) + coeffT11(i);
        end

        if(i==1)
          bEt(i,1) = -(coeffT32(i)*oldT(i+2,j)^2 + coeffT31(i)*oldT(i+2,j) + coeffT22(i)*oldT(i+1,j)^2 + coeffT21(i)*oldT(i+1,j) + coeffT12(i)*Tbc1^2 + coeffT11(i)*Tbc1 + const(i));
        elseif(i==n-2)
          bEt(i,1) = -(coeffT32(i)*Tbc2^2 + coeffT31(i)*Tbc2 + coeffT22(i)*oldT(i+1,j)^2 + coeffT21(i)*oldT(i+1,j) + coeffT12(i)*oldT(i,j)^2 + coeffT11(i)*oldT(i,j) + const(i));
        else
          bEt(i,1) = -(coeffT32(i)*oldT(i+2,j)^2 + coeffT31(i)*oldT(i+2,j) + coeffT22(i)*oldT(i+1,j)^2 + coeffT21(i)*oldT(i+1,j) + coeffT12(i)*oldT(i,j)^2 + coeffT11(i)*oldT(i,j) + const(i));
        end

      end
  
      resiterT(iterT) = norm(bEt,inf);

      delT = AEt\bEt;
      normdelT = norm(delT,inf);

     %Bank and Rose Newton method

      for k=1:4,
        newT(1,j) = Tbc1;
        t = tmat(k);
        for i=1:n-2,
           newT(i+1,j) = oldT(i+1,j) + t*delT(i);
        end
        newT(n,j) = Tbc2;

       for i=1:n-2,
           if(i==1)
              gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*Tbc1^2 + coeffT11(i)*Tbc1 + const(i));
           elseif(i==n-2)
              gofu(i,1) = -(coeffT32(i)*Tbc2^2 + coeffT31(i)*Tbc2 + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
           else
              gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
              end
       end
       Voft(k) = delT'*gofu;

       if(k>1 & ((Voft(k)>0 & Voft(k-1)<0) | (Voft(k)<0 & Voft(k-1)>0)))
         if(tmat(k) >1)
	    t = 1;
            break;
         else
            % Use bisection
            Vval = Voft(k);
            V1 = Voft(k);
            V2 = Voft(k-1);
            t1 = tmat(k);
            t2 = tmat(k-1);
            bisect = 0;
            t = 1;
            while(abs(Vval) > 1e-3 & bisect < 20 & t>1e-3)
                bisect = bisect+1;
                t = (t1+t2)/2;

                newT(1,j) = Tbc1;
                for i=1:n-2,
                   newT(i+1,j) = oldT(i+1,j) + t*delT(i);
                end
                newT(n,j) = Tbc2;

                for i=1:n-2,
                   if(i==1)
                      gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*Tbc1^2 + coeffT11(i)*Tbc1 + const(i));
                   elseif(i==n-2)
                      gofu(i,1) = -(coeffT32(i)*Tbc2^2 + coeffT31(i)*Tbc2 + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
                   else
                     gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
                   end
                end
                Vval = delT'*gofu;
            
                if((Vval>0 & V1>0) | (Vval<0 & V1<0))
                   t1 = t;
                   V1 = Vval;
                elseif((Vval>0 & V2>0) | (Vval<0 & V2<0))
                   t2 = t;
                   V2 = Vval;
                end
               end  
            end
          elseif(k==4)
              tmax = tmax/2;
              t = tmax;
              break;
          end
      end

      newT(1,j) = Tbc1;
      for i=1:n-2,
         newT(i+1,j) = oldT(i+1,j) + t*delT(i);
      end
      newT(n,j) = Tbc2;
      

      for i=1:n-2,
           if(i==1)
              gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*Tbc1^2 + coeffT11(i)*Tbc1 + const(i));
           elseif(i==n-2)
              gofu(i,1) = -(coeffT32(i)*Tbc2^2 + coeffT31(i)*Tbc2 + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
           else
              gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
              end
      end

      if(norm(gofu,inf) <= resiterT(iterT))
         tmax = 1;
         oldT(:,j) = newT(:,j);
      end
         
      %figure(5);
      %plot(resiterT);
      %figure(6);
      %plot(xvec,oldT(:,j),xvec,newT(:,j));
      %pause;

end
 
for i=1:n-1, 
   L = xvec(i+1) - xvec(i);
   newEtvec(i,j) = (kb/q)*(oldT(i+1,j) - oldT(i,j))/L;
end


%figure(1);
%plot(xvec,Jhd(:,j)./nhd(:,j)/(-q));
%title(['v vs x for iters = ' num2str(iters) ' and V = ' num2str(Va(j))]);
%figure(2);
%plot(xvec,nhd(:,j));
%title(['n vs x for iters = ' num2str(iters) ' and V = ' num2str(Va(j))]);
%figure(3);
%plot(xvec,oldT(:,j));
%title(['T vs x for iters = ' num2str(iters) ' and V = ' num2str(Va(j))]);
%figure(4);
%plot(midxvec,newEtvec(:,j));
%title(['Et vs x for iters = ' num2str(iters) ' and V = ' num2str(Va(j))]);
%pause;
   

convg = 0;

while(iters<maxiters & convg==0)
  
    iters  = iters+1;
    T(:,j) = oldT(:,j);
    Etvec(:,j) = newEtvec(:,j);
    A = sparse(zeros(2*n-2,2*n-2));
    b = sparse(zeros(2*n-2,1));

    i = 1;
    Ex1 = Exvec(i);
    Et1 = Etvec(i,j);

    %mu1 = mu0*T0/T(i,j);
    if(T(i,j)>T0)
%       mu1 = mu0/(1 + (3*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
        mu1 = mu0/(1 + (2*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
     else
        mu1 = mu0;
     end

     %Two cases if Et is smaller or larger than Ex

     if(abs(Et1) > abs(Ex1/20))
        A(2*i-1,2*i-1) = ((kb*T(1,j)/q)/(kb*T(i,j)/q))^(-(Ex1+Et1)/Et1);
     else
         A(2*i-1,2*i-1) = 1;
     end 
 
     A(2*i-1,2*i) = 1/(Ex1+Et1);

     b(2*i-1) = q*mu1*n0;
   
     for i=1:n-2,	

         L1 = (xvec(i+1) - xvec(i));
         Ex1 = Exvec(i);
         Et1 = Etvec(i,j);
         Ex2 = Exvec(i+1);
         Et2 = Etvec(i+1,j);

         A(2*i,2*i) = 1;
         A(2*i,2*i+2) = -1;
         
	 if(abs(Et1) > abs(Ex1/20))
            A(2*i+1,2*i-1) = -((kb*T(i,j)/q + Et1*L1)/(kb*T(i,j)/q))^(-(Ex1+Et1)/Et1);
	 else
            A(2*i+1,2*i-1) = -exp(-(Ex1+Et1)*q*L1/(kb*T(i,j)));
         end             

	 A(2*i+1,2*i) = -1/(Ex1+Et1);

         if(abs(Et2) > abs(Ex2/20))
            A(2*i+1,2*i+1) = ((kb*T(i+1,j)/q)/(kb*T(i+1,j)/q))^(-(Ex2+Et2)/Et2);
         else
            A(2*i+1,2*i+1) = 1;
         end

         A(2*i+1,2*i+2) = 1/(Ex2+Et2);
   end

   i = n-1;
   L1 = (xvec(i+1)-xvec(i));
   Ex1 = Exvec(i);
   Et1 = Etvec(i,j);

   if(abs(Et1) > abs(Ex1/20))
      A(2*i,2*i-1) =  ((kb*T(i,j)/q + Et1*L1)/(kb*T(i,j)/q))^(-(Ex1+Et1)/Et1);
   else
      A(2*i,2*i-1) = exp(-(Ex1+Et1)*q*L1/(kb*T(i,j)));
   end

   A(2*i,2*i) = 1/(Ex1+Et1);

   %mu1 = mu0*T0/T(i+1,j);
   if(T(i+1,j)>T0)
%     mu1 = mu0/(1 + (3*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
      mu1 = mu0/(1 + (2*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
   else
      mu1 = mu0;
   end

   b(2*i) = q*mu1*nL;

   solhd = A\b;

   c1 = solhd(1);
   c2 = solhd(2);

   for i=2:n-1,
      c1 = [c1;solhd(2*i-1)];
      c2 = [c2;solhd(2*i)];
   end


   for i=1:n-1,
       L1 = (xvec(i+1)-xvec(i));
       Ex1 = Exvec(i);
       Et1 = Etvec(i,j);     
 
       %mu1 = mu0*T0/T(i,j);

       if(T(i,j)>T0)
%         mu1 = mu0/(1 + (3*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
          mu1 = mu0/(1 + (2*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
       else
          mu1 = mu0;
       end

       Jhd(i,j) = c2(i);

       if(abs(Et1)>abs(Ex1/20))
          nhd(i,j) = (c1(i)*((kb*T(i,j)/q)/(kb*T(i,j)/q))^(-(Ex1+Et1)/Et1) + c2(i)/(Ex1+Et1))/(q*mu1);
       else
          nhd(i,j) = (c1(i) + c2(i)/(Ex1+Et1))/(q*mu1);
       end
   end

  i=n-1;
  L1 = (xvec(i+1)-xvec(i));
  Ex1 = Exvec(i);
  Et1 = Etvec(i,j);  

  %mu1 = mu0*T0/T(i+1,j);
  if(T(i+1,j)>T0)
%    mu1 = mu0/(1 + (3*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
     mu1 = mu0/(1 + (2*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
  else
     mu1 = mu0;
  end

  Jhd(i+1,j) = c2(i);

  if(abs(Et1)>abs(Ex1/20))   
     nhd(i+1,j) = (c1(i)*((kb*T(i,j)/q + Et1*L1)/(kb*T(i,j)/q))^(-(Ex1+Et1)/Et1) + c2(i)/(Ex1+Et1))/(q*mu1);
  else
     nhd(i+1,j) = (c1(i)*exp(-(Ex1+Et1)*q*L1/(kb*T(i,j))) + c2(i)/(Ex1+Et1))/(q*mu1);
  end

   %Find the new temperature
   if(pbc==1)
     Tbc1 = sum(T(:,j))/size(T(:,j),1);
     Tbc2 = Tbc1;
   end

   tmat = [0;0.5;1.0;2.0];
   AEt = sparse(zeros(n-2,n-2));
   bEt = sparse(zeros(n-2,1));
  
   %oldT(:,j) = T(:,j);

   normdelT = 1000;

   clear resiterT;   
   clear gofu;
   clear Vofk;
   clear Vval;
   iterT = 0;
   t = 1;
   tmax = 1;
   while(normdelT > 1 & iterT < 20 & t>1e-3)
       
      iterT = iterT+1;
      for i=1:n-2,

        L1 = xvec(i+1) - xvec(i);
        L2 = xvec(i+2) - xvec(i+1);
  
        if(i>1)
           mu1 = 0;
           if(T(i+1,j) > T0)
%             mu1 = mu1 + 0.5*mu0/(1 + (3*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i+1)*vt^2))^beta)^(1/beta);
              mu1 = mu1 + 0.5*mu0/(1 + (2*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i+1)*vt^2))^beta)^(1/beta);
           else
              mu1 = mu1 + 0.5*mu0;
           end

           if(T(i,j) > T0)
%            mu1 = mu1 + 0.5*mu0/(1 + (3*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
             mu1 = mu1 + 0.5*mu0/(1 + (2*kb*(T(i,j) - T0)*mu0/(2*q*te(i)*vt^2))^beta)^(1/beta);
           else
             mu1 = mu1 + 0.5*mu0;
           end 
        end

        mu2 = 0;
        if(T(i+1,j) > T0)
%          mu2 = mu2 + 0.5*mu0/(1 + (3*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i+1)*vt^2))^beta)^(1/beta);
           mu2 = mu2 + 0.5*mu0/(1 + (2*kb*(T(i+1,j) - T0)*mu0/(2*q*te(i+1)*vt^2))^beta)^(1/beta);
        else
           mu2 = mu2 + 0.5*mu0;
        end

        if(T(i+2,j) > T0)
%          mu2 = mu2 + 0.5*mu0/(1 + (3*kb*(T(i+2,j) - T0)*mu0/(2*q*te(i+2)*vt^2))^beta)^(1/beta);
           mu2 = mu2 + 0.5*mu0/(1 + (2*kb*(T(i+2,j) - T0)*mu0/(2*q*te(i+2)*vt^2))^beta)^(1/beta);
        else
           mu2 = mu2 + 0.5*mu0;
        end 
   
        n1 = 0.5*(nhd(i+1,j) + nhd(i,j));

        n2 = 0.5*(nhd(i+1,j) + nhd(i+2,j));


        coeffT32(i) = wfc*(kb^2/q)*ELECQ*mu2*n2/(2*L2);

        coeffT31(i) = wfc*(kb/q)*Jhd(i+1,j)/2;

        coeffT22(i) = -wfc*(kb^2/q)*ELECQ*(mu2*n2/(2*L2) + mu1*n1/(2*L1));

        coeffT21(i) = (ave+1)*(kb/q)*Jhd(i+1,j)*(1/2 - 1/2) - ave*(kb/te(i))*nhd(i+1,j)*(L1+L2)/2;

        coeffT12(i) = wfc*(kb^2/q)*ELECQ*mu1*n1/(2*L1);

        coeffT11(i) = -wfc*(kb/q)*Jhd(i+1,j)/2;

        const(i) = (Exvec(i+1)+Exvec(i))*Jhd(i+1,j)*(L1+L2)/4 + ave*(kb/te(i+1))*nhd(i+1,j)*T0*(L1+L2)/2;

        AEt(i,i) = 2*coeffT22(i)*oldT(i+1,j)+ coeffT21(i);
        if(i<n-2)
          AEt(i,i+1) = 2*coeffT32(i)*oldT(i+2,j) + coeffT31(i);
        end
        if(i>1)
           AEt(i,i-1) = 2*coeffT12(i)*oldT(i,j) + coeffT11(i);
        end

        if(i==1)
          bEt(i,1) = -(coeffT32(i)*oldT(i+2,j)^2 + coeffT31(i)*oldT(i+2,j) + coeffT22(i)*oldT(i+1,j)^2 + coeffT21(i)*oldT(i+1,j) + coeffT12(i)*Tbc1^2 + coeffT11(i)*Tbc1 + const(i));
        elseif(i==n-2)
          bEt(i,1) = -(coeffT32(i)*Tbc2^2 + coeffT31(i)*Tbc2 + coeffT22(i)*oldT(i+1,j)^2 + coeffT21(i)*oldT(i+1,j) + coeffT12(i)*oldT(i,j)^2 + coeffT11(i)*oldT(i,j) + const(i));
        else
          bEt(i,1) = -(coeffT32(i)*oldT(i+2,j)^2 + coeffT31(i)*oldT(i+2,j) + coeffT22(i)*oldT(i+1,j)^2 + coeffT21(i)*oldT(i+1,j) + coeffT12(i)*oldT(i,j)^2 + coeffT11(i)*oldT(i,j) + const(i));
        end

      end
  
      resiterT(iterT) = norm(bEt,inf);
      
      delT = AEt\bEt;
      normdelT = norm(delT,inf);


     %Bank and Rose Newton method

      for k=1:4,
        newT(1,j) = Tbc1;
        t = tmat(k);
        for i=1:n-2,
           newT(i+1,j) = oldT(i+1,j) + t*delT(i);
        end
        newT(n,j) = Tbc2;

       for i=1:n-2,
           if(i==1)
              gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*Tbc1^2 + coeffT11(i)*Tbc1 + const(i));
           elseif(i==n-2)
              gofu(i,1) = -(coeffT32(i)*Tbc2^2 + coeffT31(i)*Tbc2 + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
           else
              gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
              end
       end
       Voft(k) = delT'*gofu;

       if(k>1 & ((Voft(k)>0 & Voft(k-1)<0) | (Voft(k)<0 & Voft(k-1)>0)))
         if(tmat(k) >1)
	    t = 1;
            break;
         else
            % Use bisection
            Vval = Voft(k);
            V1 = Voft(k);
            V2 = Voft(k-1);
            t1 = tmat(k);
            t2 = tmat(k-1);
            bisect = 0;
            t = 1;
            while(abs(Vval) > 1e-3 & bisect < 20 & t>1e-3)
                bisect = bisect+1;
                t = (t1+t2)/2;

                newT(1,j) = Tbc1;
                for i=1:n-2,
                   newT(i+1,j) = oldT(i+1,j) + t*delT(i);
                end
                newT(n,j) = Tbc2;

                for i=1:n-2,
                   if(i==1)
                      gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*Tbc1^2 + coeffT11(i)*Tbc1 + const(i));
                   elseif(i==n-2)
                      gofu(i,1) = -(coeffT32(i)*Tbc2^2 + coeffT31(i)*Tbc2 + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
                   else
                     gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
                   end
                end
                Vval = delT'*gofu;
            
                if((Vval>0 & V1>0) | (Vval<0 & V1<0))
                   t1 = t;
                   V1 = Vval;
                elseif((Vval>0 & V2>0) | (Vval<0 & V2<0))
                   t2 = t;
                   V2 = Vval;
                end
               end  
            end
          elseif(k==4)
              tmax = tmax/2;
              t = tmax;
              break;
          end
      end

      newT(1,j) = Tbc1;
      for i=1:n-2,
         newT(i+1,j) = oldT(i+1,j) + t*delT(i);
         %if(newT(i+1,j) < Tlower)
         %   newT(i+1,j) = Tlower;
         %end
      end
      newT(n,j) = Tbc2;

      for i=1:n-2,
           if(i==1)
              gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*Tbc1^2 + coeffT11(i)*Tbc1 + const(i));
           elseif(i==n-2)
              gofu(i,1) = -(coeffT32(i)*Tbc2^2 + coeffT31(i)*Tbc2 + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
           else
              gofu(i,1) = -(coeffT32(i)*newT(i+2,j)^2 + coeffT31(i)*newT(i+2,j) + coeffT22(i)*newT(i+1,j)^2 + coeffT21(i)*newT(i+1,j) + coeffT12(i)*newT(i,j)^2 + coeffT11(i)*newT(i,j) + const(i));
              end
      end
  
      if(norm(gofu,inf) <= resiterT(iterT))
         tmax = 1;
         oldT(:,j) = newT(:,j);
      end
         
      %figure(5);
      %plot(resiterT);
      %figure(6);
      %plot(xvec,newT(:,j));
      %pause;

  end
 
  for i=1:n-1, 
     L = xvec(i+1) - xvec(i);
     newEtvec(i,j) = (kb/q)*(oldT(i+1,j) - oldT(i,j))/L;
  end


%figure(1);
%plot(xvec,Jhd(:,j)./nhd(:,j)/(-q));
%title(['v vs x for iters = ' num2str(iters) ' and V = ' num2str(Va(j))]);
%figure(2);
%plot(xvec,nhd(:,j));
%title(['n vs x for iters = ' num2str(iters) ' and V = ' num2str(Va(j))]);
%figure(3);
%plot(xvec,oldT(:,j));
%title(['T vs x for iters = ' num2str(iters) ' and V = ' num2str(Va(j))]);
%figure(4);
%plot(midxvec,newEtvec(:,j));
%title(['Et vs x for iters = ' num2str(iters) ' and V = ' num2str(Va(j))]);
%pause;

  if(norm(oldT(:,j)-T(:,j),inf)<1)
    convg = 1; 
  end
   
end

%======UNIT CONV.=====================
nhd=nhd*1e6;		% cm^-3 --> m^-3
Jhd=-Jhd*1e4;		% #/cm^2-sec --> A/m^2
