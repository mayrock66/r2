%**************************************************************************
%*A M file function to save and plot the output of nanomos ****************
%**************************************************************************
function saveoutput(Ec,Ne,Ie,Ne_sub,E_sub,Te_sub,converge)

global transport_model;
global dirname;
global Ng_step Nd_step;
global delta_x Lsd Lg_top Lg_bot t_top t_bot t_si dx dy;
global Vg1 Vg2 Vs Vd_temp Vd_initial Vg_step Ng_step Vd_step Nd_step;
global filename t_vall max_subband;
global plot_IV plot_Ec3d plot_Ne3d plot_Ecsub plot_Te plot_Nesub;
global plot_Ec_IV plot_Ne_IV;
global N_dos Trans E;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPECIFY XI, YI, Vg_bias, and Vd_bias for plotting use only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=round((2*Lsd+Lg_top)/dx)+1;
XI=linspace(-(Lsd+Lg_top/2)*1e9,(Lsd+Lg_top/2)*1e9,Nx);%in nanometers
Ny=round((t_top+t_bot+t_si)/dy)+1;
YI=linspace(-t_top*1e9,(t_si+t_bot)*1e9,Ny);
Vg_bias=linspace(Vg1,Vg1+Vg_step*Ng_step,Ng_step+1);
Vd_bias=linspace(Vd_temp,Vd_temp+Vd_step*Nd_step,Nd_step+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       POST PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MEc=reshape(Ec(:,Ng_step+1,Nd_step+1),Nx,Ny);
trMEc=MEc';
MNe=reshape(Ne(:,Ng_step+1,Nd_step+1),Nx,Ny);
trMNe=MNe';

mkdir(dirname);
cd(dirname);

%raw_data is always stored!!!!!
%--------------------------------------------------
save output_rawdata;
%Save convergence data in a file
%save DOS.dat N_dos -ascii;
save trans.dat Trans -ascii;
save E.dat E -ascii;

fid=fopen('convergence.dat','a');

count=fprintf(fid,'%s\n','****************************************************************');
count=fprintf(fid,'%s','*******You ran a Nanomos simulation on the '); 
count=fprintf(fid,'%s',date);count=fprintf(fid,'%s\n',' *********');
count=fprintf(fid,'%s\n','********************* Using Nanomos 2.5 ************************');
count=fprintf(fid,'%s\n','********** The input file you used  was the following **********');
count=fprintf(fid,'%s\n','****************************************************************');

cd ..
fid1=fopen(filename,'r');
    while 1
        tline = fgetl(fid1);
        if ~ischar(tline), break, end
        count=fprintf(fid,'%s\n',tline);  
    end
status=fclose(fid1);
cd(dirname);

count=fprintf(fid,'%s\n','****************************************************************');
count=fprintf(fid,'%s\n','********** Convergence data for simulated bias points **********');
count=fprintf(fid,'%s\n','****************************************************************');

for i=1:Ng_step+1
   for j=1:Nd_step+1
      
       count=fprintf(fid,'%s\n','');
       count=fprintf(fid,'**** Gate voltage = %e\n**** Drain voltage = %e\n',Vg_bias(i),Vd_bias(j));
       count=fprintf(fid,'%s\n','');
       count=fprintf(fid,'%e\n',converge(:,i,j));
       
   end    
end

count=fprintf(fid,'%s\n','****************************************************************');
count=fprintf(fid,'%s\n','********** The end of the convergence data file*****************');
count=fprintf(fid,'%s\n','****************************************************************');

status=fclose(fid);
%******************************************** IV plots **********************************************
if plot_IV==1
% ID-VD CHARACTERISTICS (A/m), SAVED TO "id_vd.dat"
% -------------------------------------------------
if (Ng_step>=1 & Nd_step>=1)

Ng=size(Ie,1);    
figure(1)
plot(Vd_bias,Ie(1,:),'o-');
grid on
hold on
for i=2:Ng
    plot(Vd_bias,Ie(i,:),'o-');
end
xlabel('V_{DS} [V]')
ylabel('I_{DS} [\muA/\mum]')
title('I_{DS} vs. V_{DS}');
print -depsc ID_VD.ps;

temmm=[Vd_bias',Ie(:,:)'];
save ID_VD.dat temmm -ascii;

end

% ID-VG CHARACTERISTICS (A/m), SAVED TO "id_vg.dat"
% -------------------------------------------------
if (Ng_step>=1 & Nd_step==0)

figure(1)
semilogy(Vg_bias',Ie(:,1),'o-');
grid on
xlabel('V_{GS} [V]')
ylabel('I_{DS} [\muA/\mum]')
title('I_{DS} vs. V_{GS}');
print -depsc ID_VG.ps;

temmm=[Vg_bias',Ie(:,1)];
save ID_VG.dat temmm -ascii;

end

% ID-VD CHARACTERISTICS (A/m), SAVED TO "id_vd.dat"
% -------------------------------------------------
if (Nd_step>=1 & Ng_step==0)

figure(1)
plot(Vd_bias,Ie(1,:),'o-');
grid on
xlabel('V_{DS} [V]')
ylabel('I_{DS} [\muA/\mum]')
title('I_{DS} vs. V_{DS}');
print -depsc ID_VD.ps;

temmm=[Vd_bias',Ie(1,:)'];
save ID_VD.dat temmm -ascii;
end

end %if plot_Iv==1
%***********************************************************************************************
if plot_Ne_IV==1
% SUBBAND CHARGE DENSITY (/cm^2) vs X for Diff. Vg
% --------------------------------------------------------------
if (Ng_step>=1 | (Ng_step==0 & Nd_step==0))

figure(2)
for iii=1:Ng_step+1
 Ne1=sum(squeeze(Ne_sub(:,:,:,iii,Nd_step+1)),3);
 Ne2(:,iii)=sum(Ne1,2)*1e-4;
 plot(XI',Ne2(:,iii),'r-');
 hold on
 grid on
end
if Ng_step>=1
title('2D electron density along the channel at different Vg');
elseif Ng_step==0
title('2D electron density along the channel');    
end    
xlabel('X [nm]');
ylabel('N2D [cm^{-2}]');
print -depsc N2D_X1.ps

%figure(11)
%Ne1=sum(squeeze(Ne_sub(:,:,:,1,Nd_step+1)),3);
%Ne2(:,1)=sum(Ne1,2)*1e-4;
%semilogy(XI',Ne2(:,1),'r-');
%hold on
%grid on
%if Ng_step>0
%  for iii=2:Ng_step+1
%    Ne1=sum(squeeze(Ne_sub(:,:,:,iii,Nd_step+1)),3);
%    Ne2(:,iii)=sum(Ne1,2)*1e-4;
%    semilogy(XI',Ne2(:,iii),'r-');
%  end
%end
%title('2D electron density along the channel @Diff. VG');
%xlabel('X [nm]');
%ylabel('N2D [cm^{-2}]');
%print -depsc ./output/LOG_N2D_X.ps

temmm=[XI',Ne2];
save N2D_X1.dat temmm -ascii;

end

% SUBBAND CHARGE DENSITY (/cm^2) vs X for Diff. Vd
% --------------------------------------------------------------
if Nd_step>=1

figure(3)
for iii=1:Nd_step+1
 Ne1=sum(squeeze(Ne_sub(:,:,:,Ng_step+1,iii)),3);
 Ne2(:,iii)=sum(Ne1,2)*1e-4;
 plot(XI',Ne2(:,iii),'r-');
 hold on
 grid on
end
title('2D electron density along the channel at different Vd');
xlabel('X [nm]');
ylabel('N2D [cm^{-2}]');
print -depsc N2D_X2.ps

temmm=[XI',Ne2];
save N2D_X2.dat temmm -ascii;

end

end %if plot_Ne_IV
%******************************************************************************************
if plot_Ec_IV==1
% The First SUBBAND ENERGY PROFILE vs X for Diff. Vg
% ------------------------------------------------------

if (Ng_step>=1 | (Ng_step==0 & Nd_step==0))
figure(4);
for iii=1:Ng_step+1
 plot(XI',E_sub(:,1,1,iii,Nd_step+1),'r-');
 hold on
 grid on
end

if Ng_step>=1 
title('The First Subband energy profile along the channel at different Vg');
elseif Ng_step==0
title('The First Subband energy profile along the channel');
end
xlabel('X [nm]');
ylabel('E_{SUB} [eV]');
print -depsc Ec_X1.ps

tem=E_sub(:,1,1,iii,Nd_step+1);
sq_tem=squeeze(tem);
temmm=[XI',sq_tem];
save Ec_X1.dat temmm -ascii;

end

% The First SUBBAND ENERGY PROFILE vs X for Diff. Vd
% ------------------------------------------------------

if Nd_step>=1
figure(5);
for iii=1:Nd_step+1
 plot(XI',E_sub(:,1,1,Ng_step+1,iii),'r-');
 hold on
 grid on
end

title('The First Subband energy profile along the channel at different Vd');
xlabel('X [nm]');
ylabel('E_{SUB} [eV]');
print -depsc Ec_X2.ps

tem=E_sub(:,1,1,Ng_step+1,:);
sq_tem=squeeze(tem);
temmm=[XI',sq_tem];
save Ec_X2.dat temmm -ascii;

end

end %if plot_Ec_IV
%*******************************************************************************************
if (plot_Ecsub==1 & max_subband>=1)
    figure(8)
        for iii=1:max_subband
            plot(XI',E_sub(:,iii,1,Ng_step+1,Nd_step+1),'r-');
            hold on
            grid on
            if (t_vall==3)
            plot(XI',E_sub(:,iii,2,Ng_step+1,Nd_step+1),'k-');
            plot(XI',E_sub(:,iii,3,Ng_step+1,Nd_step+1),'-');
            end
        end
        
   title('The Subbands energy profile along the channel');
   xlabel('X [nm]');
   ylabel('E_{SUB} [eV]');
   print -depsc Ec_sub_X.ps 
    
end
%*******************************************************************************************
%*******************************************************************************************
if (plot_Nesub==1 & max_subband>=1)
    figure(9)
    
            for iii=1:max_subband    
                semilogy(XI',Ne_sub(:,iii,1,Ng_step+1,Nd_step+1),'r-');
                hold on
                grid on
                if (t_vall==3)
                   semilogy(XI',Ne_sub(:,iii,2,Ng_step+1,Nd_step+1),'k-');
                   semilogy(XI',Ne_sub(:,iii,3,Ng_step+1,Nd_step+1),'-');
                end
             end
    title('2D electron density of the subbands along the channel ');
    xlabel('X [nm]');
    ylabel('N2D [cm^{-2}]');

    print -depsc Ne_sub_X.ps 
    
end
%*****************************************************************************************
if transport_model==6%ENERGY TRANSPORT MODEL *****************************************
%*****************************************************************************************
if (plot_Te==1)
    figure(10)
        plot(XI',Te_sub(:,1,1,Ng_step+1,Nd_step+1),'r-');  
        
        title('The carrier Temperature along the channel for the first sub-band');
   xlabel('X [nm]');
   ylabel('Te_{SUB} [eV]');
   grid on;
   print -depsc Te_sub_X.ps 
    
   save Te_sub_X.dat Te_sub -ascii; 
       
end 
end   
%***************************************************************************************
% Ec(X,Y) 
% -------------------------------------------
if plot_Ec3d==1
figure(6);
surf(XI,YI,trMEc);
%shading interp commented out to reduce size of the .ps file
title('3D Conduction band edge potential profile')
hx=xlabel('X [nm]');
hy=ylabel('Y [nm]');
hz=zlabel('Ec [eV]');
print -depsc Ec_X_Y.ps

XII=[0,XI];
tem1=[YI',trMEc];
tem2=[XII;tem1];
save Ec_X_Y.dat tem2 -ascii;

end
% 3D CHARGE DENSITY N(X,Y) 
% ------------------------------------------------------------
if plot_Ne3d==1
figure(7);
surf(XI,YI,trMNe);
%shading interp commented out to reduce size of the .ps file
title('3D Electron density profile')
hx=xlabel('X [nm]');
hy=ylabel('Y [nm]');
hz=zlabel('Ne [m^{-3}]');
print -depsc Ne_X_Y.ps

XII=[0,XI];
tem1=[YI',trMNe];
tem2=[XII;tem1];
save Ne_X_Y.dat tem2 -ascii;
end

if (Ng_step==0 & Nd_step==0)
%	PLOT the graph of Transmission coefficient versus Energy
%------------------------------------------------------------------
if transport_model==4

	figure(11)
	plot(Trans,E);
	xlabel('Transmission Coefficient');
	ylabel('Energy (eV)');
	print -depsc Trans.ps
end
%	PLOT the Density of states versus energy along the device
%------------------------------------------------------------------
if transport_model==4

	figure(12)
	pcolor(XI,E,sqrt(N_dos(:,1:length(E)))');
	shading interp;
	% Square root of the DOS to get better visualization
	xlabel('Distance along the device (nm)');
	ylabel('Energy (eV)');
	print -djpeg100 DOS.jpg
	%print -depsc2 DOS.ps
end
%	PLOT the Density of states versus energy along the device
%------------------------------------------------------------------
if transport_model==5

	figure(12)
	pcolor(XI,E,sqrt(abs(N_dos(:,1:length(E))')));
	shading interp;
	% Square root of the DOS to get better visualization
	xlabel('Distance along the device (nm)');
	ylabel('Energy (eV)');
	print -djpeg100 DOS.jpg
	%print -depsc2 DOS.ps
end


end	%end if Ng_step==0 and Nd_step==0

cd ..

%**************** The end of saveoutput function ***************************
