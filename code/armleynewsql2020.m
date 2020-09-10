clear all;
%
% Square-lake graphs of four Armley flood-mitigation scenarios
% with 4 options: 1) HW higher walls; 2) GRR river-bed widering and obstruction clearance;
%                 3) FP flood-pllain storage; 4) NFM (leaky dams and trees); 5) BD beaver dams.
% S1: FP=Calverley; S2: FP=Rodley; 3) F{P=Calverley+Rodley; S4: FP=Cononley Washlands+Holden Park
%
FEV = 9.34;  % Mcubic metres Armley Boxing Day 2015 flood
depth = 2;   % m
Lside = sqrt(FEV*10^6/depth);
xmax = Lside*1.5; % 950;
larr = 10
for nSc=1:4 
 switch nSc
 case 1 % S1: Calverley
  %  Vb = [0.08 0.1071]*FEV;    % FP M cubic metres 100 percent
  Vb = [0.07 0.11]*FEV;    % FP M cubic metres 100 percent
  Vgrr = [0.07 0.07]*FEV; % GRR M cubic metres 100 percent
  Vr = (1-(Vb(1)+Vgrr(1))/FEV)*[1 1]*FEV;     % HW M cubic metres 100 percent
  Vnfm = [0.01 0.05]*FEV;    % NFM M cubic metres 100 percent
  Vbeav = [0.0 0.01]*FEV;    % BD M cubic metres 100 percent
 case 2 % S2: Rodley
  Vgrr = [0.07 0.07]*FEV; % GRR M cubic metres 100 percent
  % Vb = [0.12 0.2355]*FEV;    % FP M cubic metres 100 percent
  Vb = [0.11 0.24]*FEV;    % FP M cubic metres 100 percent
  Vr = (1-(Vb(1)+Vgrr(1))/FEV)*[1 1]*FEV;     % HW M cubic metres 100 percent
  Vnfm = [0.01 0.05]*FEV;    % NFM M cubic metres 100 percent
  Vbeav = [0.0 0.01]*FEV;    % BD M cubic metres 100 percent
 case 3 % S3: Calverley+ Rodley
  Vgrr = [0.07 0.07]*FEV; % GRR M cubic metres 100 percent
  % Vb = [0.14 0.3426]*FEV;    % FP M cubic metres 100 percent
  Vb = [0.13 0.34]*FEV;    % FP M cubic metres 100 percent
  Vr = (1-(Vb(1)+Vgrr(1))/FEV)*[1 1]*FEV;
  Vnfm = [0.01 0.05]*FEV;    % NFM M cubic metres 100 percent
  Vbeav = [0.0 0.01]*FEV;    % BD M cubic metres 100 percent			      
 case 4 % S4: Cononley Washlands+Holden Park
  Vgrr = [0.07 0.07]*FEV; % GRR M cubic metres 100 percent
  % Vb = [0.504 0.79285]*FEV;  % FP M cubic metres 100 percent
  Vb = [0.45 0.45]*FEV;  % FP M cubic metres 100 percent
  Vr = (1-(Vb(1)+Vgrr(1))/FEV)*[1 1]*FEV;     % HW M cubic metres 100 percent
  Vnfm = [0.01 0.05]*FEV;    % NFM M cubic metres 100 percent
  Vbeav = [0.0 0.01]*FEV;    % BD cubic metres 100 percent
 end			       
 fac = 1.4;
 Lxv = [Vr, Vgrr, Vb, Vnfm, Vbeav]/FEV;
 res = mean(Lxv(1:2));  % HW
 grr = mean(Lxv(3:4));  % GRR
 tree = mean(Lxv(5:6)); % FPS
 nfm = mean(Lxv(7:8));  % NFM
 bea = mean(Lxv(9:10)); % BD

 lower_fev = [Lxv(1) Lxv(3) Lxv(5) Lxv(7) Lxv(9)];
 upper_fev = [Lxv(2) Lxv(4) Lxv(6) Lxv(8) Lxv(10)];
 lower_fev_sum = cumsum(lower_fev);
 upper_fev_sum = cumsum(upper_fev);
 alph = 0.4;
 fs = 16;
 %%
 figure(100+nSc); clf;
 axis('equal');
 Lx = Lside;
 Ly = Lside;
 % Box
 xr = Lx*[0,1,1,0,0];
 yr = Ly*[0,0,1,1,0];
 plot(xr,yr,'-k','linewidth',3); hold on;
 %
 xr = Lx*[0,res,res,0,0];
 yr = Ly*[0,0,1,1,0];
 plot(xr,yr,'-b','linewidth',2); hold on;
 patch(xr,yr,[0 0 0.99]); alpha(alph); hold on;
 %
 xn = Lx*[res,res+grr,res+grr,res,res];
 yn = Ly*[0,0,1,1,0];
 plot(xn,yn,'-k','linewidth',2); hold on;
 patch(xn,yn,[0.99 0 0]); alpha(alph); hold on;
 %xt = Lx*[res+nfm,res+nfm+tree,res+grr+tree, res+grr,res+grr];
 xt = [lower_fev_sum(2),lower_fev_sum(3),upper_fev_sum(3),upper_fev_sum(2),lower_fev_sum(2)]*Lside;
 yt = Ly*[0,0,1,1,0];
 plot(xt,yt,'-k','linewidth',2); hold on;
 patch(xt,yt,[0 0.99 0]); alpha(alph); hold on;
 %
 xt = [lower_fev_sum(3),lower_fev_sum(4),upper_fev_sum(4),upper_fev_sum(3),lower_fev_sum(3)]*Lside;
 yt = Ly*[0,0,1,1,0];
 plot(xt,yt,'-k','linewidth',2); hold on;
 patch(xt,yt,[0.99 0 0.99]); alpha(alph); hold on;
 %
 xt = [lower_fev_sum(4),lower_fev_sum(5),upper_fev_sum(5),upper_fev_sum(4),lower_fev_sum(4)]*Lside;
 yt = Ly*[0,0,1,1,0];
 plot(xt,yt,'-k','linewidth',2); hold on;
 patch(xt,yt,[0.4 0.6 0.7]); alpha(alph); hold on;
 %
 arr = annotation('doublearrow');
 arr.Parent = gca;
 arr.X = Lx*[0 1];
 arr.Y = [0.2*Ly 0.2*Ly];
 arr.Color = 'black'; arr.Head1Style = 'vback3'; arr.Head2Style = 'vback3';
 arr.Head1Length = larr; arr.Head2Length = larr; arr.LineWidth = 2; 
 %
 arr = annotation('doublearrow');
 arr.Parent = gca;
 arr.X = Lx*[0 Lxv(1)];
 arr.Y = [0.8*Ly 0.8*Ly];
 arr.Color = 'black'; arr.Head1Style = 'vback3'; arr.Head2Style = 'vback3';
 arr.Head1Length = larr; arr.Head2Length = larr; arr.LineWidth = 2; 
 %
 arn = annotation('doublearrow');
 arn.Parent = gca;
 arn.X = Lx*[Lxv(1) Lxv(1)+Lxv(3)];
 arn.Y = [0.6*Ly 0.6*Ly];
 arr.Color = 'black'; arr.Head1Style = 'vback3'; arr.Head2Style = 'vback3';
 arr.Head1Length = larr; arr.Head2Length = larr; arr.LineWidth = 2; 
 %
 art = annotation('doublearrow');
 art.Parent = gca;
 art.X = Lx*[res+grr res+grr+Lxv(5)];
 art.Y = [0.4*Ly 0.4*Ly];
 arr.Color = 'black'; arr.Head1Style = 'vback3'; arr.Head2Style = 'vback3';
 arr.Head1Length = larr; arr.Head2Length = larr; arr.LineWidth = 2; 
 %	 				     
 facc = (Vb(2)-Vb(1)+Vnfm(1)+Vnfm(2)+Vbeav(2))/FEV/2;
 plot(Lx*[res+grr+Lxv(5)+facc,res+grr+Lxv(5)+facc],Ly*[0 1],'--k','linewidth',2);
 arm = annotation('doublearrow');
 arm.Parent = gca;
 arm.X = Lx*[res+grr res+grr+Lxv(6)];
 arm.Y = [0.97*Ly 0.97*Ly];
 arr.Color = 'black'; arr.Head1Style = 'vback3'; arr.Head2Style = 'vback3';
 arr.Head1Length = larr; arr.Head2Length = larr; arr.LineWidth = 2; 
 %
 % arm = annotation('arrow');
 % arm.Parent = gca;
 % arm.X = Lx*[res+grr res+grr+Lxv(5)+0.25*Lxv(6)];
 % arm.Y = [0.5*Ly 0.8*Ly];
 % arr.Color = 'black'; arr.Head1Style = 'vback3'; arr.Head2Style = 'vback3';
 % arr.Head1Length = larr; arr.Head2Length = larr; arr.LineWidth = 2; 
 %
 xlabel('Sidelength (m)','fontsize',fs);
 ylabel('Sidelength (m)','fontsize',fs);
 axis([0 fac*Lside 0 Lside 0 depth]);
 %box on
 %%
 switch nSc
 case 1
   title(sprintf('S1: FEV $$\\approx %d^2$$m$$^2$$ x $$2$$m $$ \\approx %.2f$$Mm$$^3$$',...
   round(Lside,0),FEV),'Interpreter','latex','fontsize',18);
   gtext('HW 86$\%$ 2.10m','Interpreter','latex','fontsize',16);
   gtext('$\pounds$75.6M at $\pounds$0.88M$/\%$','Interpreter','latex','fontsize',16);
   gtext('GRR 7$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$10M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$1.43M$/\%$','Interpreter','latex','fontsize',16);
   gtext('FPS 7$\%$','Interpreter','latex','fontsize',16);
   gtext('$\pounds$10M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$1.43M$/\%$','Interpreter','latex','fontsize',16);
   gtext('FPS 12$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$10M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$0.83M$/\%$','Interpreter','latex','fontsize',16);
   gtext('Total 100$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$95.6M ','Interpreter','latex','fontsize',16);
   gtext('NFM $1$-$5\%$','Interpreter','latex','fontsize',16);
   gtext('beavers $0$-$1\%$','Interpreter','latex','fontsize',16);
   gtext('Extra uncertain mitigation','Interpreter','latex','fontsize',14);
   gtext('for climate-change uptake','Interpreter','latex','fontsize',14);
  case 2
   title(sprintf('S2: FEV $$\\approx %d^2$$m$$^2$$ x $$2$$m $$ \\approx %.2f$$Mm$$^3$$',...
   round(Lside,0),FEV),'Interpreter','latex','fontsize',18);
   gtext('HW 82$\%$ 2.00m','Interpreter','latex','fontsize',16);
   gtext('$\pounds$72.2M at $\pounds$0.88M$/\%$','Interpreter','latex','fontsize',16);
   gtext('GRR 7$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$10M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$1.43M$/\%$','Interpreter','latex','fontsize',16);
   gtext('FPS 11$\%$','Interpreter','latex','fontsize',16);
   gtext('$\pounds$14M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$1.27M$/\%$','Interpreter','latex','fontsize',16);
   gtext('FPS 24$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$14M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$0.58M$/\%$','Interpreter','latex','fontsize',16);
   gtext('Total 100$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$96.2M ','Interpreter','latex','fontsize',16);
   gtext('NFM $1$-$5\%$','Interpreter','latex','fontsize',16);
   gtext('beavers $0$-$1\%$','Interpreter','latex','fontsize',16);
   gtext('Extra uncertain mitigation','Interpreter','latex','fontsize',14);
   gtext('for climate-change uptake','Interpreter','latex','fontsize',14);
 case 3
   title(sprintf('S3: FEV $$\\approx %d^2$$m$$^2$$ x $$2$$m $$ \\approx %.2f$$Mm$$^3$$',...
   round(Lside,0),FEV),'Interpreter','latex','fontsize',18);
   gtext('HW 80$\%$ 1.95m','Interpreter','latex','fontsize',16);
   gtext('$\pounds$96.2M at $\pounds$0.88M$/\%$','Interpreter','latex','fontsize',16);
   gtext('GRR 7$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$10M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$1.43M$/\%$','Interpreter','latex','fontsize',16);
   gtext('FPS 13$\%$','Interpreter','latex','fontsize',16);
   gtext('$\pounds$24M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$1.85M$/\%$','Interpreter','latex','fontsize',16);
   gtext('FPS 34$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$24M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$0.71M$/\%$','Interpreter','latex','fontsize',16);
   gtext('Total 100$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$104.4M ','Interpreter','latex','fontsize',16);
   gtext('NFM $1$-$5\%$','Interpreter','latex','fontsize',16);
   gtext('beavers $0$-$1\%$','Interpreter','latex','fontsize',16);
   gtext('Extra uncertain mitigation','Interpreter','latex','fontsize',14);
   gtext('for climate-change uptake','Interpreter','latex','fontsize',14);
case 4
   title(sprintf('S4: FEV $$\\approx %d^2$$m$$^2$$ x $$2$$m $$ \\approx %.2f$$Mm$$^3$$',...
   round(Lside,0),FEV),'Interpreter','latex','fontsize',18);
   gtext('HW 55$\%$ 1.34m','Interpreter','latex','fontsize',16);
   gtext('$\pounds$48.4M at $\pounds$0.88M$/\%$','Interpreter','latex','fontsize',16);
   gtext('GRR 7$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$10M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$1.43M$/\%$','Interpreter','latex','fontsize',16);
   gtext('FPS 45$\%$','Interpreter','latex','fontsize',16);
   gtext('$\pounds$35M','Interpreter','latex','fontsize',16);
   gtext('$\pounds$0.78M$/\%$','Interpreter','latex','fontsize',16);
   %gtext('FPS 79.3$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   %gtext('$\pounds$35M','Interpreter','latex','fontsize',16);
   %gtext('$\pounds$0.44M$/\%$','Interpreter','latex','fontsize',16);
   gtext('Total 100$\%$','fontsize',20,'Interpreter','latex','fontsize',16);
   gtext('$\pounds$93.4M ','Interpreter','latex','fontsize',16);
   gtext('NFM $1$-$5\%$','Interpreter','latex','fontsize',16);
   gtext('beavers $0$-$1\%$','Interpreter','latex','fontsize',16);
   gtext('Extra uncertain mitigation','Interpreter','latex','fontsize',14);
   gtext('for climate-change uptake','Interpreter','latex','fontsize',14);
 end
end
%
%
%
set(101, 'PaperPositionMode', 'auto');
set(102, 'PaperPositionMode', 'auto');
set(103, 'PaperPositionMode', 'auto');
set(104, 'PaperPositionMode', 'auto');
% figure(101); print -djpeg /users/amtob/dropbox/S1armgrrv3.jpg
% figure(102); print -djpeg /users/amtob/dropbox/S2armgrrv3.jpg
% figure(103); print -djpeg /users/amtob/dropbox/S3armgrrv3.jpg
% figure(104); print -djpeg /users/amtob/dropbox/S4armgrrv3.jpg
