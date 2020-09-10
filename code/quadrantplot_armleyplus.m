%%% Script <quadrantplot_armley.m> carries out the following:
%
% > loads and extracts relevant flow-data for Armley
% > calculates rating curve from empirical formula
% > computes FEV and duration etc.
% > plots data using subroutines:
%   <plot3panel.m> 
%   <plotFEVhT.m>
%   <plot_h_year.m>
%   <plot_ratingcurve.m>

% TK, August 2018 - adapted from OBs original <flowdatafloods.m>
% TK later updated to deal with errors; OB amended April 2019 to include
% GRR.

% Run all, or section-by-section.

clear;

%% DATA
% from OB's original <flowdatafloods.m>
% Flood data analysis: Q = C(h-A)^B; A in m, B dimensionless (?), C in m * m^(-B)
%
% Armley Aire River excel files: 22 to 32181
%
armleyexc = [22, 32181]; % rows in .csv file where data lies
nlenarmleyfi = armleyexc(2)-armleyexc(1)+1;
armleymeimaartdayspmo = [31,30,31,31,30,31,30,31,31,28,31];
days2may = 31+28+31+30;
daysfromay = 31+30+31+31+30+31+30.0-1.0;
day15min = 24*4;
no15minarmley = (sum(armleymeimaartdayspmo))*day15min;
fprintf('number 15min intervals May 2015 to March 2016 Armley/Calder flood gauge: %g %g\n',no15minarmley,nlenarmleyfi);
ndecarmley(1:2) = [(sum(armleymeimaartdayspmo(1:7))+24), (sum(armleymeimaartdayspmo(1:7))+29)]*day15min;


%% open csv file and extract

fid = fopen('Armley F1707 Flow 15min May 15 to Mar 16.csv'); % Q data
dataflow = textscan(fid, '%s %f %s %s %s %s', 'Delimiter', ',', 'HeaderLines', 21);
fclose(fid);

fid2 = fopen('Armley F1707 Stage 15min May 15 to Mar 16.csv'); % h data
datastage = textscan(fid2, '%s %f %s %f %f %s %f %s %s', 'Delimiter', ',', 'HeaderLines', 21);
fclose(fid2);


% extract Q data
timestampQ = datetime(dataflow{1},'InputFormat','dd/MM/yyyy HH:mm:SS');
nq = length(dataflow{2});
timeQ = zeros(1,nq);
dischargeQ = zeros(1,nq);
for ii = 1:nq
    timeQ(ii) = 0.25*ii/24; % Times 0.25 gives hours, 1/24 gives day
    dischargeQ(ii) = dataflow{2}(ii);
end

% extract h data
timestamph = datetime(datastage{1},'InputFormat','dd/MM/yyyy HH:mm:SS');
nh = length(datastage{2});
timeh = zeros(1,nh);
stage = zeros(1,nh);
for ii = 1:nh
    timeh(ii) = 0.25*ii/24;
    stage(ii) = datastage{2}(ii);
end

%% GRR ONNO 29-04-2019
%
Slope = 0.0005; % Need DTM OSmap data Tom? Now: between Forge and Armley 3 weirs of 2.5m each; 
                % but drop is ~10m (from 50m contour to 40m one over 5km), so 2.5m over 5000m
wgrr = 5;  % 5m
zgrr = 2; % corrected later below
Cm = 0.04; % Manning coefficient
bankslope = 2; % transverse bank slope
% (sqrt(Slope)/Cm)*((stage-zgrr)*wgrr).^(5/3)./(wgrr+(stage-zgrr)*sqrt(1+bankslope^2)).^(2/3)
%

%% Rating curve: from EA report for Armley gauge
arml = [0.2, 0.685, 1.917]; % lower stage limit
zgrr = arml(3);
armu = [0.685, 1.917, 4.17]; % upper stage limit
se = [0.0542, 0.0344, 0.0528]; % SE in m +-(5.42, 3.44, 5.28)% winter rating max deviation circa -20% and +16%
armll = (1.0-se).*arml;
armlu = (1.0+se).*arml;
armul = (1.0-se).*armu;
armuu = (1.0+se).*armu;

%rc coeffs
Crc = [30.69, 27.884, 30.127];
brc = [1.115, 1.462, 1.502];
arc = [0.156, 0.028, 0.153];

Ns = 10000;
[armlevc,ic] = max(stage);
rcstage = arml(1):(armlevc-arml(1))/Ns:armlevc; % linear h from hmin to hmax
[armlevc,ic] = max(rcstage);

rcdischarge = rcstage; % Q for rating curve
rcdischargeL = rcstage;
rcdischargeU = rcstage;

discharge = stage;
dischargeL = stage;
dischargeU = stage;

% GRR ONNO
rcdischargegrr = rcstage;
rcdischargeLgrr = rcstage;
rcdischargeUgrr = rcstage;


%%
% for ii = 1:nh % 
%     if (stage(ii) < armu(1)) && (stage(ii) > arml(1))
%         armleyflow2(ii) = Crc(1)*(stage(ii)-arc(1))^brc(1);
%     else
%         if (stage(ii) < armu(2))
%             armleyflow2(ii) = Crc(2)*(stage(ii)-arc(2))^brc(2);
%         else
%             %if (armleystage(ii) < armu(3))
%             armleyflow2(ii) = Crc(3)*(stage(ii)-arc(3))^brc(3);
%             %end
%         end
%     end
% end
% 
% %
% % For lower 
% %
% for ii = 1:nh
%     if (stage(ii) < armul(1)) && (stage(ii) > armll(1))
%         armleyflow2l(ii) = Crc(1)*(stage(ii)-arc(1))^brc(1);
%     else
%         if (stage(ii) < armul(2))
%             armleyflow2l(ii) = Crc(2)*(stage(ii)-arc(2))^brc(2);
%         else
%             %if (armleystage(ii) < armul(3))
%             armleyflow2l(ii) = Crc(3)*(stage(ii)-arc(3))^brc(3);
%             %end
%         end
%     end
% end
% 
% %
% % For upper 
% %
% for ii = 1:nh
%     
%     if (stage(ii) < armuu(1)) && (stage(ii) > armlu(1))
%         
%         armleyflow2u(ii) = Crc(1)*(stage(ii)-arc(1))^brc(1);
%         
%     elseif (stage(ii) < armuu(2)) && (stage(ii) > armlu(2))
%         
%             armleyflow2u(ii) = Crc(2)*(stage(ii)-arc(2))^brc(2);
%             
%     elseif (stage(ii) > armuu(2))
%         
%             armleyflow2u(ii) = Crc(3)*(stage(ii)-arc(3))^brc(3);
%             
%     end
% end


% Q = Q(h)
for ii = 1:Ns+1
    
    if (rcstage(ii) < armu(1)) && (rcstage(ii) > arml(1))
        
        rcdischarge(ii) = Crc(1)*(rcstage(ii)-arc(1))^brc(1);
        rcdischargeL(ii) = (1.0-se(1))*rcdischarge(ii); % -SE
        rcdischargeU(ii) = (1.0+se(1))*rcdischarge(ii); % +SE

        rcdischargegrr(ii) = Crc(1)*(rcstage(ii)-arc(1))^brc(1);
        
    elseif (rcstage(ii) < armu(2)) && (rcstage(ii) > arml(2))
        
        rcdischarge(ii) = Crc(2)*(rcstage(ii)-arc(2))^brc(2);
        rcdischargeL(ii) = (1.0-se(2))*rcdischarge(ii); % -SE
        rcdischargeU(ii) = (1.0+se(2))*rcdischarge(ii); % +SE

        rcdischargegrr(ii) = Crc(2)*(rcstage(ii)-arc(2))^brc(2);
        
        elseif (rcstage(ii) > armu(2)) % (sqrt(Slope)/Cm)*((stage-zgrr)*wgrr).^(5/3)./(wgrr+(stage-zgrr)*sqrt(1+bankslope^2)).^(2/3)
        
        rcdischarge(ii) = Crc(3)*(rcstage(ii)-arc(3))^brc(3);
        rcdischargeL(ii) = (1.0-se(3))*rcdischarge(ii); % -SE
        rcdischargeU(ii) = (1.0+se(3))*rcdischarge(ii); % +SE

        rcdischargegrr(ii) = Crc(3)*(rcstage(ii)-arc(3))^brc(3)+(sqrt(Slope)/Cm)*((rcstage(ii)-zgrr)*wgrr).^(5/3)./(wgrr+(rcstage(ii)-zgrr)*sqrt(1+bankslope^2)).^(2/3);
        
    end 
end

% lower and upper Q = Q(h(t)) = Q(t)
for ii = 1:nh
    
    if (stage(ii) < armu(1)) && (stage(ii) >= arml(1))
        
        discharge(ii) = Crc(1)*(stage(ii)-arc(1))^brc(1);
        dischargeL(ii) = (1.0-se(1))*discharge(ii); % -SE
        dischargeU(ii) = (1.0+se(1))*discharge(ii); % +SE

    elseif (stage(ii) < armu(2)) && (stage(ii) >= arml(2))
        
        discharge(ii) = Crc(2)*(stage(ii)-arc(2))^brc(2);
        dischargeL(ii) = (1.0-se(2))*discharge(ii); % -SE
        dischargeU(ii) = (1.0+se(2))*discharge(ii); % +SE

    elseif (stage(ii) >= armu(2))
        
        discharge(ii) = Crc(3)*(stage(ii)-arc(3))^brc(3);
        dischargeL(ii) = (1.0-se(3))*discharge(ii); % -SE
        dischargeU(ii) = (1.0+se(3))*discharge(ii); % +SE

        % dischargegrr(ii) = Crc(3)*(stage(ii)-arc(3))^brc(3)+(sqrt(Slope)/Cm)*((stage(ii)-zgrr)*wgrr).^(5/3)./(wgrr+(stage(ii)-zgrr)*sqrt(1+bankslope^2)).^(2/3);
    end
    
end

%% calculate FEV etc

armin = 2.7;
armax = 5.2;
Na = 20;
da = 0.1; % (armax-armin)/Na;
armstep = armin:da:armax;
Na = size(armstep,2);
nle = size(stage(ndecarmley(1):ndecarmley(2)),2);
arnle = zeros(1,nle);

for nna = 1:Na
    armcrit = armstep(nna); % m
    % [harmax,iarmax] = max(armleystage(ndecarmley(1):ndecarmley(2)));
    % [harmup,iarmup] = min(abs(armleystage(ndecarmley(1):ndecarmley(1)+iarmax)-armcrit));
    % [harmdo,iarmdo] = min(abs(armleystage(ndecarmley(1)+iarmax:ndecarmley(2))-armcrit));
    [harmax,iarmax] = max(stage(ndecarmley(1):ndecarmley(2)));
    [harmup,iarmup] = min(abs(stage(ndecarmley(1):ndecarmley(1)+iarmax-1)-armcrit));
    [harmdo,iarmdo] = min(abs(stage(ndecarmley(1)+iarmax-1:ndecarmley(2))-armcrit));
    iarmdo = iarmax+iarmdo;
    %
    %
    flowarmcrit = discharge(ndecarmley(1)+iarmup-1);
    armexcessvol = 15*60*sum( discharge(ndecarmley(1)+iarmup-1:ndecarmley(1)+iarmdo-1)-flowarmcrit ); % 15min
    flowarmcrit = Crc(3)*(armcrit-arc(3))^brc(3);
    posarm = max(discharge(ndecarmley(1):ndecarmley(2))-flowarmcrit,arnle);
    armexcessvol = 15*60*sum(posarm);
    armexcvol(nna) = armexcessvol;
end

%%
armcrit = 3.9; %ht
QT = Crc(3)*(armcrit-arc(3))^brc(3);
QTminus = QT*(1-0.055);
QTplus = QT*(1+0.055);

QTgrr = Crc(3)*(armcrit-arc(3))^brc(3)+(sqrt(Slope)/Cm)*((armcrit-zgrr)*wgrr).^(5/3)./(wgrr+(armcrit-zgrr)*sqrt(1+bankslope^2)).^(2/3);

posarm = max(discharge(ndecarmley(1):ndecarmley(2))-QT,arnle);
posarmL = max(dischargeL(ndecarmley(1):ndecarmley(2))-QTplus,arnle);
posarmU = max(dischargeU(ndecarmley(1):ndecarmley(2))-QTminus,arnle);
posarmgrr = max(discharge(ndecarmley(1):ndecarmley(2))-QTgrr,arnle);

% FEVs
armexcessvol = 15*60*sum(posarm);
armexcessvolL = 15*60*sum(posarmL); 
armexcessvolU = 15*60*sum(posarmU);
armexcessvolgrr = 15*60*sum(posarmgrr);

% flood durations
Tf = 15*60*nnz(posarm);
TfL = 15*60*nnz(posarmL);
TfU = 15*60*nnz(posarmU);

Qm = QT+armexcessvol/Tf;
hm = (Qm/Crc(3))^(1/brc(3))+arc(3);

%
%
[harmax,iarmax] = max(stage(ndecarmley(1):ndecarmley(2)));
[harmup,iarmup] = min(abs(stage(ndecarmley(1):ndecarmley(1)+iarmax)-armcrit));
[harmdo,iarmdo] = min(abs(stage(ndecarmley(1)+iarmax:ndecarmley(2))-armcrit));
iarmdo = iarmax+iarmdo;
%
%
flowarmcrit = discharge(ndecarmley(1)+iarmup-1);
flowarmcrit = Crc(3)*(armcrit-arc(3))^brc(3);
posarm = max(discharge(ndecarmley(1):ndecarmley(2))-flowarmcrit,arnle);
armexcessvol = 15*60*sum( discharge(ndecarmley(1)+iarmup-1:ndecarmley(1)+iarmdo-1)-flowarmcrit ); % 15min
armexcessvol22 = 15*60*sum(posarm);


%% data for plotting

t = timeh(ndecarmley(1):ndecarmley(2))-daysfromay; % time
h = stage(ndecarmley(1):ndecarmley(2)); % h data
% q = dischargeQ(ndecarmley(1):ndecarmley(2)); % q data
q = discharge(ndecarmley(1):ndecarmley(2)); % q data % Check ONNO 29-04-2019

q2 = discharge(ndecarmley(1):ndecarmley(2));
q2L = dischargeL(ndecarmley(1):ndecarmley(2));
q2U = dischargeU(ndecarmley(1):ndecarmley(2));

hrc = rcstage; % h rating curve
qrc = rcdischarge; % q rating curve
qrcgrr = rcdischargegrr; % q rating curve

% ONNO GRR

%% plotting routines

% 3 panel with h(t), Q(h), Q(t)
plot3panelgrr;

% ONNO GRR
% hrc_scaled = (hrc - h_min)./(h_max - h_min);
qrcgrr_scaled = (qrcgrr - q_min)./(q_max - q_min);
QTgrr_sc = (QTgrr - q_min)./(q_max - q_min);

figure(100);
plot(-hrc_scaled,qrcgrr_scaled,'--b','linewidth',1.5); hold on; % rating curve top left
plot([-hT_sc,1],[QTgrr_sc,QTgrr_sc],':b','linewidth',1);hold on;
text(0.8,1.08*QTgrr_sc,'$Q_{TGRR}$','fontsize',fs,'VerticalAlignment', 'middle','Interpreter','latex');

plot3panelerrgrr;
figure(1001);
plot(-hrc_scaled,qrcgrr_scaled,'--b','linewidth',1.5); hold on; % rating curve top left
plot([-hT_sc,1],[QTgrr_sc,QTgrr_sc],':b','linewidth',1);hold on;
text(0.8,1.08*QTgrr_sc,'$Q_{TGRR}$','fontsize',fs,'VerticalAlignment', 'middle','Interpreter','latex');

% NOTE on exporting figure from File > Export setup...:
% > Size: width = 25, height = 25, units = centimeters, expand axes to fill
% figure.
% > Rendering: resolution (dpi) = 300


% FEV and sidelength as a function of threshold
plotFEVhT;

% h(t) for whole year
plot_h_year;

% rating curve and discharge with errors
plot_ratingcurve;

