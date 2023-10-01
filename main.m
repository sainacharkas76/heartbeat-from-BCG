close all
    clear all
    clc 

    subjectsmap = readtable("subjectsmap.xlsx");

    subject = "S26";
    level = "3"; 
 


    %% Initialization of the various notable points of the signal for the specified subject

    sessionid = gsesisonwithguser(subject);

    ECG_L1_SJ = str2num(string(subjectsmap{1, subject}));
    ECG_L1_EJ = str2num(string(subjectsmap{2, subject}));
    ECG_L2_SJ = str2num(string(subjectsmap{3, subject}));
    ECG_L2_EJ = str2num(string(subjectsmap{4, subject}));
    ECG_L3_SJ = str2num(string(subjectsmap{5, subject}));
    ECG_L3_EJ = str2num(string(subjectsmap{6, subject}));
    ECG_L4_SJ = str2num(string(subjectsmap{7, subject}));
    ECG_L4_EJ = str2num(string(subjectsmap{8, subject}));
    VIS_L1_SJ = str2num(string(subjectsmap{9, subject}));
    VIS_L1_EJ = str2num(string(subjectsmap{10, subject}));
    VIS_L2_SJ = str2num(string(subjectsmap{11, subject}));
    VIS_L2_EJ = str2num(string(subjectsmap{12, subject}));
    VIS_L3_SJ = str2num(string(subjectsmap{13, subject}));
    VIS_L3_EJ = str2num(string(subjectsmap{14, subject}));
    VIS_L4_SJ = str2num(string(subjectsmap{15, subject}));
    VIS_L4_EJ = str2num(string(subjectsmap{16, subject}));
    VIS_L1_P2 = str2num(string(subjectsmap{17, subject}));
    VIS_L1_P4 = str2num(string(subjectsmap{18, subject}));
    VIS_L2_P2 = str2num(string(subjectsmap{19, subject}));
    VIS_L2_P4 = str2num(string(subjectsmap{20, subject}));
    VIS_L3_P2 = str2num(string(subjectsmap{21, subject}));
    VIS_L3_P4 = str2num(string(subjectsmap{22, subject}));
    VIS_L4_P2 = str2num(string(subjectsmap{23, subject}));
    VIS_L4_P4 = str2num(string(subjectsmap{24, subject}));
   
    ECG_L_SJ = 0;
    ECG_L_EJ = 0;
    VIS_L_SJ = 0;
    VIS_L_EJ = 0;
    VIS_L_P2 = 0;
    VIS_L_P4 = 0;
    
    if (level == "1")
        ECG_L_SJ = ECG_L1_SJ;
        ECG_L_EJ = ECG_L1_EJ;
        VIS_L_SJ = VIS_L1_SJ;
        VIS_L_EJ = VIS_L1_EJ;
        VIS_L_P2 = VIS_L1_P2;
        VIS_L_P4 = VIS_L1_P4;
    elseif (level == "2")
        ECG_L_SJ = ECG_L2_SJ;
        ECG_L_EJ = ECG_L2_EJ;
        VIS_L_SJ = VIS_L2_SJ;
        VIS_L_EJ = VIS_L2_EJ;
        VIS_L_P2 = VIS_L2_P2;
        VIS_L_P4 = VIS_L2_P4;
    elseif (level == "3")
        ECG_L_SJ = ECG_L3_SJ;
        ECG_L_EJ = ECG_L3_EJ;
        VIS_L_SJ = VIS_L3_SJ;
        VIS_L_EJ = VIS_L3_EJ;
        VIS_L_P2 = VIS_L3_P2;
        VIS_L_P4 = VIS_L3_P4;
    elseif (level == "4")
        ECG_L_SJ = ECG_L4_SJ;
        ECG_L_EJ = ECG_L4_EJ;
        VIS_L_SJ = VIS_L4_SJ;
        VIS_L_EJ = VIS_L4_EJ;
        VIS_L_P2 = VIS_L4_P2;
        VIS_L_P4 = VIS_L4_P4;
    end

    %% Sampling frequencies of the various sensors

    % Oculus Quest VR
    vis_sf = 71.43; 
    
    % Movisense
    acc_sf = 64;
    ecg_sf = 1024;

    
    %% Signal recovery
    
    % signals obtained from Movisense

    ecg = readtable('grecords/'+subject+'/ecg.csv');
    ecg = table2array(ecg);
    acc = readtable('grecords/'+subject+'/acc.csv');
    acc = table2array(acc);

    % signals obtained from Oculus Quest

    vdata = [];
    if exist('grecords/'+subject+'/records-game_'+sessionid+'_'+level+'.csv', 'file')==2
        vdata = readtable('grecords/'+subject+'/records-game_'+sessionid+'_'+level+'.csv');
    else
        vdata = readtable('grecords/'+subject+'/vis.csv');
    end
    
    PosX = vdata.PosX; 
    PosY = vdata.PosY;
    PosZ = vdata.PosZ;
    
    PosAngX = vdata.PosAngX;
    PosAngY = vdata.PosAngY;
    PosAngZ = vdata.PosAngZ;
    
 
    %% Start
    
   
vdata.Timestamp = vdata.Timestamp - vdata.Timestamp(1,1);
vdata.Timestamp = (vdata.Timestamp).*(10^(-3));
t_nuovo = (0:1/vis_sf:vdata.Timestamp(end))';

VR_PosX = interp1(vdata.Timestamp, vdata.PosX, t_nuovo);
VR_PosY = interp1(vdata.Timestamp, vdata.PosY, t_nuovo);
VR_PosZ = interp1(vdata.Timestamp, vdata.PosZ, t_nuovo);   
VR_PosAngX = interp1(vdata.Timestamp, vdata.PosAngX, t_nuovo);
VR_PosAngY = interp1(vdata.Timestamp, vdata.PosAngY, t_nuovo);
VR_PosAngZ = interp1(vdata.Timestamp, vdata.PosAngZ, t_nuovo);

timi =t_nuovo(10*vis_sf:end-10*vis_sf);
tima = timi-timi(1,1);

VR_PosX = VR_PosX(10*vis_sf:end-10*vis_sf);
VR_PosY = VR_PosY(10*vis_sf:end-10*vis_sf);
VR_PosZ = VR_PosZ(10*vis_sf:end-10*vis_sf);
VR_PosAngX = VR_PosAngX(10*vis_sf:end-10*vis_sf);
VR_PosAngY = VR_PosAngY(10*vis_sf:end-10*vis_sf);
VR_PosAngZ = VR_PosAngZ(10*vis_sf:end-10*vis_sf);

% timeless =table(VR_PosX, VR_PosY, VR_PosZ, VR_PosAngX, VR_PosAngY, VR_PosAngZ);

dur_vr=length(VR_PosX)./vis_sf;
len_ecg=round(dur_vr.*ecg_sf);
ecg_new = ecg(10*ecg_sf:len_ecg);

%filtering BCG

VR_PosX_filtered = filtHP(VR_PosX,vis_sf);
VR_PosY_filtered = filtHP(VR_PosY,vis_sf);
VR_PosZ_filtered = filtHP(VR_PosZ,vis_sf);
VR_PosAngX_filtered = filtHP(VR_PosAngX,vis_sf);
VR_PosAngY_filtered = filtHP(VR_PosAngY,vis_sf);
VR_PosAngZ_filtered = filtHP(VR_PosAngZ,vis_sf);

% filtering ECG

 [amp pos] = pan_tompkin(ecg_new,ecg_sf,0);

 VR_new =table(tima,VR_PosX_filtered, VR_PosY_filtered, VR_PosZ_filtered, VR_PosAngX_filtered, VR_PosAngY_filtered, VR_PosAngZ_filtered);
 

%%
%Template matching

[HRX, pos_picchiX, amp_picchiX,TemplateX] = SCG_template_matching_corr(diff(diff(VR_PosX_filtered)),vis_sf,'X',[]);
[HRY, pos_picchiY, amp_picchiY,TemplateY] = SCG_template_matching_corr(diff(diff(VR_PosY_filtered)),vis_sf,'Y',[]);
[HRZ, pos_picchiZ, amp_picchiZ,TemplateZ] = SCG_template_matching_corr(diff(diff(VR_PosZ_filtered)),vis_sf,'Z',[]);
[HRRX, pos_picchiRX, amp_picchiRX,TemplateRX] = SCG_template_matching_corr(diff(diff(VR_PosAngX_filtered)),vis_sf,'RX',[]);
[HRRY, pos_picchiRY, amp_picchiRY,TemplateRY] = SCG_template_matching_corr(diff(diff(VR_PosAngY_filtered)),vis_sf,'RY',[]);
[HRRZ, pos_picchiRZ, amp_picchiRZ,TemplateRZ] = SCG_template_matching_corr(diff(diff(VR_PosAngZ_filtered)),vis_sf,'RZ',[]);


%%
% figure
%  subplot(2,2,1)
%  plot(t_nuovo(1:end-2), diff(vdata.PosAngX))
% 
% subplot(2,2,2)
% plot(t_nuovo(1:end-2), diff(vr_filtered(:,1)));
% 
% subplot(2,2,3)
% plot(t_nuovo(1:end-2), diff(vdata.PosAngY));
% 
% subplot(2,2,4)
% plot(t_nuovo(1:end-2), diff(vr_filtered(:,2)));

%%
% Kinetic energy
W=76;%in KG
R=(0.0682.*W+430.124)./(2.*pi);%R is radius of the head
H=2.*R;%H is the head hight
SW=0.079.*W;%W is the body mass & SW is the head mass
SM=SW./32.2;
Ix=0.2*(R.^2+H.^2);
Iy=0.4.*SM.*H.^2;
Iz=Ix;

Vx=diff(VR_PosX_filtered);
Vy=diff(VR_PosY_filtered);
Vz=diff(VR_PosZ_filtered);

Wx=diff(VR_PosAngX_filtered);
Wy=diff(VR_PosAngY_filtered);
Wz=diff(VR_PosAngZ_filtered);

K_lin=0.5.*W.*(Vx.^2+Vy.^2+Vz.^2);
K_rot=0.5.*(Ix.*Wx.^2+Iy.*Wy.^2 +Iz.*Wz.^2);

[HRK_lin, pos_picchiK_lin, amp_picchiK_lin,TemplateK_lin] = SCG_template_matching_corr(K_lin,vis_sf,'K_lin',[]);
[HRK_rot, pos_picchiK_rot, amp_picchiK_rot,TemplateK_rot] = SCG_template_matching_corr(K_rot,vis_sf,'K_rot',[]);








