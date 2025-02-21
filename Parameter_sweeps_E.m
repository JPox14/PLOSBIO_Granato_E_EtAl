close(figure(4))
clear;clc;
%cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here

%-------------Modify these values according to desired figure generation---

%Save the Figure? 1 = Yes; 0 = No
SaveFig = 0;
SaveData = 1;

%To generate Figure 1: Private_Attacker = 0 Private_Target = 0;
%To generate Figure 3: Private_Attacker = 1 Private_Target = 1;
%To generate Figure SY: Private_Attacker = 0 Private_Target = 1;  
%To generate Figure SZ: Private_Attacker = 1 Private_Target = 0;  


% Set 1 = True. 0 = False. If both = 1 Then Fig_f
conjugation = 1; %if 1 and toxins = 0. Then Fig_e
toxins =      1; %if 1 and conjugation = 0. Then Fig_d

% Set 1 = True. 0 = False
%Private Nutrients. Both 0 = Fig 1.
Private_Attacker = 0;
Private_Target = 0;

%------------Do not modify anything below this line------------------------


r = [1,1,1];       %Growth rate for 3 strains
KN = [5,5,5];        %Monod
if Private_Attacker + Private_Target == 1
    N = [1,0.25,0.25];      %Nutrients
else
    N = [1,1,1];      %Nutrients
end
E = 10;
Km = 1;            %Toxin affinity for target
b = 0;
c = 0;          %Toxin production rate - Key value and not sure what to put
A0 = 0.01;        %Initial donor abundance
T0 = 0.01;        %Initial recipient abundance
Tr0 = 0;          %Initial transconjugant abundance  
NO_D = 0;       %Donor Private Nutrient Access? 0 1 (False True)
NO_Tr = 0;        %Target Private Nutrient Access? 0 1      
tend = 500000;
k=0;


if Private_Attacker == 1
    NO_D = 1;       %Donor Private Nutrient Access? 0 1 (False True)
end
if Private_Target == 1
    NO_Tr = 1;        %Target Private Nutrient Access? 0 1
end
if conjugation == 1
    b = 0.25;   %conjugation rate
end
if toxins == 1
    c = 0.15;   %toxin production rate
end

min_val = 0;
max_val = 50;
steps = 1000;  
bindex = min_val:(max_val-min_val)/(steps-1):max_val;
don_v=zeros(1,steps-1);trans_v=zeros(1,steps-1);rec_v=zeros(1,steps-1);

for i = min_val:((max_val-min_val)/(steps-1)):max_val
    k = k + 1;
    E = i;
    y = [A0 Tr0 T0 0 N(1) N(2) N(3)];
    y0= y;
    eventfunc = @(t,y) HGT_ss_3(t, y, r, KN, Km, c, b, E,NO_D, NO_Tr);
    optionsode=odeset('Events',eventfunc,'NonNegative',1:7);
    [t,y,te,ye,ie] = ode45(@(t,y) HGT_func_3(t, y, r, KN, Km, c, b, E, NO_D, NO_Tr), [0 tend], y0,optionsode);
    m = [t,y];
    don_v(1,k) = m(end,2);
    trans_v(1,k) = m(end,3);
    rec_v(1,k) = m(end,4);
end

tot = don_v+trans_v+rec_v;
figure(4)
%plot(bindex, rec_v./tot,'color',[0, 0.6196, 0.4510],'LineWidth',1.66)
hold on
%plot(bindex, don_v./tot, 'color',[0.8275, 0.3765, 0.1529],'LineWidth',1.66)
plot(bindex, trans_v./tot, '--','color',[0, 0.6196, 0.4510],'LineWidth',3.33)
%title('Starting Frequency')
xlabel('Toxin Potency')
ylabel('Transconjugant Frequency')
ylim([0 1])
x0=500;
y0=500;
width=500;
height=200;
set(gcf,'position',[x0,y0,width,height])
%yticklabels({'0', '0.5', '1'})
yticks([0 0.25 0.5 0.75 1])

if SaveFig == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/_Updated_Figures/
    A = num2str(Private_Attacker);
    B = num2str(Private_Target);
    C = num2str(conjugation);
    D = num2str(toxins);
    % if N(2) < 0.5
  filnam = append('Parameters_E',A,B,C,D);
  saveas(gcf,filnam,'svg')
  cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
    % else
    %     filnam = append('Initial_Ratio',A,B,C,D,'_High-N');
    %     saveas(gcf,filnam,'svg')
    %     cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
    % end
end

mat = [bindex' (trans_v./tot)'];
if SaveData == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/RawData/
    FigNum = 'FigS9c';
    if Private_Attacker == 0 & Private_Target == 0
        FigLet = '_Left';
    elseif Private_Attacker == 1 & Private_Target == 1
        FigLet = '_Right';
    end
    filnam = append(FigNum,FigLet,'.csv');
    csvwrite(filnam,mat);
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/
end