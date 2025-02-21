%HGT Model - JPalmer - April 14, 23 start. 
clear;clc;
close(figure(1))
%cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here

%-------------Modify these values according to desired figure generation---

%Save the Figure? 1 = Yes; 0 = No
SaveFig = 0;
SaveData = 0;

%To generate Figure 1: Private_Attacker = 0 Private_Target = 0;
%To generate Figure 3: Private_Attacker = 1 Private_Target = 1;
%To generate Figure S5: Private_Attacker = 0 Private_Target = 1;    

% Set 1 = True. 0 = False. If both = 1 Then Fig_c
conjugation =   0; %if 1 and toxins = 0. Then Fig_b
toxins =        1; %if 1 and conjugation = 0. Then Fig_a

% Set 1 = True. 0 = False
%Private Nutrients. Both 0 = Fig 1.
Private_Attacker = 0;
Private_Target = 0;

%------------Do not modify anything below this line------------------------

r = [1,1,1];      %Growth rate for 3 strains
KN = [5,5,5];       %Monod
if Private_Attacker + Private_Target == 1
    N = [1,0.25,0.25];      %Nutrients
else
    N = [1,1,1];      %Nutrients
end
E = 10;           %Toxin effectiveness
Km = 1;           %Toxin affinity for target
A0 = 0.01;        %Initial donor abundance
T0 = 0.01;        %Initial recipient abundance
Tr0 = 0;          %Initial transconjugant abundance  
NO_D = 0;
NO_Tr = 0;
b = 0;
c = 0;
%---------------Alter starting abundances----------------------------------
A = 0;          % percent change of current starting abundance
A = A/100;
A0 = A0 + A0*A;
T0 = T0 - T0*A;
%--------------------------------------------------------------------------

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


tend = 5000;
y = [A0 Tr0 T0 0 N(1) N(2) N(3)];
y0= y;

eventfunc = @(t,y) HGT_ss_3_Dynamics(t, y, r, KN, Km, c, b, E,NO_D, NO_Tr);
optionsode=odeset('Events',eventfunc,'NonNegative',1:7);
[t,y,te,ye,ie] = ode45(@(t,y) HGT_func_3(t, y, r, KN, Km, c, b, E, NO_D, NO_Tr), [0 tend], y0,optionsode);
M = [t, y];

Time = M(:,1);
Donor = M(:,2);
Transconjugant = M(:,3);
Target = M(:,4);
Toxin = M(:,5);
Nutrient1 = M(:,6);
Nutrient2 = M(:,7);
Nutrient3 = M(:,8);

Freq = Transconjugant(end)/(Donor(end) + Transconjugant(end) + Target(end));

figure(1)
plot(Time, Donor, 'color',[0.8353, 0.3686, 0],'LineWidth',1.66)
hold on
plot(Time,Transconjugant,'color',[0, 0.6196, 0.4510],'LineWidth',3.33,'LineStyle','--')
plot(Time,Target,'color',[0, 0.6196, 0.4510],'LineWidth',1.66)
%plot(Time,Toxin,'-','color',[0, 0, 0],'LineWidth',1)
%plot(Time,Nutrient1,'color','black','LineWidth',1,'linestyle',':')
if NO_D > 0 || NO_Tr > 0
    %plot(Time,Nutrient2,'color',[0.8353, 0.3686, 0],'LineWidth',1,'linestyle',':')
    %plot(Time,Nutrient3,'color',[0, 0.6196, 0.4510],'linewidth',1,'linestyle',':')
end
%title(num2str(Freq))
%title(num2str(Donor(end)))
%legend('Donor','Transconjugant','Target','Shared','Private_D','Private_T','Location','northeast')

if Private_Attacker + Private_Target < 1
    xlim([0 80])
    ylim([0 1.0])
elseif Private_Attacker + Private_Target == 1
    xlim([0 150])
    ylim([0 1.0])
else
    xlim([0 80])
    ylim([0 1.8])
end
if SaveFig == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/_Updated_Figures/
    A = num2str(Private_Attacker);
    B = num2str(Private_Target);
    C = num2str(conjugation);
    D = num2str(toxins);
    % if N(2) < 0.5
  filnam = append('Dynamics',A,B,C,D);
  saveas(gcf,filnam,'svg')
  cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
    % else
    %     filnam = append('Dynamics',A,B,C,D,'_High-N');
    %     saveas(gcf,filnam,'svg')
    %     cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
    % end
end

if SaveData == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/RawData/
    if Private_Attacker == 0 & Private_Target == 0
        FigNum = 'Fig1';
    elseif Private_Attacker == 1 & Private_Target == 1
        FigNum = 'Fig3';
    elseif Private_Attacker == 0 & Private_Target == 1
        FigNum = 'FigS5';
    end

    if conjugation == 1 & toxins == 1
        FigLet = 'c';
    elseif conjugation == 0 & toxins == 1
        FigLet = 'a';
    elseif conjugation == 1 & toxins == 0
        FigLet = 'b';
    end
    filnam = append(FigNum,FigLet,'.csv');
    csvwrite(filnam,M)
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/
end