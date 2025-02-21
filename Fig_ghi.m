
%Change commenting of b and c in for loop to generate figXghi
%figXi = keep both commented
close(figure(2))
clear;clc;
%cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here

%-------------Modify these values according to desired figure generation---

%Save the Figure? 1 = Yes; 0 = No
SaveFig = 0;
SaveData = 0;

%To generate Figure 1: Private_Attacker = 0 Private_Target = 0;
%To generate Figure 3: Private_Attacker = 1 Private_Target = 1;
%To generate Figure S5: Private_Attacker = 0 Private_Target = 1;  


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
A1 = 0.01;         %average initial abundance
Tr0 = 0;            %Transconjugant initial abundance
NO_D = 0;       %Donor Private Nutrient Access? 0 1 (False True)
NO_Tr = 0;        %Target Private Nutrient Access? 0 1
steps = 100;        
bindex = -steps+(steps*2)/(steps):(steps*2)/(steps):steps-(steps*2)/(steps);
tend = 500000;
k=0;
don_v=zeros(1,steps-1);trans_v=zeros(1,steps-1);rec_v=zeros(1,steps-1);

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

for i = -steps+(steps*2)/(steps):(steps*2)/(steps):steps-(steps*2)/(steps)
    A0 = A1 + A1*i/steps;
    T0 = A1 - A1*i/steps;
    k = k + 1;
    y = [A0 Tr0 T0 0 N(1) N(2) N(3)];
    y0= y;
    eventfunc = @(t,y) HGT_ss_3(t, y, r, KN, Km, c, b, E,NO_D, NO_Tr);
    optionsode=odeset('Events',eventfunc,'NonNegative',1:7);
    [t,y,te,ye,ie] = ode45(@(t,y) HGT_func_3(t, y, r, KN, Km, c, b, E, NO_D, NO_Tr), [0 tend], y0,optionsode);
    m = [t,y];
    don_v(1,k) = m(end,2);
    trans_v(1,k) = m(end,3);
    rec_v(1,k) = m(end,4);
    trackA(1,k) = A0;
    trackT(1,k) = T0;
end

tot = don_v+trans_v+rec_v;
matrix = [bindex' (rec_v./tot)' (don_v./tot)' (trans_v./tot)'];

figure(2)
plot(bindex, rec_v./tot,'color',[0, 0.6196, 0.4510],'LineWidth',1.66)
hold on
plot(bindex, don_v./tot, 'color',[0.8275, 0.3765, 0.1529],'LineWidth',1.66)
plot(bindex, trans_v./tot, '--','color',[0, 0.6196, 0.4510],'LineWidth',3.33)
%title('Starting Frequency')
xlabel('Initial Frequency (Target : Attacker)')
ylabel('Final Frequency')
ylim([0 1.0])
xticklabels({'99:1', '1:1', '1:99'})
xticks([-steps 0 steps])

if SaveFig == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/_Updated_Figures/
    A = num2str(Private_Attacker);
    B = num2str(Private_Target);
    C = num2str(conjugation);
    D = num2str(toxins);
    % if N(2) < 0.5
  filnam = append('Initial_Ratio',A,B,C,D);
  saveas(gcf,filnam,'svg')
  cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
    % else
    %     filnam = append('Initial_Ratio',A,B,C,D,'_High-N');
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
        FigLet = 'i';
    elseif conjugation == 0 & toxins == 1
        FigLet = 'g';
    elseif conjugation == 1 & toxins == 0
        FigLet = 'h';
    end
    filnam = append(FigNum,FigLet,'.csv');
    csvwrite(filnam,matrix)
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/
end