
close(figure(3));
clear;clc;
%cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here

%-------------Modify these values according to desired figure generation---

%Save the Figure? 1 = Yes; 0 = No
SaveFig = 0;
SaveData = 1;

%To generate Figure 1: Private_Attacker = 0 Private_Target = 0;
%To generate Figure 3: Private_Attacker = 1 Private_Target = 1;
%To generate Figure S5: Private_Attacker = 0 Private_Target = 1;  


% Set 1 = True. 0 = False. If both = 1 Then Fig_f
conjugation = 1; %if 1 and toxins = 0. Then Fig_e
toxins = 1; %if 1 and conjugation = 0. Then Fig_d

% Set 1 = True. 0 = False
%Private Nutrients. Both 0 = Fig 1.
Private_Attacker = 1;
Private_Target = 1;

%------------Do not modify anything below this line------------------------

r = [1,1,1];       %Growth rate for 3 strains
KN = [5,5,5];        %Monod
if Private_Attacker + Private_Target == 1
    N = [1,0.25,0.25];      %Nutrients
else
    N = [1,1,1];      %Nutrients
end
Km = 1;            %Toxin affinity for target
A0 = 0.01;       %Initial attacker abundance
T0 = 0.01;       %Initial target abundance
Tr0 = 0;         %Initial tranconjugant abundance
NO_D = 0;       %attacker Private Nutrient Access? 0 1 (False True)
NO_Tr = 0;        %Target Private Nutrient Access? 0 1
steps = 249;
Emax = 50;          %Max toxin potency
bmax = 1.0;        %Max transfer rate
b1=(bmax/steps);
E1=(Emax/steps);
bindex = 0:b1:bmax; 
eindex = 0:E1:Emax;
c = 0;

if Private_Attacker == 1
    NO_D = 1;       %Donor Private Nutrient Access? 0 1 (False True)
end
if Private_Target == 1
    NO_Tr = 1;        %Target Private Nutrient Access? 0 1
end


tend = 500000;
y = [A0 0 T0 0 N(1) N(2) N(3)];
y0= y;k=0;
don_v=zeros(1,steps+1);trans_v=zeros(1,steps+1);rec_v=zeros(1,steps+1);
DON_V=zeros(1,steps+1);TRANS_V=zeros(1,steps+1);REC_V=zeros(1,steps+1);

for b = 0:b1:bmax
    %comment out both b = 0 and c = 0 to generate Fig 1f.
    if conjugation == 0
        b = 0;   %conjugation rate
    end
    if toxins == 1
        c = 0.15;   %toxin production rate
    end
    %space (Fig _e)
    k = k + 1;
    j = 0;
    for E = 0:E1:Emax
        j = j + 1;
        eventfunc = @(t,y) HGT_ss_3(t, y, r, KN, Km, c, b, E,NO_D, NO_Tr);
        optionsode=odeset('Events',eventfunc,'NonNegative',1:7);
        [t,y,te,ye,ie] = ode45(@(t,y) HGT_func_3(t, y, r, KN, Km, c, b, E, NO_D, NO_Tr), [0 tend], y0,optionsode);
        m = [t,y];
        don_v(k,j) = m(end,2);
        trans_v(k,j) = m(end,3);
        rec_v(k,j) = m(end,4);
    end
end

time = 1:1:k;
tot = don_v+trans_v+rec_v;
matrix=trans_v./tot;

colorMap = zeros(101, 3);
colorMap(:,3) = linspace(1, 0.4510, 101);
colorMap(:,2) = linspace(1, 0.6196, 101);
colorMap(:,1) = linspace(1, 0, 101);


figure(3)
imagesc(trans_v./tot)
set(gca,'YDir','normal')
colormap(colorMap)
caxis([0 1])
colorbar
%title('Transconjugant Frequency')
xlabel('Toxin Potency')
xticklabels({'0', Emax*0.2, Emax*0.4, Emax*0.6,Emax*0.8,Emax})
xticks([1 (steps+1)*0.2 (steps+1)*0.4 (steps+1)*0.6 (steps+1)*0.8 (steps+1)])
yticklabels({0 bmax*0.2 bmax*0.4 bmax*0.6 bmax*0.8 bmax})
yticks([1 (steps+1)*0.2 (steps+1)*0.4 (steps+1)*0.6 (steps+1)*0.8 steps+1])
ylabel('HGT Rate')
%yticks([0, bmax/2, bmax])
%xticks([0 Emax/2 Emax])


if SaveFig == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/_Updated_Figures/
    A = num2str(Private_Attacker);
    B = num2str(Private_Target);
    C = num2str(conjugation);
    D = num2str(toxins);
    % if N(2) < 0.5
  filnam = append('Heat_plot',A,B,C,D);
  saveas(gcf,filnam,'svg')
  cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
    % else
    %     filnam = append('Heat_plot',A,B,C,D,'_High-N');
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
        FigLet = 'f';
    elseif conjugation == 0 & toxins == 1
        FigLet = 'd';
    elseif conjugation == 1 & toxins == 0
        FigLet = 'e';
    end
    filnam = append(FigNum,FigLet,'.csv');
    csvwrite(filnam,matrix)
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/
end