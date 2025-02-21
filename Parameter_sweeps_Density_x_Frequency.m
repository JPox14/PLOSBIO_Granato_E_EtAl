close(figure(4))
close(figure(5))
close(figure(6))
clear;clc;
%cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here

%-------------Modify these values according to desired figure generation---

%Save the Figure? 1 = Yes; 0 = No
SaveFig = 0;
SaveData = 0;

%To generate Figure 1: Private_Attacker = 0 Private_Target = 0;
%To generate Figure 3: Private_Attacker = 1 Private_Target = 1;
%To generate Figure S5: Private_Attacker = 0 Private_Target = 1;   


% Set 1 = True. 0 = False. 
conjugation = 1; 
toxins =      1; 

% Set 1 = True. 0 = False
%Private Nutrients. Both 0 = Fig S8ace. Both 1 = FigS8bdf
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
E = 10;
Km = 1;            %Toxin affinity for target
b = 0;
c = 0;          %Toxin production rate - Key value and not sure what to put
ID = 0.01;        %Initial Density
Tr0 = 0;          %Initial transconjugant abundance  
NO_D = 0;       %Donor Private Nutrient Access? 0 1 (False True)
NO_Tr = 0;        %Target Private Nutrient Access? 0 1      
tend = 500000;
k=0;
count = 0;


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
steps = 100;  

minj_val = 0;
maxj_val = 1;
bjndex = minj_val:(maxj_val-minj_val)/(steps-1):maxj_val;

mini_val = 0.01;
maxi_val = 1;
bindex = mini_val:(maxi_val-mini_val)/(steps-1):maxi_val;
don_v=zeros(steps,steps-1);trans_v=zeros(steps,steps-1);rec_v=zeros(steps,steps-1);
trackA=zeros(steps,steps-1);trackT=zeros(steps,steps-1);

for i = mini_val:((maxi_val-mini_val)/(steps-1)):maxi_val
    count = count + 1;
    k = 0;
    for j = -steps+(steps*2)/(steps):(steps*2)/(steps):steps-(steps*2)/(steps)
        k = k + 1;
        %These vals are not being set correctly. 
        A0 = i/2 + i/2*(j/steps);
        T0 = i/2 - i/2*(j/steps);
        %
        y = [A0 Tr0 T0 0 N(1) N(2) N(3)];
        y0= y;
        eventfunc = @(t,y) HGT_ss_3(t, y, r, KN, Km, c, b, E,NO_D, NO_Tr);
        optionsode=odeset('Events',eventfunc,'NonNegative',1:7);
        [t,y,te,ye,ie] = ode45(@(t,y) HGT_func_3(t, y, r, KN, Km, c, b, E, NO_D, NO_Tr), [0 tend], y0,optionsode);
        m = [t,y];
        don_v(count,k) = m(end,2);
        trans_v(count,k) = m(end,3);
        rec_v(count,k) = m(end,4);
        %track(count,k) = A0/T0;
        trackA(count,k) = A0;
        trackT(count,k) = T0;
    end
end
trackA(end, end);
trackT(end,end);
track = A0./T0;
time = 1:1:count;
tot = don_v+trans_v+rec_v;
ab_mat = trans_v./tot;
cd_mat = don_v./tot;
ef_mat = rec_v./tot;

colorMap = zeros(101, 3); %Transconjugant green
colorMap(:,3) = linspace(1, 0.4510, 101);
colorMap(:,2) = linspace(1, 0.6196, 101);
colorMap(:,1) = linspace(1, 0, 101);

colorMap1 = zeros(101, 3); %Brown Attacker/Donor
colorMap1(:,3) = linspace(1, 0, 101);
colorMap1(:,2) = linspace(1, 0.3686, 101);
colorMap1(:,1) = linspace(1, 0.8353, 101);

colorMap2 = zeros(101, 3); %Black Target
colorMap2(:,3) = linspace(1, 0, 101);
colorMap2(:,2) = linspace(1, 0, 101);
colorMap2(:,1) = linspace(1, 0, 101);

figure(4)
imagesc(trans_v./tot)
set(gca,'YDir','normal')
colormap(colorMap)
caxis([0 1])
colorbar
%title('Transconjugant Frequency')
xlabel('Initial Ratio')
xticks([1 steps/2 steps-1])
yticklabels({mini_val maxi_val*.25 maxi_val*.5 maxi_val*.75 maxi_val})
xlabel('Initial Frequency (Target : Attacker)')
xticklabels({'99:1', '1:1', '1:99'})
ylabel('Initial Total Abundance')
yticks([1, steps*.25,steps/2,steps*.75, steps])
%xticks([0 Emax/2 Emax])

if SaveFig == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/_Updated_Figures/
    A = num2str(Private_Attacker);
    B = num2str(Private_Target);
    C = num2str(conjugation);
    D = num2str(toxins);
    % if N(2) < 0.5
  filnam = append('Parameters_dens_x_freq_Trans',A,B,C,D);
  saveas(gcf,filnam,'svg')
  cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
end

figure(5)
imagesc(don_v./tot)
set(gca,'YDir','normal')
colormap(colorMap1)
caxis([0 1])
colorbar
xlabel('Initial Ratio')
xticks([1 steps/2 steps-1])
yticklabels({mini_val maxi_val*.25 maxi_val*.5 maxi_val*.75 maxi_val})
xlabel('Initial Frequency (Target : Attacker)')
xticklabels({'99:1', '1:1', '1:99'})
ylabel('Initial Total Abundance')
yticks([1, steps*.25,steps/2,steps*.75, steps])

if SaveFig == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/_Updated_Figures/
    A = num2str(Private_Attacker);
    B = num2str(Private_Target);
    C = num2str(conjugation);
    D = num2str(toxins);
    % if N(2) < 0.5
  filnam = append('Parameters_dens_x_freq_Attack',A,B,C,D);
  saveas(gcf,filnam,'svg')
  cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
end

figure(6)
imagesc(rec_v./tot)
set(gca,'YDir','normal')
colormap(colorMap2)
caxis([0 1])
colorbar
xlabel('Initial Ratio')
xticks([1 steps/2 steps-1])
yticklabels({mini_val maxi_val*.25 maxi_val*.5 maxi_val*.75 maxi_val})
xlabel('Initial Frequency (Target : Attacker)')
xticklabels({'99:1', '1:1', '1:99'})
ylabel('Initial Total Abundance')
yticks([1, steps*.25,steps/2,steps*.75, steps])

if SaveFig == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/_Updated_Figures/
    A = num2str(Private_Attacker);
    B = num2str(Private_Target);
    C = num2str(conjugation);
    D = num2str(toxins);
    % if N(2) < 0.5
  filnam = append('Parameters_dens_x_freq_Target',A,B,C,D);
  saveas(gcf,filnam,'svg')
  cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
end




if SaveFig == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/_Updated_Figures/
    A = num2str(Private_Attacker);
    B = num2str(Private_Target);
    C = num2str(conjugation);
    D = num2str(toxins);
    % if N(2) < 0.5
  filnam = append('Parameters_dens_x_freq',A,B,C,D);
  saveas(gcf,filnam,'svg')
  cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/ %Put your directory here
end

if SaveData == 1
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/RawData/
    FigNum = 'FigS8';
    if Private_Attacker == 0 & Private_Target == 0
        FigLet_ab = 'a';
        FigLet_cd = 'c';
        FigLet_ef = 'e';
    elseif Private_Attacker == 1 & Private_Target == 1
        FigLet_ab = 'b';
        FigLet_cd = 'd';
        FigLet_ef = 'f';
    end

    filnam = append(FigNum,FigLet_ab,'.csv');
    csvwrite(filnam,ab_mat);
    filnam = append(FigNum,FigLet_cd,'.csv');
    csvwrite(filnam,cd_mat);
    filnam = append(FigNum,FigLet_ef,'.csv');
    csvwrite(filnam,ef_mat);
    cd /Users/jpalmer10/Documents/MATLAB/HGT/2024_Granato/Updated_code/
end