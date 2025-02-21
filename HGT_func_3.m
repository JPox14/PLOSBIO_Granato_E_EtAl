%For a unique nutrient supplementation 

function Y = HGT_func_3(t, y, r, KN, Km, c, b, E, NO_D, NO_Tr)

Y = zeros(7,1);

%Note that the manuscript text has C2 = Targets and C3 = Transconjugants,
%which is opposite as it is encoded here. 

%Attacker
Y(1) = y(1)*((y(5)/(y(5)+KN(1)))+(y(6)/(y(6)+KN(2)))*NO_D)*r(1)*(1-c);
%Transconjugant
Y(2) = y(2)*((y(5)/(y(5)+KN(1)))+(y(7)/(y(7)+KN(3)))*NO_Tr)*r(2)*(1-c) + (b * (y(1)*((y(5)/(y(5)+KN(1)))+(y(6)/(y(6)+KN(2)))*NO_D) + y(2)*((y(5)/(y(5)+KN(1)))+(y(7)/(y(7)+KN(3)))*NO_Tr)) * y(3));
%Target
Y(3) = (y(3)*((y(5)/(y(5)+KN(1)))+(y(7)/(y(7)+KN(3)))*NO_Tr)*r(3)) - (b * (y(1)*((y(5)/(y(5)+KN(1)))+(y(6)/(y(6)+KN(2)))*NO_D) + y(2)*((y(5)/(y(5)+KN(1)))+(y(7)/(y(7)+KN(3)))*NO_Tr)) * y(3)) - E*(y(4)/(y(4)+Km))*y(3);
%Toxin 
Y(4) = c*y(1)*((y(5)/(y(5)+KN(1)))+(y(6)/(y(6)+KN(2)))*NO_D) + c*y(2)*((y(5)/(y(5)+KN(1)))+(y(7)/(y(7)+KN(3)))*NO_Tr) - (y(4)/(y(4)+Km))*y(3);
%Nutrient 1 (Shared)
Y(5) = -(y(5)/(y(5)+KN(1))) * (y(1)+y(2)+y(3));
%Nutrient 2 (Attacker Private)
Y(6) = -(y(6)/(y(6)+KN(2))) * (y(1)*NO_D);
%Nutrient 3 (Target/Transconjugant Private)
Y(7) = -(y(7)/(y(7)+KN(3))) * (y(2)*NO_Tr + y(3)*NO_Tr);
end
