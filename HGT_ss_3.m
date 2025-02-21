function [value,isterminal,direction] = HGT_ss_3(t, y, r, KN, Km, c, b, E, NO_D, NO_Tr)

dy = HGT_func_3(t, y, r, KN, Km, c, b, E,NO_D, NO_Tr);

SS = (abs(dy(5)) + abs(dy(6)) + abs(dy(7)));
val1 = SS - 0.000001;
value = val1;
isterminal = 1;
direction = -1;