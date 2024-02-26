function xp=RP_eval_fct(t,R,inputfile)
xp = zeros(2, length(t));
% xp(1)=R(2);
% Read the JSON file into a string
str = fileread(inputfile);

% Decode the JSON string into a structure array
data = jsondecode(str);
% Extract the values from the structure array
P_v = data.P_v;
P_inf = data.P_inf;
S = data.S;
nu_l = data.nu_l;
rho_l = data.rho_l;
r_inf=data.r_inf;
for i = 1:length(t)
    xp(1, i) = R(2,i);
    % Calculate xp(2) for each t value corresponding to P_v
    P_v_value = P_v(i);  % Use P_v corresponding to t(i)       
    xp(2,i)=((P_v_value-P_inf-S/R(1))/rho_l-2*nu_l*R(2)/R(1)+(1-(R(1)/r_inf)^2)/2*R(2)*R(2)-log(r_inf/R(1))*R(2)*R(2))/(log(r_inf/R(1))*R(1));
% Rdd==-((log(R) - log(rbound) - R^2/(2*rbound^2) + 1/2)*Rd^2 - (2*nu*Rd)/(R*rho) + deltap/rho - (2*gamma)/(R*rho))/(R*(log(R) - log(rbound)))
% xp(2)=((P_v-P_inf-S*R(1)-2*nu_l*R(2)/R(1))/rho_l+(1-(R(1)/r_inf)^2)/2*R(2)*R(2)-log(r_inf/R(1))*R(2)*R(2))/(log(r_inf/R(1))*R(1));
end