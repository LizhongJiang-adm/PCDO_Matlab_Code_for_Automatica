function [gmax,gmax_u] = PloyMax_Gra(tEnd,gEnd,dotgEnd,gEnd_u,dotgEnd_u)%a0,a1,a2,a3,tmax,
%UNTITLED3 此处提供此函数的摘要

%   此处提供详细说明
tLE = tEnd(:,1);            tRE = tEnd(:,2);
gLE = gEnd(:,1);            gRE = gEnd(:,2);
dotgLE = dotgEnd(:,1);      dotgRE = dotgEnd(:,2);

gLE_u = gEnd_u(:,1);          gRE_u = gEnd_u(:,2);
dotgLE_u = dotgEnd_u(:,1);    dotgRE_u = dotgEnd_u(:,2);


a3 = (2*(gLE - gRE) + (tRE - tLE).*(dotgLE + dotgRE)) / (tRE - tLE).^3;
a2 = (3*(gRE - gLE) - (tRE - tLE).*(2*dotgLE + dotgRE)) / (tRE - tLE).^2;
a1 = dotgLE;
a0 = gLE;


a3_u = -(2*gLE_u - 2*gRE_u - dotgLE_u*tLE - dotgRE_u*tLE + dotgLE_u*tRE + dotgRE_u*tRE)/(tLE - tRE)^3;
a2_u = -(3*gLE_u - 3*gRE_u - 2*dotgLE_u*tLE - dotgRE_u*tLE + 2*dotgLE_u*tRE + dotgRE_u*tRE)/(tLE - tRE)^2;
a1_u = dotgLE_u;
a0_u = gLE_u;



% tmax = tLE + (-2*a2 - sqrt(4*a2.^2 - 12*a3.*a1) )/6*a3;

if a3==0
    tmax = tLE;
end
if 4*a2.^2 <= 12*a3.*a1 & (a3~=0)
    tmax = tLE - a2./(3.*a3);
end
if (4*a2.^2 > 12*a3.*a1) & (a3~=0)
    tmax = tLE + (-2*a2 - sqrt(4*a2.^2 - 12*a3.*a1) )./(6*a3);
end


if tmax<=tLE

    gmax = a0;
    gmax_u = a0_u;

end
if tRE<=tmax
    gmax = gRE;
    gmax_u = gRE_u;
%     gmax = a3.*(tRE-tLE).^3 + a2.*(tRE-tLE).^2 + a1.*(tRE-tLE) + a0;
%     gmax_u = a3_u.*(tRE-tLE).^3 + a2_u.*(tRE-tLE).^2 + a1_u.*(tRE-tLE) + a0_u;

end
if tLE<tmax & tmax<=tRE

    if 4*a2.^2 <= 12*a3.*a1
        gmax = (2*a2^3 - 9*a1*a2*a3 + 27*a0*a3^2)/(27*a3^2);
        gmax_u = a0_u - (a2_u*(- 2*a2^2 + 3*a1*a3))/(9*a3^2) - (a2*a1_u)/(3*a3) + (a2*a3_u*(- 4*a2^2 + 9*a1*a3))/(27*a3^3);
    else
        gmax = (27*a0*a3^2 + 2*(a2^2 - 3*a1*a3)^(3/2) + 2*a2^3 - 9*a1*a2*a3)/(27*a3^2);
        gmax_u = a0_u - (a3_u*(4*(a2^2 - 3*a1*a3)^(3/2) + 4*a2^3 - 9*a1*a2*a3 + 9*a1*a3*(a2^2 - 3*a1*a3)^(1/2)))/(27*a3^3) + (a2_u*(2*a2*(a2^2 - 3*a1*a3)^(1/2) - 3*a1*a3 + 2*a2^2))/(9*a3^2) - (a1_u*(a2 + (a2^2 - 3*a1*a3)^(1/2)))/(3*a3);
    end

end

% gmax = gRE;
% gmax_u = gRE_u;
% gmax = dotgLE;
% gmax_u = dotgLE_u;
% gmax = a0;
% gmax_u = a0_u;

if ~isreal(gmax_u)
    error(num2str(tmax));
end

end