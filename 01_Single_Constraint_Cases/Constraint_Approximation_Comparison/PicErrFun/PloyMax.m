function [a0,a1,a2,a3,tmax,gmax] = PloyMax(tEnd,gEnd,dotgEnd)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
tLE = tEnd(:,1);            tRE = tEnd(:,2);
gLE = gEnd(:,1);            gRE = gEnd(:,2);
dotgLE = dotgEnd(:,1);      dotgRE = dotgEnd(:,2);

a3 = (2*(gLE - gRE) + (tRE - tLE).*(dotgLE + dotgRE)) / (tRE - tLE).^3;
a2 = (3*(gRE - gLE) - (tRE - tLE).*(2*dotgLE + dotgRE)) / (tRE - tLE).^2;
a1 = dotgLE;
a0 = gLE;
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
end
if tRE<=tmax
    gmax = gRE;
end
if tLE<tmax & tmax<=tRE

    if 4*a2.^2 <= 12*a3.*a1
        gmax = (2*a2^3 - 9*a1*a2*a3 + 27*a0*a3^2)/(27*a3^2);
    else
        gmax = (27*a0*a3^2 + 2*(a2^2 - 3*a1*a3)^(3/2) + 2*a2^3 - 9*a1*a2*a3)/(27*a3^2);        
    end

end

end