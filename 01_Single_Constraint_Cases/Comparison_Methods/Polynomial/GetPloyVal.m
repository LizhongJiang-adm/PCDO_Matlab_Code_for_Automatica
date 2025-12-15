function PolyVal = GetPloyVal(tEnd,gEnd,dotgEnd,t)
% 取出近似多项式在t上的值

%   此处提供详细说明
tLE = tEnd(:,1);            tRE = tEnd(:,2);
gLE = gEnd(:,1);            gRE = gEnd(:,2);
dotgLE = dotgEnd(:,1);      dotgRE = dotgEnd(:,2);

a3 = (2*(gLE - gRE) + (tRE - tLE).*(dotgLE + dotgRE)) / (tRE - tLE).^3;
a2 = (3*(gRE - gLE) - (tRE - tLE).*(2*dotgLE + dotgRE)) / (tRE - tLE).^2;
a1 = dotgLE;
a0 = gLE;

PolyVal = a3*(t-tLE).^3 + a2*(t-tLE).^2 + a1*(t-tLE) + a0;

% PloyApx = @(t) a3.*t.^3 + a2.*t.^2 + a1.*t + a0;

end