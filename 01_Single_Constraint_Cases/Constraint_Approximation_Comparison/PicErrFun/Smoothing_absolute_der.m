function [sms_rlst] = Smoothing_absolute_der(input)
%Actually function is different from version described in loxton's article,
%because inequality in problem reformulation is in the opposite direction.
%%global epislon

epislon=1e-3;

% epislon.^2 + input.^2
sms_rlst=input./(epislon.^2 + input.^2).^(1/2);
end



