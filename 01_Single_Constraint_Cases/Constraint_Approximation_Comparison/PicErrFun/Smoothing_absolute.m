function [sms_rlst] = Smoothing_absolute(input)
%Actually function is different from version described in loxton's article,
%because inequality in problem reformulation is in the opposite direction.
%%global epislon

epislon=1e-3;


sms_rlst=sqrt(input.^2+epislon^2);
end



