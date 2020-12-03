% takes values between 0 and 1, will be multiplied by factor later
function [y] = switch_rate_of_glucose(G)
y=[];
for g=G
    if (g > 5) || (g < 0.01)
    %if (g > 5)
        y = [y;0];
    else
        y = [y;1];
    end
end

