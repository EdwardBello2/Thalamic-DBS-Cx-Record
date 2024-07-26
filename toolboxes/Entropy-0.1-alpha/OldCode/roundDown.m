function out = roundDown(in)
% Function to round down to the first decimal place

if in>=1
    out = floor(in);
else
    count = 0;
    test = in;
    while test<1
        test = test*10;
        count = count + 1;
    end
    out = floor(in*10^(count))/10^(count);
end
if out == in % In this case there was nothing to round down
    out = 0.9*in; % Decrease size by 10% to get a slightly smaller number
end
end