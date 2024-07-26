function out = roundUp(in)
% Function to round up to the first decimal place

if in>=1
    out = ceil(in);
else
    count = 0;
    test = in;
    while test<1
        test = test*10;
        count = count + 1;
    end
    out = ceil(in*10^(count))/10^(count);
end
if out == in % In this case there was nothing to round up
    out = 1.1*in; % increase size by 10% to get a slightly larger number
end
end 


