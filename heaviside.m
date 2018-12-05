function [h] = heaviside(x)
    if(x<0)
        h = 0;
    end
    if(x == 0)
        h = 0.5;
    end
    if(x>0)
        h = 1;
    end