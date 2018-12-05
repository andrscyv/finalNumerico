function d = deriva(v, h)
    d(1) = (v(2)-v(1))/h;
    tam = length(v);
    for i = 2:tam-1
        d(i) = (v(i+1)-v(i-1))/h;
    end
    d(tam) = (v(tam)-v(tam-1))/h;
