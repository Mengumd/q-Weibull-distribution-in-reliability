u = 0;
p = 0;

while(u < 0.99)
    u = rand;
    p = p + u;
    
    if (p > 5)
        break;
    end
end

sol = p;