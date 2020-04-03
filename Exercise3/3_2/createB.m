function B = createB(index,b)
    if(index == 1)
        B = [0,0;0,0];
    elseif(index ==2)
        B = [0,0;b,0];
    else 
        B = [0,0;0,b];
    end
end

