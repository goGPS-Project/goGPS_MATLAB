function mbr = mbrot()
    i_list = 1:1024;
    mbr = nan(1024,1024,3);
    for i = i_list(:)'
        for j = i_list(:)'
            
            mbr(i,j,1) = red(i,j);
            mbr(i,j,2) = green(i,j);
            mbr(i,j,3) = blue(i,j);
        end
    end
    
    figure; image(mbr./max(mbr(:)));
    
    function out = red(i, j)
        a = 0;
        b = 0;
        d = 0;
        n = 0;
        while (a*a+(b*b)<4 && n<1024)
            b=2*a*b+(j-1)/5e4+.06;
            a=a*a-d+(i-1)/5e4+.34;
            d = b*b;
            n = n+1;
        end
        out = n/4;
    end
    
    function out = green(i, j)
        out = 4 * red(i, j);
    end
    
    function out = blue(i, j)
        out = 6 * red(i, j);
    end
end
