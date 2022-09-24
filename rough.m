%Lbsh definition
    if h_b > h_r
        l_bsh = -18*log10(1 + dh_b);
    else
        l_bsh = 0;
    end

%ka definition
    k_a = ones(1,length(f));
    for i=1:length(f)
        if h_b > h_r && f(i) > 2000
            k_a(i) = 71.4;
        elseif (h_b < h_r || h_b == h_r) && f(i) > 2000 && (d > 500 || d == 500)
            k_a(i) = 73 - 0.8*dh_b;
        elseif (h_b < h_r || h_b == h_r) && f(i) > 2000 && d < 500
            k_a(i) = 73 - 1.6*dh_b*(d/1000);
        elseif h_b > h_r && (f(i) < 2000 || f(i) == 2000)
            k_a(i) = 54;
        elseif (h_b < h_r || h_b == h_r) && (f(i) < 2000 || f(i) == 2000) && (d > 500 || d == 500)
            k_a(i) = 54 - 0.8*dh_b;
        elseif (h_b < h_r || h_b == h_r) && (f(i) < 2000 || f(i) == 2000) && d < 500
            k_a(i) = 54 - 1.6*dh_b*(d/1000);
        end
    end

%kd definition
    if h_b > h_r
        k_d = 18;
    else
        k_d = 18 - 15*(dh_b/h_r);
    end

%kf definition
    k_f = ones(1,length(f));
    for i=1:length(f)
        if f(i) > 2000
            k_f(i) = -8;
        else 
            k_f(i) = -4 + 0.7*((f(i)/925) -1);
        end
    end

    l1_msd = ones(1,length(f));
    for i=1:length(f)
        l1_msd = l_bsh + k_a + k_d*log10(d/1000) + k_f*log10(f(i)) - 9*log10(b);
    end