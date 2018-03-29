function [next_index, max_diff] = max_difference(x, z_f, z_h1, z_h2, g_f, g_h1, g_h2, A, b, lambda1, R, lambda2, index, update_type)

switch index 
    case 0
        h1_diff = yh1(x, lambda1) - (yh1(z_h1, lambda1) + (g_h1')*(x - z_h1));
        h2_diff = yh2(x, R, lambda2) - (yh2(z_h2, R, lambda2) + (g_h2')*(x - z_h2));                  
        if (h1_diff >= h2_diff)
            next_index = 1;           
            max_diff = h1_diff;
        else
            next_index = 2;
            max_diff = h2_diff;
        end
        
        if (update_type == 1)
            max_diff = h1_diff + h2_diff;
        end
    case 1
        f_diff = yf(x, A, b) - (yf(z_f, A, b) + (g_f')*(x - z_f));
        h2_diff = yh2(x, R, lambda2) - (yh2(z_h2, R, lambda2) + (g_h2')*(x - z_h2));
        if (f_diff >= h2_diff)
            next_index = 0;
            max_diff = f_diff;
        else
            next_index = 2;
            max_diff = h2_diff;
        end
        
        if (update_type == 1)
            max_diff = f_diff + h2_diff;
        end
    case 2
        f_diff = yf(x, A, b) - (yf(z_f, A, b) + (g_f')*(x - z_f));
        h1_diff = yh1(x, lambda1) - (yh1(z_h1, lambda1) + (g_h1')*(x - z_h1));
        if (f_diff >= h1_diff)
            next_index = 0;
            max_diff = f_diff;
        else
            next_index = 1;
            max_diff = h1_diff;
        end      
        
        if (update_type == 1)
            max_diff = f_diff + h1_diff;
        end
end
end