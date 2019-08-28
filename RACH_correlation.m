clear all
close all

% In LTE RACH, there are used ZC sequences with the length 839, and prime
% numbers are used as root numbers
prime_numbers = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829];
Nprime = length(prime_numbers);

corr_mean_value = zeros(Nprime);
min_mean_value = 1000;
corr_min_value = zeros(Nprime);
min_min_value = 1000;
for temp_ID = 1: Nprime-1
    
    [rach, seq] = AF_func_PRACH(temp_ID);
    
    for l = temp_ID+1:Nprime
        [temp_rach, temp_seq] = AF_func_PRACH(l);
        [temp_out_s, temp_out_t] = AF_func_correlator(seq, transpose(temp_rach));
        mean_corr = mean(abs(temp_out_s));
        corr_mean_value(temp_ID, l) = mean_corr;
        if mean_corr < min_mean_value
            min_mean_value = mean_corr;
            min_mean_ID = temp_ID;
            min_mean_l = l;
        end
        
        min_corr = min(abs(temp_out_s));
        corr_min_value(temp_ID, l) = min_corr;
        if min_corr < min_min_value
            min_min_value = min_corr;
            min_min_ID = temp_ID;
            min_min_l = l;
        end
        
    end
    
end

disp(['min_min_value = ', num2str(min_min_value), '; y position = ', num2str(min_min_ID), '; x position = ', num2str(min_min_l)])
disp(['min_mean_value = ', num2str(min_mean_value), '; y position = ', num2str(min_mean_ID), '; x position = ', num2str(min_mean_l)])


figure
surf(corr_mean_value);

figure
surf(corr_min_value);

% the result of the correlation shows that sequences [126, 12, 97],
% [144,54, 81], and [123, 41, 115] give the smallest mutual
% crosscorrelation.
