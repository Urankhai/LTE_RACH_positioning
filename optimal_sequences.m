clear all
close all

prime_numbers = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829];

N = length(prime_numbers);

Max_corr_matrix = zeros(N);
figure
hold on
step = 0;
for i = 1:N
    seq1 = lteZadoffChuSeq(prime_numbers(i),839);
    fft_seq1 = fft(seq1);
    for j = i:N
        
        seq2 = lteZadoffChuSeq(prime_numbers(j),839);
        fft_seq2 = fft(seq2);
        fft_corr = fft_seq1.*conj(fft_seq2);
        out_s = ifft(fft_corr);
        if j > i
            step = step + 1;
        Max_corr_matrix(i,j) = max(abs(out_s));
        plot(step, max(abs(out_s)),'o')
        end
    end
end

figure
surf(Max_corr_matrix)

        