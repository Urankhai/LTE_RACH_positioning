function [preamble, seq] = AF_func_PRACH2(temp_ID)
N_ZC = 839;     % according to LTE specification
% we assume, the format of PRACH is 0, which means the length of sequence
% 24576, the length of CP is 3168, and GT length is 2976.
length_seq = 24576;
length_CP = 3168;
length_GT = 2976;

GF_below = 12;
GF_above = 13;
% all prime numbers up to 839
prime_numbers = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829];

% if we need to use one sequence for all users
ue_seq = 2;
if ue_seq == 1
    % generation of Zadoff Chu sequence
    pn = prime_numbers(1);
    seq = lteZadoffChuSeq(pn, N_ZC);
    % align for 6 RBs
    time_shift = 132*(temp_ID - 1);
    temp_seq = [seq(end-time_shift+1:end);seq(1:end-time_shift)];
%     keyboard
else
%     keyboard
    pn = prime_numbers(temp_ID);%randi(144,1));
    seq = lteZadoffChuSeq(pn, N_ZC);
    % align for 6 RBs
    temp_seq = seq;
end
% keyboard
fseq = [zeros(GF_above,1); fft(temp_seq)/sqrt(839); zeros(GF_below,1)]; % it takes 6 RBs
preamble = sqrt(length_seq)*ifft(fseq, length_seq);

% rach = [preamble(end-length_CP+1:end);preamble; zeros(length_GT,1)];
