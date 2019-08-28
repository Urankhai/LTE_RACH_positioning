function [UL_r, discret_shift, nondcrt_part, Nondscrt_shift] = AF_func_MIMO_rec(w, uprach, udt, time_delay)

discret_shift = ceil(time_delay/udt);
discret_r = [zeros(1,discret_shift), uprach(1:end - discret_shift)];

discret_part = udt*discret_shift;

nondcrt_part = discret_part - time_delay;
Nondscrt_shift = exp(1j * w * nondcrt_part);
% keyboard
%                  exp(1j * w * udt * samples)

UL_r = Nondscrt_shift*discret_r;