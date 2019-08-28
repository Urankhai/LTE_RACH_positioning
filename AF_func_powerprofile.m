function P = AF_func_powerprofile(Shifts, impulse_signal, udt, c, w)
Num_ant = size(impulse_signal,1);
Geom_shift = Shifts/c;
discret_part = ceil(Geom_shift/udt); % dt*N, where N is integer

discret_r = zeros(size(impulse_signal));
for lk = 1:Num_ant
    discret_r(lk,:) = [zeros(1,discret_part(lk)), impulse_signal(lk,1:end-discret_part(lk))];
end
nondcrt_part = udt*discret_part - Geom_shift;
Nondscrt_shift = exp(1j*w*nondcrt_part);
ant_sum = transpose(Nondscrt_shift)*discret_r;
P = ant_sum*ant_sum';