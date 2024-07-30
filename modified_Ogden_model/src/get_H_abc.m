% refer to eqn. 2.15 of Liu, Guan, Zhao & Luo 2024 preprint
function out = get_H_abc(Na, Nb, Nc)
out = zeros(3,3,3,3,3,3);
out = out + cross_otimes_1d_to_6d(Na, Nb, Nc, Na, Nb, Nc);
out = out + cross_otimes_1d_to_6d(Na, Nb, Nc, Na, Nc, Nb);

out = out + cross_otimes_1d_to_6d(Na, Nb, Na, Nc, Nb, Nc);
out = out + cross_otimes_1d_to_6d(Na, Nb, Na, Nc, Nc, Nb);

out = out + cross_otimes_1d_to_6d(Nb, Na, Nc, Na, Nb, Nc);
out = out + cross_otimes_1d_to_6d(Nb, Na, Nc, Na, Nc, Nb);

out = out + cross_otimes_1d_to_6d(Nb, Na, Na, Nc, Nb, Nc);
out = out + cross_otimes_1d_to_6d(Nb, Na, Na, Nc, Nc, Nb);
end