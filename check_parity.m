function [parity, decoded_word] = check_parity(precbits, word)

parity = 0;
decoded_word = '';
precbits_dec = zeros(1,2);
word_dec = zeros(1,24);
computed_parity = zeros(1,6);

for p = 1 : length(precbits)
    precbits_dec(p) = bin2dec(precbits(p));
end

for q = 1 : length(word)
    word_dec(q) = bin2dec(word(q));
end

seq25 = [1,2,3,5,6,10,11,12,13,14,17,18,20,23];
seq26 = [2,3,4,6,7,11,12,13,14,15,18,19,21,24];
seq27 = [1,3,4,5,7,8,12,13,14,15,16,19,20,22];
seq28 = [2,4,5,6,8,9,13,14,15,16,17,20,21,23];
seq29 = [1,3,5,6,7,9,10,14,15,16,17,18,21,22,24];
seq30 = [3,5,6,8,9,10,11,13,15,19,22,23,24];

for i = 1 : 24
    decoded_word_dec(i) = double(xor(word_dec(i),precbits_dec(end)));
    decoded_word = [decoded_word dec2bin(double(xor(word_dec(i),precbits_dec(end))))];
end

computed_parity(1) = precbits_dec(end-1);
for j = 1 : length(seq25)
    computed_parity(1) = double(xor(computed_parity(1),decoded_word_dec(seq25(j))));
end

computed_parity(2) = precbits_dec(end);
for j = 1 : length(seq26)
    computed_parity(2) = double(xor(computed_parity(2),decoded_word_dec(seq26(j))));
end

computed_parity(3) = precbits_dec(end-1);
for j = 1 : length(seq27)
    computed_parity(3) = double(xor(computed_parity(3),decoded_word_dec(seq27(j))));
end

computed_parity(4) = precbits_dec(end);
for j = 1 : length(seq28)
    computed_parity(4) = double(xor(computed_parity(4),decoded_word_dec(seq28(j))));
end

computed_parity(5) = precbits_dec(end);
for j = 1 : length(seq29)
    computed_parity(5) = double(xor(computed_parity(5),decoded_word_dec(seq29(j))));
end

computed_parity(6) = precbits_dec(end-1);
for j = 1 : length(seq30)
    computed_parity(6) = double(xor(computed_parity(6),decoded_word_dec(seq30(j))));
end

if (word_dec(25:30) == computed_parity(1:6))
    parity = 1;
end