import math
from decimal import Decimal

r_s = 0.8523
L_rbd = 200
amino_subs = 19

sum = 0
for k in range(0,4):
    N_k = math.pow(r_s,k)*(math.pow(19,k))*(math.comb(L_rbd-k,k))
    sum += N_k

print('%.2E' % Decimal(sum))