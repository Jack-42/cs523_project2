import math
from decimal import Decimal

'''
calculate the size of a k-step neutral network
'''

r_s = 0.8523
L_rbd = 200
subs = 19

steps = 3

sum = 0
for k in range(0,steps+1):
    N_k = math.pow(r_s,k)*(math.pow(subs,k))*(math.comb(L_rbd-k,k))
    sum += N_k

print('%.2E' % Decimal(sum))