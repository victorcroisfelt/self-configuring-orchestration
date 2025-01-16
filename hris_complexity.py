C = 8
N = 32
K = 4

pow = K * (2 * N + C + 3) + C
dsp = N * (K**3 + 5 * K + C) + 3*K


print(pow)
print(dsp)

pow_pd = 0.782805 
dsp_pd = 0.818438 