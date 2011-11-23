
nStates = 3
absorbing_cond = [0,0.1,0.1,inf] # 0-fully absorb, inf- reflective
std_state = [0.33 ,0.33,0.33]

for i in arange(nStates)
  koni = smol(absorbing_cond[i])
  kon  = kon + (std_state[i] * koni)

