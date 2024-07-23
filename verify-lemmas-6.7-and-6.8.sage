from polylog.sage import *

O = ThetaSharpOperator(6, 2)
R = O.OU_phi_algebra
S = O.OU_S_ring
gens = ['log'] + [f'Li{n}' for n in range(1,7)] + list(R.gens())
degs = len(list(R.gens()))*[1]+ [1,1,2,3,4,5,6]
LR = PolynomialRing(S, len(gens), gens, order=TermOrder('wdegrevlex', degs))

Li = [LR(O.log())] + [LR(O.Li(i)) for i in range(1,7)]

#Elimanate phisgma5 and phisgma3 from Li6
f6 = Li[5].coefficient(LR('phisigma5'))*(LR('Li6') - Li[6]) - Li[6].coefficient(LR('phisigma5'))*(LR('Li5')-Li[5])
f6 = f6.coefficient(LR('phisigma3'))*(LR('Li3') - Li[3]) - Li[3].coefficient(LR('phisigma3'))*f6

#Elimanate phisgma3 from Li4
f4 = Li[4].coefficient(LR('phisigma3'))*(LR('Li3')-Li[3]) - Li[3].coefficient(LR('phisigma3'))*(LR('Li4')-Li[4])
# [Li_3 - theta#(Li3)] * (phi3 coefficient of Li4) - [Li_4 - theta#(Li4)] * (phi3 coefficient of Li3)

#Find the correct substitutions for phi1t0 and phi1t1
M = Matrix(LR, [[Li[1].coefficient(LR('phi1t0')),Li[1].coefficient(LR('phi1t1'))],[Li[2].coefficient(LR('phi1t0')),Li[2].coefficient(LR('phi1t1'))]])
(phi0sub, phi1sub) = M.adjugate()*vector(LR, [LR('Li1'), LR('Li2')])/M.determinant()

#Substitute into f4 and f6. Clear M.determinant() from the denominator at the same time.
f6 = M.determinant()*f6.subs({LR('phi1t0'):LR.fraction_field()(phi0sub), LR('phi1t1'):LR.fraction_field()(phi1sub)})
f4 = M.determinant()*f4.subs({LR('phi1t0'):LR.fraction_field()(phi0sub), LR('phi1t1'):LR.fraction_field()(phi1sub)})

#Make phi0t1 the subject of log
phi0t1sub = LR.fraction_field()((LR('log') - LR('Sa*phi0t0'))/LR('Sb'))

#Substitue into f4 and f6, and clear denominators
f6 = LR(LR('Sb')^4*f6.subs({LR('phi0t1'):phi0t1sub})) #Degree 4 in phi0t0
f4 = LR(LR('Sb')^2*f4.subs({LR('phi0t1'):phi0t1sub})) #Degree 2 in phi0t0

print('Degree of v4: ',f4.polynomial(LR('phi0t0')).degree()) #2
print('Degree of v6: ',f6.polynomial(LR('phi0t0')).degree()) #4

print('Is leading coefficient of v4 divisible by log*Li1 - 2*Li2: ', 
      f4.polynomial(LR('phi0t0')).coefficients()[-1] % LR('log*Li1 - 2*Li2') == 0) #True

print('Is leading coefficient of v6 divisible by log*Li1 - 2*Li2: ', 
      f6.polynomial(LR('phi0t0')).coefficients()[-1] % LR('log*Li1 - 2*Li2') == 0) #True
