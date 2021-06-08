# **Kessler94**

# symbols
- T
- qv 
- qc 
- qp 

# coefficients
- k1 = coeffs_["condensation"]
- k2 = coeffs_["autoconversion"]
- k3 = coeffs_["accretion"]
- k4 = coeffs_["evaporation"]

# verbatim
~~~C++
Thermodynamics *pthermo = pmy_block->pthermo;
pthermo->SaturationSurplus(dqsat_.data(), q, VariableType::chem);
Real dq = dqsat_[iqv];
Real qs = qv - dq;
Real Rv = pthermo->GetRd()/pthermo->GetMassRatio(iqv);
Real lf = pthermo->GetLatent(iqc,T)/(Rv*T) - 1.;
Real dqsdt = lf*qs/T;
~~~

# simple reactions
- qc -> qp ; k2
- qc + qp -> 2qp ; k3
- qp -> qv ; k4 ; -H1 | qv < qs

# custom reactions
- qv -> qc ; k1*(qv - qs(T)) ; H1 | qv > qs
- qc -> qv ; k1*qc*(qs(T) - qv) ; -H1 | qv < qs

# relations
- Derivative(qs(T), T) -> dqsdt
- qs(T) -> qs
