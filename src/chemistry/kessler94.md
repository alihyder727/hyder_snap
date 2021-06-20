# **Kessler94**

# variables
- T
- qv
- qc
- qp 

# coefficients
- k1 = coeffs_["condensation"]
- k2 = coeffs_["autoconversion"]
- k3 = coeffs_["accretion"]
- k4 = coeffs_["evaporation"]

# simple reactions
- qc -> qp ; k2
- qc + qp -> 2qp ; k3
- qc -> qv ; k1*(qs(T) - qv)/mols | qv < qs
- qp -> qv ; k4*(qs(T) - qv)/mols | qv < qs

# custom reactions
- qv -> qc ; k1*(qv - qs(T)) | qv > qs

# relations
- Derivative(qs(T), T) -> dqsdt
- qs(T) -> qs

# verbatim
~~~C++
pthermo->SaturationSurplus(dqsat_.data(), q, VariableType::chem);
Real dq = dqsat_[iqv];
Real qs = qv - dq;
Real Rv = pthermo->GetRd()/pthermo->GetMassRatio(iqv);
Real lf = pthermo->GetLatent(iqc,T) - Rv*T;
Real dqsdt = qs/T*lf/(Rv*T);
Real mols = q[IPR]/(Thermodynamics::Rgas*q[IDN]);
~~~
