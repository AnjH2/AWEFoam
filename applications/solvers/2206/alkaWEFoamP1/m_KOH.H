dimensionedScalar H1("H1",(dimMoles/dimMass)/(pow(dimMoles/dimVolume,3)*pow(dimTemperature,2)),1.28756e-9);
dimensionedScalar H2("H2",(dimMoles/dimMass)/(pow(dimMoles/dimVolume,3)*dimTemperature),4.57346e-7);
dimensionedScalar H3("H3",(dimMoles/dimMass)/(pow(dimMoles/dimVolume,3)),8.51811e-4);
dimensionedScalar I1("I1",(dimMoles/dimMass)/(pow(dimMoles/dimVolume,2)*pow(dimTemperature,2)),-1.87305e-8);
dimensionedScalar I2("I2",(dimMoles/dimMass)/(pow(dimMoles/dimVolume,2)*dimTemperature),3.52482e-6);
dimensionedScalar I3("I3",(dimMoles/dimMass)/(pow(dimMoles/dimVolume,2)),0.00764237);
dimensionedScalar J1("J1",(dimMoles/dimMass)/((dimMoles/dimVolume)*pow(dimTemperature,2)),3.09226e-7);
dimensionedScalar J2("J2",(dimMoles/dimMass)/((dimMoles/dimVolume)*dimTemperature),3.67376e-4);
dimensionedScalar J3("J3",(dimMoles/dimMass)/((dimMoles/dimVolume)),0.869801);

volScalarField cal_mKOH(const volScalarField& T, const volScalarField& M) {
	return (H1*sqr(T)+H2*T+H3)*pow(M,3)+(I1*sqr(T)+I2*T+I3)*pow(M,2)+(J1*sqr(T)+J2*T+J3)*M;
}

