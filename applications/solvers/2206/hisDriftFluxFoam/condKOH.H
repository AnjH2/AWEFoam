dimensionedScalar A("A",dimless/(dimMoles/dimVolume),-2.041);
dimensionedScalar B("B",dimless/pow(dimMoles/dimVolume,2),-0.0028);
dimensionedScalar C("C",dimless/(dimMoles/dimVolume*dimTemperature),0.005332);
dimensionedScalar D("D",dimTemperature/(dimMoles/dimVolume),207.2);
dimensionedScalar E("E",dimless/pow3(dimMoles/dimVolume),0.001043);
dimensionedScalar F("F",dimless/(pow2(dimMoles/dimVolume)*pow2(dimTemperature)),-0.0000003);
dimensionedScalar G("G",dimensionSet(-1,-3,3,0,0,2,0),1);

volScalarField calKappa(volScalarField T, volScalarField M) {
	return (A*M+B*pow(M,2)+C*M*T+D*M/T+E*pow(M,3)+F*pow(M,2)*pow(T,2))*G*100;
}

