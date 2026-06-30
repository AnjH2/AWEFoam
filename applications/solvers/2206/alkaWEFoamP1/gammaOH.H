dimensionedScalar K1("K1",(dimless),-0.0151895);
dimensionedScalar K2("K2",(dimless),-46.8743);
dimensionedScalar K3("K3",(dimless),343.766);
dimensionedScalar K4("K4",(dimless),-566.491);
dimensionedScalar K5("K5",(dimless),472.221);
dimensionedScalar L1("L1",(dimless)/dimTemperature,4.34773e-5);
dimensionedScalar L2("L2",(dimless)/dimTemperature,0.258126);
dimensionedScalar L3("L3",(dimless)/dimTemperature,-1.38799);
dimensionedScalar L4("L4",(dimless)/dimTemperature,1.89206);
dimensionedScalar L5("L5",(dimless)/dimTemperature,-1.39535);
dimensionedScalar M1("M1",(dimless)/pow(dimTemperature,2),-1.10721e-7);
dimensionedScalar M2("M2",(dimless)/pow(dimTemperature,2),-4.37513e-4);
dimensionedScalar M3("M3",(dimless)/pow(dimTemperature,2),0.00228340);
dimensionedScalar M4("M4",(dimless)/pow(dimTemperature,2),-0.00398424);
dimensionedScalar M5("M5",(dimless)/pow(dimTemperature,2),0.00344981);

volScalarField cal_gammaOH(const volScalarField& T, const volScalarField& m_KOH) {
	volScalarField x_KOH=2*m_KOH/(dimensionedScalar(dimMoles/dimMass,55.5093)+2*m_KOH);
	return exp(
		   K1 + K2*x_KOH + K3*pow(x_KOH,2) + K4*pow(x_KOH,3) + K5*pow(x_KOH,4) 
		+( L1 + L2*x_KOH + L3*pow(x_KOH,2) + L4*pow(x_KOH,3) + L5*pow(x_KOH,4) )*T
		+( M1 + M2*x_KOH + M3*pow(x_KOH,2) + M4*pow(x_KOH,3) + M5*pow(x_KOH,4) )*pow(T,2) 
		)/ (dimensionedScalar(dimMoles/dimMass,1) + 0.03603*m_KOH );
}

