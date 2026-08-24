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
};
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
};
volScalarField dlnGammaOH_dlnCOH
(
    const volScalarField& T,
    const volScalarField& C_OH   // [mol m^-3]
)
{
    // ---- 3.1 convenience constants ----------------------------------------
    const dimensionedScalar w(dimMoles/dimMass, 55.5093);     // mol water / kg
    const dimensionedScalar b(dimless,          0.03603);
    const dimensionedScalar one(dimless,        1.0);
    const dimensionedScalar oneOver1000(dimless,1.0/1000.0);

    // ---- 3.2 basic substitutions ------------------------------------------
    volScalarField M   = C_OH * oneOver1000;           // mol kg^-1
    volScalarField a3  = H1*sqr(T) + H2*T + H3;
    volScalarField a2  = I1*sqr(T) + I2*T + I3;
    volScalarField a1  = J1*sqr(T) + J2*T + J3;

    volScalarField m   = a3*pow3(M) + a2*sqr(M) + a1*M;                // f(M,T)
    volScalarField dmdM = 3.0*a3*sqr(M) + 2.0*a2*M + a1;               // f′(M,T)

    volScalarField x   = (2.0*m)/(w + 2.0*m);
    volScalarField dxdm= (2.0*w)/sqr(w + 2.0*m);

    // ---- 3.3 dg/dx ---------------------------------------------------------
    volScalarField dgdx =
          (K2 + 2.0*K3*x + 3.0*K4*sqr(x) + 4.0*K5*pow3(x))
        + (L2 + 2.0*L3*x + 3.0*L4*sqr(x) + 4.0*L5*pow3(x))*T
        + (M2 + 2.0*M3*x + 3.0*M4*sqr(x) + 4.0*M5*pow3(x))*sqr(T);

    // ---- 3.4 final factor inside the brackets in (1) ----------------------
    volScalarField bracket = dgdx*dxdm - b / (one + b*m);

    // ---- 3.5 elasticity ----------------------------------------------------
    return M * dmdM * bracket;          // d ln γ  / d ln C
};
volScalarField dlnGammaOH_dlnmKOH
(
    const volScalarField& T,
    const volScalarField& m_KOH   // mol kg^-1  (already molality!)
)
{
    // ---- constants --------------------------------------------------------
    const dimensionedScalar w(dimMoles/dimMass, 55.5093);
    const dimensionedScalar b(dimless,          0.03603);
    const dimensionedScalar one(dimMoles/dimMass,        1.0);

    // ---- mole fraction of KOH --------------------------------------------
    volScalarField x    = (2.0*m_KOH)/(w + 2.0*m_KOH);
    volScalarField x_m  = (2.0*w)/sqr(w + 2.0*m_KOH);      // ∂x/∂m

    // ---- dg/dx (identical to previous) ------------------------------------
    volScalarField dgdx =
          (K2 + 2.0*K3*x + 3.0*K4*sqr(x) + 4.0*K5*pow3(x))
        + (L2 + 2.0*L3*x + 3.0*L4*sqr(x) + 4.0*L5*pow3(x))*T
        + (M2 + 2.0*M3*x + 3.0*M4*sqr(x) + 4.0*M5*pow3(x))*sqr(T);

    // ---- bracket term -----------------------------------------------------
    volScalarField bracket = dgdx * x_m - b / (one + b*m_KOH);

    // ---- final elasticity -------------------------------------------------

    return m_KOH * bracket;          // d ln γ  / d ln m
};
