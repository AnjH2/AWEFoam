IOdictionary transportProperties
(
	IOobject
	(
		"transportProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);
/*Info << "Assigning value from dictionary sigma_s_ref \n" << endl;
dimensionedScalar sigma_s_ref
(
	"sigma_s_ref",
	dimensionSet( -1, -3, 3, 0, 0, 2,0),
	transportProperties
);
Info << "Assigning value from dictionary sigma_e \n" << endl;
dimensionedScalar sigma_e
(
	"sigma_e_ref",
	dimensionSet( -1, -3, 3, 0, 0, 2, 0),
	transportProperties
); 

Info << "Assigning value from dictionary a_Ne \n" << endl;
dimensionedScalar a_Ne
(
	"a_Ne",
	dimensionSet( 0, -1, 0, 0, 0, 0, 0),
	transportProperties
);

Info << "Assigning value from dictionary a_Pe \n" << endl;
dimensionedScalar a_Pe
(
	"a_Pe",
	dimensionSet( 0, -1, 0, 0, 0, 0, 0),
	transportProperties
);
Info << "Assigning value from dictionary i0_Ne_ref \n" << endl;
dimensionedScalar i0_ref_Ne
(
	"i0_ref_Ne",
	dimensionSet( 0, -2, 0, 0, 0, 1, 0),
	transportProperties
);
Info << "Assigning value from dictionary i0_Pe_ref \n" << endl;
dimensionedScalar i0_ref_Pe
(
	"i0_ref_Pe",
	dimensionSet( 0, -2, 0, 0, 0, 1, 0),
	transportProperties
);
Info << "Assigning value from dictionary alpha_a_Ne \n" << endl;
dimensionedScalar alpha_a_Ne
(
	"alpha_a_Ne",
	dimensionSet( 0, 0, 0, 0, 0, 0, 0),
	transportProperties
);
Info << "Assigning value from dictionary alpha_a_Ne \n" << endl;
dimensionedScalar alpha_c_Ne
(
	"alpha_c_Ne",
	dimensionSet( 0, 0, 0, 0, 0, 0, 0),
	transportProperties
);
Info << "Assigning value from dictionary alpha_a_Pe \n" << endl;
dimensionedScalar alpha_a_Pe
(
	"alpha_a_Pe",
	dimensionSet( 0, 0, 0, 0, 0, 0, 0),
	transportProperties
);
Info << "Assigning value from dictionary alpha_a_Pe \n" << endl;
dimensionedScalar alpha_c_Pe
(
	"alpha_c_Pe",
	dimensionSet( 0, 0, 0, 0, 0, 0, 0),
	transportProperties
);*/
/*Info << "Assigning value form dictionary F \n" << endl;
dimensionedScalar F
(
	"F",
	dimensionSet( 0, 0, 1, 0, -1, 1, 0),
	transportProperties
);
Info << "Assigning value from dictionary R \n" << endl;
dimensionedScalar R
(
	"R",
	dimensionSet( 1, 2, -2, -1, -1, 0, 0),
	transportProperties
);*/
//dimensionedScalar R (Foam::constant::physicoChemical::R);
//dimensionedScalar F (Foam::constant::physicoChemical::F);
const dimensionedScalar& R = Foam::constant::physicoChemical::R;
const dimensionedScalar& F = Foam::constant::physicoChemical::F;
/*
Info << "Assigning value from dictionary ENe \n" << endl;
dimensionedScalar E0Ne
(
	"E0Ne",
	dimensionSet( 1, 2, -3, 0, 0, -1, 0),
	transportProperties
);
Info << "Assigning value from dictionary EPe \n" << endl;
dimensionedScalar E0Pe
(
	"E0Pe",
	dimensionSet( 1, 2, -3, 0, 0, -1, 0),
	transportProperties
);

Info << "Assigning value from dictionary Ea_Ne \n" << endl;
dimensionedScalar Ea_Ne
(
	"Ea_Ne",
	dimensionSet( 1, 2, -2, 0, -1, 0, 0),//[J/(mol)]
	transportProperties
);
Info << "Assigning value from dictionary Ea_Pe \n" << endl;
dimensionedScalar Ea_Pe
(
	"Ea_Pe",
	dimensionSet( 1, 2, -2, 0, -1, 0, 0),//[J/(mol)]
	transportProperties
);
*/
Info << "Assigning value from dictionary T0 \n" << endl;
dimensionedScalar T0
(
	"T0",
	dimensionSet( 0, 0, 0, 1, 0, 0, 0),
	transportProperties
);
/*Info << "Assigning value from dictionary T_ref_Ne \n" << endl;
dimensionedScalar T_ref_Ne
(
	"T_ref_Ne",
	dimensionSet( 0, 0, 0, 1, 0, 0, 0),
	transportProperties
);
Info << "Assigning value from dictionary T_ref_Pe \n" << endl;
dimensionedScalar T_ref_Pe
(
	"T_ref_Pe",
	dimensionSet( 0, 0, 0, 1, 0, 0, 0),
	transportProperties
);


Info << "Assigning value from dictionary Ds_Ne \n" << endl;
dimensionedScalar Ds_Ne
(
	"Ds_Ne",
	dimensionSet( 1, 2, -2, -1, -1, 0, 0),//[J/(mol*K)]
	transportProperties
);
Info << "Assigning value from dictionary Ds_Pe \n" << endl;
dimensionedScalar Ds_Pe
(
	"Ds_Pe",
	dimensionSet( 1, 2, -2, -1, -1, 0, 0),//[J/(mol*K)]
	transportProperties
);


Info << "Assigning value from dictionary C_KC \n" << endl;
dimensionedScalar C_KC
(
	"C_KC",
	dimensionSet( 0, 0, 0, 0, 0, 0, 0),//[m]
	transportProperties
);*/
Info << "Assigning value from dictionary mu_KOH \n" << endl;
dimensionedScalar mu_KOH
(
	"mu_KOH",
	dimensionSet( 1, -1, -1, 0, 0, 0, 0),//[m]
	transportProperties
);
/*dimensionedScalar nu
(
    	"nu_KOH",
    	dimViscosity,
    	transportProperties
);*/



