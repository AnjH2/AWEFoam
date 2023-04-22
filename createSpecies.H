IOdictionary speciesProperties
(
	IOobject
	(
		"speciesProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

//Determine number of species concentration fields***********************************


scalar nc=0;
std::vector<std::string> species;

// Search for list of objects at startTime
IOobjectList objects(mesh,"0");

// Search list of objects for volScalarFields
IOobjectList scalarFields
(
	objects.lookupClass("volScalarField")
);
for
(
	IOobjectList::iterator scalarFieldIter=scalarFields.begin();
	scalarFieldIter !=scalarFields.end();
	++scalarFieldIter
)
{
	//Read field
	volScalarField field
	(
		*scalarFieldIter(),
		mesh
	);
	word fieldname=field.name();
	if (fieldname.find("C")==0)
	{
		nc++;
		Info << "Found " << field.name() << endl;
		species.push_back(field.name());
	}
}
Info << "Number of Species/" << nc << endl;
//Create species concentration fields
PtrList <volScalarField>C(nc);
for (label i=0; i<C.size(); i++)
{
	//erase C from species list to make it universal
	species[i].erase(species[i].begin());
	// start creating fields.
	word fieldName="C"+species[i];
	Info << "Reading field " << fieldName << endl;
	
	C.set //-HM‘set or hook’(OF-version dependent)
	(
	i,
	new volScalarField
		(
		IOobject
			(
			fieldName,
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
			),
		mesh
		)
	);
}
//start reading scalars
PtrList <dimensionedScalar>D_ref(nc);
for (label i=0; i<D_ref.size(); i++)
{
	word fieldName="D"+species[i]+"_ref";
	Info << "Reading field " << fieldName << endl;
	D_ref.set
	
	(
		i,
		new dimensionedScalar
		(
			fieldName,
			dimensionSet( 0, 2, -1, 0, 0, 0, 0),//[J/(mol)]
			speciesProperties
		)
	);
}

PtrList <dimensionedScalar>z(nc);
for (label i=0; i<z.size(); i++)
{
	word fieldName="z"+species[i];
	Info << "Reading field " << fieldName << endl;
	z.set
	
	(
		i,
		new dimensionedScalar
		(
			fieldName,
			dimensionSet( 0, 0, 0, 0, 0, 0, 0),//[J/(mol)]
			speciesProperties
		)
	);
}

PtrList <dimensionedScalar>Ea(nc);
for (label i=0; i<Ea.size(); i++)
{
	word fieldName="Ea"+species[i];
	Info << "Reading field " << fieldName << endl;
	Ea.set
	
	(
		i,
		new dimensionedScalar
		(
			fieldName,
			dimensionSet( 1, 2, -2, 0, -1, 0, 0),//[J/(mol)]
			speciesProperties
		)
	);
}

PtrList <dimensionedScalar>T_ref(nc);
for (label i=0; i<T_ref.size(); i++)
{
	word fieldName="T_ref"+species[i];
	Info << "Reading field " << fieldName << endl;
	T_ref.set
	
	(
		i,
		new dimensionedScalar
		(
			fieldName,
			dimensionSet( 0, 0, 0, 1, 0, 0, 0),//[J/(mol)]
			speciesProperties
		)
	);
}

PtrList <dimensionedScalar>C_ref(nc);
for (label i=0; i<C_ref.size(); i++)
{
	word fieldName="C"+species[i]+"_ref";
	Info << "Reading field " << fieldName << endl;
	C_ref.set
	
	(
		i,
		new dimensionedScalar
		(
			fieldName,
			dimensionSet( 0, -3, 0, 0, 1, 0, 0),//[mol/(m3)]
			speciesProperties
		)
	);
}

PtrList <dimensionedScalar>MW(nc);
for (label i=0; i<MW.size(); i++)
{
	word fieldName="MW"+species[i];
	Info << "Reading field " << fieldName << endl;
	MW.set
	
	(
		i,
		new dimensionedScalar
		(
			fieldName,
			dimensionSet( 1, 0, 0, 0, -1, 0, 0),//[kg/mol]
			speciesProperties
		)
	);
}
Info << "Assigning value from dictionary MW_K \n" << endl;
dimensionedScalar MW_K
(
	"MW_K",
	dimensionSet( 1, 0, 0, 0, -1, 0, 0),
	speciesProperties
);
//start reading volume scalar fields.
PtrList <volScalarField>D(nc);
for (label i=0; i<D.size(); i++)
{
	// start creating fields.
	word fieldName="D"+species[i];
	Info << "Reading field " << fieldName << endl;
	
	D.set //-HM‘set or hook’(OF-version dependent)
	(
	i,
	new volScalarField
		(
		IOobject
			(
			fieldName,
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
			),
		D_ref[i]*exp(-Ea[i]/R*(1/T-1/T_ref[i]))
		)
	);
}

//Create species effective diffusion fields
PtrList <volScalarField>D_eff(nc);
for (label i=0; i<D_eff.size(); i++)
{
	// start creating fields.
	word fieldName="D"+species[i]+"_eff";
	Info << "Reading field " << fieldName << endl;
	
	D_eff.set //-HM‘set or hook’(OF-version dependent)
	(
	i,
	new volScalarField
		(
		IOobject
			(
			fieldName,
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
			),
		Foam::pow(epsilon1,tau)*Foam::pow(1,1-tau)*D[i]
		)
	);
}
PtrList <volScalarField>s(nc);
for (label i=0; i<s.size(); i++)
{
	// start creating fields.
	word fieldName="s"+species[i];
	Info << "Reading field " << fieldName << endl;
	
	s.set //-HM‘set or hook’(OF-version dependent)
	(
	i,
	new volScalarField
		(
		IOobject
			(
			fieldName,
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			),
		mesh
		)
	);
}
PtrList <volScalarField>S(nc);
for (label i=0; i<S.size(); i++)
{
	// start creating fields.
	word fieldName="S"+species[i];
	Info << "Reading field " << fieldName << endl;
	
	S.set //-HM‘set or hook’(OF-version dependent)
	(
	i,
	new volScalarField
		(
		IOobject
			(
			fieldName,
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
			),
		(J)/F
		)
	);
}
Info << "Assigning value from dictionary u_OH \n" << endl;
dimensionedScalar u_OH
(
	"u_OH",
	dimensionSet( 0, -1, 0, 0, 0, 0, 0),
	speciesProperties
);
Info << "Assigning value from dictionary u_K \n" << endl;
dimensionedScalar u_K
(
	"u_K",
	dimensionSet( 0, -1, 0, 0, 0, 0, 0),
	speciesProperties
);
