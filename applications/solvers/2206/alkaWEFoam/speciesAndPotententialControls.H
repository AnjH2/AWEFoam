
        const dictionary& alka =
            mesh.solutionDict().subDict("AlkaWE");

        speciesCorrections =
            alka.lookupOrDefault<scalar>("speciesCorrections", speciesCorrections);

        potentialCorrections =
            alka.lookupOrDefault<scalar>("potentialCorrections", potentialCorrections);

        outerChemicalCorrections =
            alka.lookupOrDefault<scalar>("outerChemicalCorrections", outerChemicalCorrections);

