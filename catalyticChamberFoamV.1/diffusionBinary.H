
Info<< "Solving Binary Diffusion and MDiji\n" << endl;

//Mass diffusion coeficient calculation

SquareMatrix<scalar> MDij(Y1.size());  //matrix multiplication i, j species
 
forAll(Y1, i)
{ 
    
        forAll(Y1, j)
        {
            if(i!=j)
            {
            
                scalar sigmaij = (sigma[i].value() + sigma[j].value())/2;
                                            
                scalar epsilonij = Foam::pow(epsilon[i].value()*epsilon[j].value(), 0.5);
                               
                scalar mWij = (1/specieData1[i].W())+(1/specieData1[j].W());
                
                scalar TN = 273.15/epsilonij;
                               
                scalar A = 1.06036; scalar B = 0.15610; scalar C = 0.19300; scalar D = 0.47635;  
                scalar E = 1.03587; scalar F = 1.52996; scalar G = 1.76474; scalar H = 3.89411; 
                
                scalar Ohmij = (A/Foam::pow(TN,B)) + (C/Foam::exp(D*TN)) + (E/Foam::exp(F*TN)) + G/Foam::exp(H*TN);

                scalar NUM = 0.0188825813*Foam::pow(273.15, 1.5)*Foam::pow(mWij, 0.5);       

                scalar DEN = 101325*Foam::pow(sigmaij, 2)*Ohmij; 

                scalar Dij0 = NUM/DEN;
   
                MDij(i, j) = 1/Dij0;
                 
            }        
        }
}

