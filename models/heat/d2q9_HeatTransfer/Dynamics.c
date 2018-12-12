
CudaDeviceFunction void     Init()                  //initialising function - used only once
{
	real_t u[2]     = {InitVelocityX, InitVelocityY};
	real_t density  = Density;
	real_t rhoT     = density*InitTemperature;

	if(IamWall)
	{
		u[0]    = 0.0;
		u[1]    = 0.0;
	}

	SetEquilibrium_f(density, u);

	if ((NodeType & NODE_TEMPBOUNDARY) == NODE_Heater)
	{
		rhoT = density*SourceTemperature;
	}

	SetEquilibrium_g(rhoT, u);

}

CudaDeviceFunction void     Run()                   //main function - acts every iteration
{
	switch (NodeType & NODE_BOUNDARY)
	{
		case NODE_WallSouth:
			Heating_S();
			break;
		case NODE_WallNorth:
			Cooling_N();
			break;
		case NODE_Wall:
			BounceBack();
			break;
	}
	if ((NodeType & NODE_COLLISION))
	{
		CollisionEDM();
	}

	AddToTotalHeat( getE() );

}

CudaDeviceFunction float2   Color()                 //does nothing - no CUDA
{

	float2 ret{0.0, 0.0};

	return ret;
}

CudaDeviceFunction real_t   getRho()                //gets density at the current node.
{
    return f[8]+f[7]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0];
}

CudaDeviceFunction real_t   getT()                  //gets temperature at the current node.
{
    return (g[8]+g[7]+g[6]+g[5]+g[4]+g[3]+g[2]+g[1]+g[0])/getRho();
}

CudaDeviceFunction real_t   getGr()                 //gets value of Grashof number at the current node.
{
	real_t  g = sqrt(G_Boussinesq_X*G_Boussinesq_X + G_Boussinesq_Y*G_Boussinesq_Y),
			nu = Nu_dup;

	return ( (1/Tref) * (getT() - Tref) * (g*Ldim*Ldim*Ldim)/(nu*nu) );
}

CudaDeviceFunction real_t   getPr()                 //gets value of Prandtl number at the current node.
{
	real_t nu = Nu_dup;

	return ( nu / Adiffusity );
}

CudaDeviceFunction real_t   getRa()                 //gets value of Rayleigh number at the current node.
{
	return ( getGr()*getPr() );
}

CudaDeviceFunction real_t   getE()                  //gets Energy at the current node.
{
	return g[8]+g[7]+g[6]+g[5]+g[4]+g[3]+g[2]+g[1]+g[0];
}

CudaDeviceFunction vector_t getU()                  //gets velocity vector at the current node
{
	real_t 		density = getRho(),
				g_x     = G_X,
				g_y     = G_Y,
				t       = getT();
	vector_t 	u;

	g_x += - (1/Tref) * G_Boussinesq_X * (t - Tref);
	g_y += - (1/Tref) * G_Boussinesq_Y * (t - Tref);

	if(not IamWall)
	{
		// pu' = pu + G/2
		u.x = (( f[8] - f[7] - f[6] + f[5] - f[3] + f[1]) / density + g_x * 0.5);
		u.y = ((-f[8] - f[7] + f[6] + f[5] - f[4] + f[2]) / density + g_y * 0.5);
		u.z = 0.0;
	}
	else
	{
		u.x = u.y = u.z = 0.0;
	}

	return u;
}


//======================


CudaDeviceFunction void     BounceBack()            //bouncing back f - field
{
    real_t temp, u[2]={0.0, 0.0}, uf;

	temp	= f[3];
	f[3]	= f[1];
	f[1]	= temp;
	
	temp 	= f[4];
	f[4]	= f[2];
	f[2]	= temp;

	temp 	= f[7];
	f[7]	= f[5];
	f[5]	= temp;

	temp 	= f[8];
	f[8]	= f[6];
	f[6]	= temp;

	// 0 1 2 3 4 5 6 7 8
	// 1 5 2 6 3 7 4 8 0

	//constant temperature boundary condition
	//SetEquilibrium_g(getRho()*SourceTemperature, u);

	//adiabatic wall
	uf 		= g[3];
	g[3]	= g[1];
	g[1]	= uf;

	uf 		= g[4];
	g[4]	= g[2];
	g[2]	= uf;

	uf 		= g[7];
	g[7]	= g[5];
	g[5]	= uf;

	uf 		= g[8];
	g[8]	= g[6];
	g[6]	= uf;
}

CudaDeviceFunction void     AdiabaticSouth()        //boundary condition for south wall - TODO: complete south wall
{
	/*real_t temp;

	temp	= f[3];
	f[3]	= f[1];
	f[1]	= temp;

	temp 	= f[4];
	f[4]	= f[2];
	f[2]	= temp;

	temp 	= f[7];
	f[7]	= f[5];
	f[5]	= temp;

	temp 	= f[8];
	f[8]	= f[6];
	f[6]	= temp;


	temp 	= g[3];
	g[3]	= g[1];
	g[1]	= temp;

	temp 	= g[4];
	g[4]	= 0;
	g[2]	= temp;

	temp 	= g[7];
	g[7]	= 0;
	g[5]	= temp;

	temp 	= g[8];
	g[8]	= 0;
	g[6]	= temp;*/
}

CudaDeviceFunction void     AdiabaticNorth()        //boundary condition for north wall - TODO: complete north wall
{
	/*real_t temp;

	temp	= f[3];
	f[3]	= f[1];
	f[1]	= temp;

	temp 	= f[4];
	f[4]	= f[2];
	f[2]	= temp;

	temp 	= f[7];
	f[7]	= f[5];
	f[5]	= temp;

	temp 	= f[8];
	f[8]	= f[6];
	f[6]	= temp;


	temp 	= g[3];
	g[3]	= g[1];
	g[1]	= temp;

	temp 	= g[4];
	g[4]	= g[2];
	g[2]	= 0;

	temp 	= g[7];
	g[7]	= g[5];
	g[5]	= 0;

	temp 	= g[8];
	g[8]	= g[6];
	g[6]	= 0;*/
}

CudaDeviceFunction void     Heating_S()             //boundary Zou He like condition for heating on south wall
{

	real_t	uf 		= g[3],
			density = getRho();

	//first: adiabatic wall
	g[3]	= g[1];
	g[1]	= uf;

	uf 		= g[4];
	g[4]	= g[2];
	g[2]	= uf;

	uf 		= g[7];
	g[7]	= g[5];
	g[5]	= uf;

	uf 		= g[8];
	g[8]	= g[6];
	g[6]	= uf;

	//then
	g[2] += Q*density * 2.0/3.0;
	g[5] += Q*density * 1.0/6.0;
	g[6] += Q*density * 1.0/6.0;


}

CudaDeviceFunction void     Cooling_N()             //boundary Zou He like condition for cooling on north wall
{

	real_t	uf 		= g[3],
			density = getRho();

	//first: adiabatic wall
	g[3]	= g[1];
	g[1]	= uf;

	uf 		= g[4];
	g[4]	= g[2];
	g[2]	= uf;

	uf 		= g[7];
	g[7]	= g[5];
	g[5]	= uf;

	uf 		= g[8];
	g[8]	= g[6];
	g[6]	= uf;

	//then

	g[4] -= Q*density * 2.0/3.0;
	g[7] -= Q*density * 1.0/6.0;
	g[8] -= Q*density * 1.0/6.0;


}


//======================

													//calculates the equilibrium distribution of f - field
CudaDeviceFunction void     SetEquilibrium_f(const real_t density, const real_t *u)
{
	//  relaxation factor
	static constexpr real_t  S[9] = { 4.0 /  9.0,
									  1.0 /  9.0,
									  1.0 /  9.0,
									  1.0 /  9.0,
									  1.0 /  9.0,
									  1.0 / 36.0,
									  1.0 / 36.0,
									  1.0 / 36.0,
									  1.0 / 36.0  };

	//d2q9 - 9 lattice directions
	static constexpr real_t  c[9][2] = { { 0.0, 0.0},
	                                     { 1.0, 0.0},
	                                     { 0.0, 1.0},
	                                     {-1.0, 0.0},
	                                     { 0.0,-1.0},
	                                     { 1.0, 1.0},
	                                     {-1.0, 1.0},
	                                     {-1.0,-1.0},
	                                     { 1.0,-1.0}  };

	real_t cu_temp, u2_temp;                        //cu_temp = c*u, u2_temp = u^2
	constexpr real_t c2_s = 1.0 / 3.0;      //c2_s = c_s^2 <==> lattice speed of sound




	for(int i=0; i<9; i++)
	{
		cu_temp = c[i][0]*u[0] + c[i][1]*u[1];
		u2_temp = u[0]*u[0] + u[1]*u[1];

		//f_eq = ...
		f[i] = S[i] * density * ( 1 + cu_temp/c2_s + (cu_temp*cu_temp)/(2.0*c2_s*c2_s) - u2_temp/(2.0*c2_s) );
	}

	//below the same effect but from tutorial
	/*
	real_t d=density;
	f[0] = ( 2. + ( -u[1]*u[1] - u[0]*u[0] )*3. )*d*2./9.;
	f[1] = ( 2. + ( -u[1]*u[1] + ( 1 + u[0] )*u[0]*2. )*3. )*d/18.;
	f[2] = ( 2. + ( -u[0]*u[0] + ( 1 + u[1] )*u[1]*2. )*3. )*d/18.;
	f[3] = ( 2. + ( -u[1]*u[1] + ( -1 + u[0] )*u[0]*2. )*3. )*d/18.;
	f[4] = ( 2. + ( -u[0]*u[0] + ( -1 + u[1] )*u[1]*2. )*3. )*d/18.;
	f[5] = ( 1. + ( ( 1 + u[1] )*u[1] + ( 1 + u[0] + u[1]*3. )*u[0] )*3. )*d/36.;
	f[6] = ( 1. + ( ( 1 + u[1] )*u[1] + ( -1 + u[0] - u[1]*3. )*u[0] )*3. )*d/36.;
	f[7] = ( 1. + ( ( -1 + u[1] )*u[1] + ( -1 + u[0] + u[1]*3. )*u[0] )*3. )*d/36.;
	f[8] = ( 1. + ( ( -1 + u[1] )*u[1] + ( 1 + u[0] - u[1]*3. )*u[0] )*3. )*d/36.;
	*/
}

													//calculates the equilibrium distribution of g - field
CudaDeviceFunction void     SetEquilibrium_g(const real_t rhoT, const real_t u[2])
{
	//  relaxation factor
	static constexpr real_t  S[9] = { 4.0 /  9.0,
	                                  1.0 /  9.0,
	                                  1.0 /  9.0,
	                                  1.0 /  9.0,
	                                  1.0 /  9.0,
	                                  1.0 / 36.0,
	                                  1.0 / 36.0,
	                                  1.0 / 36.0,
	                                  1.0 / 36.0  };

	//d2q9 - 9 lattice directions
	static constexpr real_t  c[9][2] = { { 0.0, 0.0},
	                                     { 1.0, 0.0},
	                                     { 0.0, 1.0},
	                                     {-1.0, 0.0},
	                                     { 0.0,-1.0},
	                                     { 1.0, 1.0},
	                                     {-1.0, 1.0},
	                                     {-1.0,-1.0},
	                                     { 1.0,-1.0}  };

	real_t cu_temp, u2_temp;                        //cu_temp = c*u, u2_temp = u^2
	constexpr real_t c2_s = 1.0 / 3.0;      //c2_s = c_s^2 <==> lattice speed of sound




	for(int i=0; i<9; i++)
	{
		cu_temp = c[i][0]*u[0] + c[i][1]*u[1];
		u2_temp = u[0]*u[0] + u[1]*u[1];

		//f_eq = ...
		g[i] = S[i] * rhoT * ( 1 + cu_temp/c2_s + (cu_temp*cu_temp)/(2.0*c2_s*c2_s) - u2_temp/(2.0*c2_s) );
	}
}


CudaDeviceFunction void     CollisionEDM()          //physics of the collision (based on Exact Difference Method)
{
	real_t      density                 = getRho(),
                u[2]                    = { (( f[8]-f[7]-f[6]+f[5]-f[3]+f[1] )/density ),
	                                        ((-f[8]-f[7]+f[6]+f[5]-f[4]+f[2] )/density )  },
				f_before_collision[9]   = { f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8] };

	SetEquilibrium_f(density, u);
	real_t      f_unforced[9]           = { f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8] };

	u[0] = (( f[8]-f[7]-f[6]+f[5]-f[3]+f[1] )/density + G_X/Omega );
	u[1] = ((-f[8]-f[7]+f[6]+f[5]-f[4]+f[2] )/density + G_Y/Omega );
	SetEquilibrium_f(density, u);

	for(int i=0; i<9; i++) {
		f[i] = (f_before_collision[i] - f_unforced[i]) * (1 - Omega) + f[i];
	}

	//=========== HEAT ===========
	//saving memory by using f-variables

	real_t  fluidAlfa;
	if ((NodeType & NODE_TEMPBOUNDARY)==NODE_DefaultAlfa)
	{
		fluidAlfa = FluidAlfa;
	}
	else if ((NodeType & NODE_TEMPBOUNDARY)==NODE_OtherAlfa)
	{
		fluidAlfa = FluidAlfa2;
	}


	real_t  omegaT          = 1.0/(3*fluidAlfa + 0.5),
			rhoT            = density*getT();

	u[0]                    = getU().x ;
	u[1]                    = getU().y ;
	f_before_collision[0] = g[0];
	f_before_collision[1] = g[1];
	f_before_collision[2] = g[2];
	f_before_collision[3] = g[3];
	f_before_collision[4] = g[4];
	f_before_collision[5] = g[5];
	f_before_collision[6] = g[6];
	f_before_collision[7] = g[7];
	f_before_collision[8] = g[8];

	SetEquilibrium_g(rhoT, u);
	f_unforced[0] = g[0];
	f_unforced[1] = g[1];
	f_unforced[2] = g[2];
	f_unforced[3] = g[3];
	f_unforced[4] = g[4];
	f_unforced[5] = g[5];
	f_unforced[6] = g[6];
	f_unforced[7] = g[7];
	f_unforced[8] = g[8];

	for(int i=0; i<9; i++) {
		g[i] = (f_before_collision[i] - f_unforced[i]) * (1 - omegaT) + g[i];
	}

}

