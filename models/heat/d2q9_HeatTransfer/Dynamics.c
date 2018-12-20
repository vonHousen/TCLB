
CudaDeviceFunction void     Init()                  //initialising function - used only once
{
	real_t  u[2]     = {InitVelocityX, InitVelocityY},
			density  = Density,
			rhoT     = density*InitTemperature,
			w_init   = 1.0;

	if(IamWall)
	{
		u[0]    = 0.0;
		u[1]    = 0.0;
	}

	if ((NodeType & NODE_TEMPBOUNDARY) == NODE_InitHeater)
		rhoT = density*InitSourceTemperature;

	if ((NodeType & NODE_POROUSWALL) == NODE_PseudoWall)
		w_init = W_init;

	w = w_init;
	SetEquilibrium_f(density, u);
	SetEquilibrium_g(rhoT, u);


}

CudaDeviceFunction void     Run()                   //main function - acts every iteration
{
	if ( (NodeType & NODE_BOUNDARY) == NODE_Wall )
		BounceBack();

	switch (NodeType & NODE_TEMPBOUNDARY)
	{
		case NODE_WallSouth:
			Heating_S();
			break;
		case NODE_WallNorth:
			Cooling_N();
			break;
		case NODE_Heater:
			real_t u[2]     = {0.0, 0.0};
			SetEquilibrium_g( Density*SourceTemperature1, u);
			break;
	}
	switch (NodeType & NODE_INLET)
	{
		case NODE_InletW:
			VelocityInlet_W();
			break;
	}

	if ((NodeType & NODE_COLLISION))
		CollisionEDM();

	w = w(0,0);


	AddToTotalHeat( getE() );
	AddToTotalMass( getRho() );
}

CudaDeviceFunction float2   Color()                 //does nothing - no CUDA
{
	float2 ret;
	ret.x = (getT()-Tref)/Tref;
	ret.y = 1;

	return ret;
}

CudaDeviceFunction real_t   getRho()                //gets density at the current node.
{
    return ( f[8]+f[7]+f[6]+f[5]+f[4]+f[3]+f[2]+f[1]+f[0] );
}

CudaDeviceFunction real_t   getT()                  //gets temperature at the current node.
{
    return ( g[8]+g[7]+g[6]+g[5]+g[4]+g[3]+g[2]+g[1]+g[0] )/getRho();
}

CudaDeviceFunction real_t   getGr()                 //gets value of Grashof number at the current node.
{
	real_t  g = sqrt(G_Boussinesq_X*G_Boussinesq_X + G_Boussinesq_Y*G_Boussinesq_Y),
			nu = Nu_dup;

	return ( (1/Tref) * fabs(getT() - Tref) * (g*Ldim*Ldim*Ldim)/(nu*nu) );
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
	return ( g[8]+g[7]+g[6]+g[5]+g[4]+g[3]+g[2]+g[1]+g[0] );
}

CudaDeviceFunction real_t   getW()                  //gets Porosity factor at the current node.
{
	return ( w(0,0) );
}

CudaDeviceFunction vector_t getG()                  //gets acceleration vector at the current node
{
	vector_t    g;

	real_t 		density = getRho(),
	            t       = getT();

	g.x     = G_X - (1/Tref) * G_Boussinesq_X * (t - Tref);
	g.y     = G_Y - (1/Tref) * G_Boussinesq_Y * (t - Tref);

	return g;
}

CudaDeviceFunction vector_t getU()                  //gets velocity vector at the current node
{
	real_t 		density = getRho();
	vector_t 	u,
				g = getG();

	if(not IamWall)
	{
		// pu' = pu + G/2
		u.x = (( f[8] - f[7] - f[6] + f[5] - f[3] + f[1]) / density + g.x * 0.5);
		u.y = ((-f[8] - f[7] + f[6] + f[5] - f[4] + f[2]) / density + g.y * 0.5);
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

CudaDeviceFunction void     Heating_S()             //boundary Zou He like condition for heating on south wall
{

	real_t	density = getRho();

	g[2] += Q*density * 2.0/3.0;
	g[5] += Q*density * 1.0/6.0;
	g[6] += Q*density * 1.0/6.0;

}

CudaDeviceFunction void     Cooling_N()             //boundary Zou He like condition for cooling on north wall
{

	real_t	density = getRho();

	g[4] -= Q*density * 2.0/3.0;
	g[7] -= Q*density * 1.0/6.0;
	g[8] -= Q*density * 1.0/6.0;

}

CudaDeviceFunction void     VelocityInlet_W()       //boundary Zou He like condition for velocity inlet on western wall
{
	//not correct
	/*
	real_t	density = getRho();

	f[1] = InletVelocity*density * 2.0/3.0;
	f[5] = InletVelocity*density * 1.0/6.0;
	f[8] = InletVelocity*density * 1.0/6.0;
	f[0] = 0.0;
	f[2] = 0.0;
	f[3] = 0.0;
	f[4] = 0.0;
	f[6] = 0.0;
	f[7] = 0.0;
	*/
}


//======================

													//calculates the equilibrium distribution of field
CudaDeviceFunction void     SetEquilibrium_f(const real_t density, const real_t *u)
{
	//  relaxation factor
	real_t  S[9] =  { 4.0 /  9.0,
					  1.0 /  9.0,
					  1.0 /  9.0,
					  1.0 /  9.0,
					  1.0 /  9.0,
					  1.0 / 36.0,
					  1.0 / 36.0,
					  1.0 / 36.0,
					  1.0 / 36.0  };

	//d2q9 - 9 lattice directions
	real_t  c[9][2] = {  { 0.0, 0.0},
                         { 1.0, 0.0},
                         { 0.0, 1.0},
                         {-1.0, 0.0},
                         { 0.0,-1.0},
                         { 1.0, 1.0},
                         {-1.0, 1.0},
                         {-1.0,-1.0},
                         { 1.0,-1.0}  };

	real_t cu_temp, u2_temp;      //cu_temp = c*u, u2_temp = u^2
	real_t c2_s = 1.0 / 3.0;      //c2_s = c_s^2 <==> lattice speed of sound




	for(int i=0; i<9; i++)
	{
		cu_temp = c[i][0]*u[0] + c[i][1]*u[1];
		u2_temp = u[0]*u[0] + u[1]*u[1];

		//f_eq = ...
		f[i] = S[i] * density * ( 1 + cu_temp/c2_s + (cu_temp*cu_temp)/(2.0*c2_s*c2_s) - u2_temp/(2.0*c2_s) );
	}
}

													//calculates the equilibrium distribution of g - field
CudaDeviceFunction void     SetEquilibrium_g(const real_t rhoT, const real_t u[2])
{
	//  relaxation factor
	real_t  S[9] =  { 4.0 /  9.0,
	                  1.0 /  9.0,
	                  1.0 /  9.0,
	                  1.0 /  9.0,
	                  1.0 /  9.0,
	                  1.0 / 36.0,
	                  1.0 / 36.0,
	                  1.0 / 36.0,
	                  1.0 / 36.0  };

	//d2q9 - 9 lattice directions
	real_t  c[9][2] = {  { 0.0, 0.0},
	                     { 1.0, 0.0},
	                     { 0.0, 1.0},
	                     {-1.0, 0.0},
	                     { 0.0,-1.0},
	                     { 1.0, 1.0},
	                     {-1.0, 1.0},
	                     {-1.0,-1.0},
	                     { 1.0,-1.0}  };

	real_t cu_temp, u2_temp;      //cu_temp = c*u, u2_temp = u^2
	real_t c2_s = 1.0 / 3.0;      //c2_s = c_s^2 <==> lattice speed of sound




	for(int i=0; i<9; i++)
	{
		cu_temp = c[i][0]*u[0] + c[i][1]*u[1];
		u2_temp = u[0]*u[0] + u[1]*u[1];

		//f_eq = ...
		g[i] = S[i] * rhoT * ( 1 + cu_temp/c2_s + (cu_temp*cu_temp)/(2.0*c2_s*c2_s) - u2_temp/(2.0*c2_s) );
	}
}


CudaDeviceFunction vector_t   G(const real_t w, const real_t* u)       //function for calculating Darcy's acceler.
{
	real_t      w_temp = w;
	vector_t    u_temp;

	u_temp.x = u[0];
	u_temp.y = u[1];

	if(w > 1.0)
		w_temp = 1.0;
	else if(w<0.0)
		w_temp = 0.0;

	w_temp = -pow(1.0-w_temp, 0.01);
	u_temp.x *= w_temp;
	u_temp.y *= w_temp;

	return u_temp;
}

CudaDeviceFunction real_t   AlfaT(const real_t w)   //function for interpolating AlfaT
{
	real_t alfa;

	if( w >= 1.0 )
		alfa = AlfaFluid;
	else if( w<=0.0 )
		alfa = AlfaSolid;
	else
		alfa = AlfaSolid*(1.0-w) + AlfaFluid*w;


	return alfa;
}

CudaDeviceFunction void     CollisionEDM()          //physics of the collision (based on Exact Difference Method)
{
	//before collision:
	real_t      density                 = getRho(),
                u_before_collision[2]   = { (( f[8]-f[7]-f[6]+f[5]-f[3]+f[1] )/density ),
	                                        ((-f[8]-f[7]+f[6]+f[5]-f[4]+f[2] )/density )  },
				f_before_collision[9]   = { f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8] },
				u_temp[2],
				u[2];
	vector_t    acceleration            = getG(),
				Darcy;


	SetEquilibrium_f(density, u_before_collision);

	//after collision:
	real_t      f_unforced[9]           = { f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8] };
	u_temp[0] = u_before_collision[0] + acceleration.x/Omega ;
	u_temp[1] = u_before_collision[1] + acceleration.y/Omega ;

	Darcy           = G( w(0,0), u_temp );
	acceleration.x +=  + Darcy.x;
	acceleration.y +=  + Darcy.y;
	u[0] = u_temp[0] + acceleration.x/Omega ;
	u[1] = u_temp[1] + acceleration.y/Omega ;

	SetEquilibrium_f(getRho(), u);

	for(int i=0; i<9; i++) {
		f[i] = (f_before_collision[i] - f_unforced[i]) * (1 - Omega) + f[i];
	}

	//=========== HEAT ===========
	//saving memory by using f-variables

	real_t  omega_T         = 1.0/(3* AlfaT( w(0,0) ) + 0.5),
			g_before_collision[9] = { g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8] };


	SetEquilibrium_g(density, u_before_collision);

	//after collision
	real_t      g_unforced[9]           = { g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8] };

	SetEquilibrium_g(getRho(), u);

	for(int i=0; i<9; i++) {
		g[i] = (g_before_collision[i] - g_unforced[i]) * (1 - omega_T) + g[i];
	}

	//

}

