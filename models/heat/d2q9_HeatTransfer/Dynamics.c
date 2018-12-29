
CudaDeviceFunction void     Init()                  //initialising function - used only once
{
	real_t  u[2]     = {InitVelocityX, InitVelocityY},
			density  = Density,
			rhoT     = density*InitTemperature,
			w_init   = 1.0;

	if((NodeType & NODE_BOUNDARY) == NODE_Wall)
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
	real_t      u[2];


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
			u[0] = 0.0;
			u[1] = 0.0;
			SetEquilibrium_g( getRho()*SourceTemperature1, u);
			break;
	}
	switch (NodeType & NODE_INLET)
	{
		case NODE_InletW:
			VelocityInlet_W();
			break;
	}

	if ((NodeType & NODE_COLLISION) == NODE_BGK )
		CollisionEDM();



	switch (NodeType & NODE_GAUGE)
	{
		case NODE_InletGauge:
			AddToMassFlowIn(  u_x() * getRho()/1.0 );
			break;
		case NODE_OutletGauge:
			AddToMassFlowOut( u_x() * getRho()/1.0 );
			break;
	}

	AddToTotalHeat( getE() );
	AddToTotalMass( getRho() );
	AddToPenalty( w*(1-w) );
}

CudaDeviceFunction float2   Color()                 //does nothing - no CUDA
{
	float2 ret;
	ret.x = (float)((getT()-Tref)/Tref);
	ret.y = (float)1;

	return ret;
}

CudaDeviceFunction real_t   getRho()                //gets density at the current node.
{
    return ( f8+f7+f6+f5+f4+f3+f2+f1+f0 );
}

CudaDeviceFunction real_t   getT()                  //gets temperature at the current node.
{
    return ( g8+g7+g6+g5+g4+g3+g2+g1+g0 )/getRho();
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
	return ( getT()*getRho() );
}

CudaDeviceFunction real_t   getW()                  //gets Porosity factor at the current node.
{
	return ( w );
}

CudaDeviceFunction vector_t getG()                  //gets acceleration vector at the current node
{
	vector_t    g;

	g.x     = acceleration_x();
	g.y     = acceleration_y();
	g.z     = 0.0;

	return g;
}

CudaDeviceFunction vector_t getU()                  //gets velocity vector at the current node
{
	vector_t 	u;

	u.x = u_x();
	u.y = u_y();
	u.z = 0.0;

	return u;
}


//======================


CudaDeviceFunction void     BounceBack()            //bouncing back on walls
{
    real_t temp, uf;

	temp	= f3;
	f3	    = f1;
	f1	    = temp;
	
	temp 	= f4;
	f4	    = f2;
	f2	    = temp;

	temp 	= f7;
	f7	    = f5;
	f5	    = temp;

	temp 	= f8;
	f8	    = f6;
	f6	    = temp;

	// 0 1 2 3 4 5 6 7 8
	// 1 5 2 6 3 7 4 8 0

	//constant temperature boundary condition
	//SetEquilibrium_g(getRho()*SourceTemperature, u);

	//adiabatic wall
	uf 		= g3;
	g3  	= g1;
	g1  	= uf;

	uf 		= g4;
	g4	    = g2;
	g2	    = uf;

	uf 		= g7;
	g7	    = g5;
	g5  	= uf;

	uf 		= g8;
	g8	    = g6;
	g6	    = uf;
}

CudaDeviceFunction void     Heating_S()             //boundary Zou He like condition for heating on south wall
{

	real_t	density = getRho();

	g2 += Q*density * 2.0/3.0;
	g5 += Q*density * 1.0/6.0;
	g6 += Q*density * 1.0/6.0;

}

CudaDeviceFunction void     Cooling_N()             //boundary Zou He like condition for cooling on north wall
{

	real_t	density = getRho();

	g4 -= Q*density * 2.0/3.0;
	g7 -= Q*density * 1.0/6.0;
	g8 -= Q*density * 1.0/6.0;

}

CudaDeviceFunction void     VelocityInlet_W()       //boundary Zou He like condition for velocity inlet on western wall
{
	//not correct
}


//======================

CudaDeviceFunction void     f_assign(real_t f_temp[9])   //assigns f - field based on f_temp[9]
{
	f0 = f_temp[0];
	f1 = f_temp[1];
	f2 = f_temp[2];
	f3 = f_temp[3];
	f4 = f_temp[4];
	f5 = f_temp[5];
	f6 = f_temp[6];
	f7 = f_temp[7];
	f8 = f_temp[8];
}


CudaDeviceFunction void     g_assign(real_t g[9])   //assigns g - field based on g[9]
{
	g0 = g[0];
	g1 = g[1];
	g2 = g[2];
	g3 = g[3];
	g4 = g[4];
	g5 = g[5];
	g6 = g[6];
	g7 = g[7];
	g8 = g[8];
}
													//calculates the equilibrium distribution of field
CudaDeviceFunction void     SetEquilibrium_f(real_t density, real_t u[2])
{
	int i=0;
	real_t f_temp[9];

	//  relaxation factor
	real_t  S[9];
			S[0] = (4.0 /  9.0);
			S[1] = (1.0 /  9.0);
			S[2] = (1.0 /  9.0);
			S[3] = (1.0 /  9.0);
			S[4] = (1.0 /  9.0);
			S[5] = (1.0 / 36.0);
			S[6] = (1.0 / 36.0);
			S[7] = (1.0 / 36.0);
			S[8] = (1.0 / 36.0);

	//d2q9 - 9 lattice directions
	real_t  c[9][2];
			c[0][0] =  0.0;
			c[0][1] =  0.0;
			c[1][0] =  1.0;
			c[1][1] =  0.0;
			c[2][0] =  0.0;
			c[2][1] =  1.0;
			c[3][0] = -1.0;
			c[3][1] =  0.0;
			c[4][0] =  0.0;
			c[4][1] = -1.0;
			c[5][0] =  1.0;
			c[5][1] =  1.0;
			c[6][0] = -1.0;
			c[6][1] =  1.0;
			c[7][0] = -1.0;
			c[7][1] = -1.0;
			c[8][0] =  1.0;
			c[8][1] = -1.0;

	real_t cu_temp, u2_temp;      //cu_temp = c*u, u2_temp = u^2
	real_t c2_s = 1.0 / 3.0;      //c2_s = c_s^2 <==> lattice speed of sound




	for(i=0; i<9; i++)
	{
		cu_temp = c[i][0]*u[0] + c[i][1]*u[1];
		u2_temp = u[0]*u[0] + u[1]*u[1];

		//f_eq = ...
		f_temp[i] = S[i] * density * ( 1 + cu_temp/c2_s + (cu_temp*cu_temp)/(2.0*c2_s*c2_s) - u2_temp/(2.0*c2_s) );
	}

	f_assign(f_temp);
}

													//calculates the equilibrium distribution of g - field
CudaDeviceFunction void     SetEquilibrium_g(real_t rhoT, real_t u[2])
{
	int i=0;
	real_t g_temp[9];

	//  relaxation factor
	real_t  S[9];
			S[0] = (4.0 /  9.0);
			S[1] = (1.0 /  9.0);
			S[2] = (1.0 /  9.0);
			S[3] = (1.0 /  9.0);
			S[4] = (1.0 /  9.0);
			S[5] = (1.0 / 36.0);
			S[6] = (1.0 / 36.0);
			S[7] = (1.0 / 36.0);
			S[8] = (1.0 / 36.0);

	//d2q9 - 9 lattice directions
	real_t  c[9][2];
			c[0][0] =  0.0;
			c[0][1] =  0.0;
			c[1][0] =  1.0;
			c[1][1] =  0.0;
			c[2][0] =  0.0;
			c[2][1] =  1.0;
			c[3][0] = -1.0;
			c[3][1] =  0.0;
			c[4][0] =  0.0;
			c[4][1] = -1.0;
			c[5][0] =  1.0;
			c[5][1] =  1.0;
			c[6][0] = -1.0;
			c[6][1] =  1.0;
			c[7][0] = -1.0;
			c[7][1] = -1.0;
			c[8][0] =  1.0;
			c[8][1] = -1.0;

	real_t cu_temp, u2_temp;      //cu_temp = c*u, u2_temp = u^2
	real_t c2_s = 1.0 / 3.0;      //c2_s = c_s^2 <==> lattice speed of sound




	for(i=0; i<9; i++)
	{
		cu_temp = c[i][0]*u[0] + c[i][1]*u[1];
		u2_temp = u[0]*u[0] + u[1]*u[1];

		//f_eq = ...
		g_temp[i] = S[i] * rhoT * ( 1 + cu_temp/c2_s + (cu_temp*cu_temp)/(2.0*c2_s*c2_s) - u2_temp/(2.0*c2_s) );
	}

	g_assign(g_temp);
}

													//returns Darcy's acceleration on one direction
CudaDeviceFunction real_t   G_darcy(real_t w, real_t u)
{
	real_t      w_temp  = w,
				g       = u;

	if(w > 1.0)
		w_temp = 1.0;
	else if(w<0.0)
		w_temp = 0.0;

	w_temp = (-1.0)*exp(1000.0*(w_temp-1.0));
	g *= w_temp;

	return g;
}

CudaDeviceFunction real_t   acceleration_x()    //returns acceleration_x
{
	return ( G_X - (1/Tref) * G_Boussinesq_X * (getT() - Tref) );
}


CudaDeviceFunction real_t   acceleration_y()    //returns acceleration_y
{
	return ( G_Y - (1/Tref) * G_Boussinesq_Y * (getT() - Tref) );
}


CudaDeviceFunction real_t   u_x()               //returns velocity_x
{
	real_t 		u;

	if( ((NodeType & NODE_BOUNDARY) == NODE_Wall) )
	{
		u = 0.0;
	}
	else
	{
		u = (( f8 - f7 - f6 + f5 - f3 + f1) / getRho() + acceleration_x() * 0.5);
	}

	return u;
}


CudaDeviceFunction real_t   u_y()               //returns velocity_y
{
	real_t 		u;

	if( ((NodeType & NODE_BOUNDARY) == NODE_Wall) )
	{
		u = 0.0;
	}
	else
	{
		u = ((-f8 - f7 + f6 + f5 - f4 + f2) / getRho() + acceleration_y() * 0.5);
	}

	return u;
}


CudaDeviceFunction real_t   AlfaT(real_t w)     //function for interpolating AlfaT
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


CudaDeviceFunction void     CollisionEDM()      //physics of the collision (based on Exact Difference Method)
{
	//before collision:
	int i=0;
	real_t      density                 = getRho();
	real_t      u_before_collision[2],
				f_before_collision[9],
				u_temp[2],
				u[2],
				f_temp[9],
				acceleration[2];
				acceleration[0] = acceleration_x();
				acceleration[1] = acceleration_y();
		f_temp[0] = f0;
		f_temp[1] = f1;
		f_temp[2] = f2;
		f_temp[3] = f3;
		f_temp[4] = f4;
		f_temp[5] = f5;
		f_temp[6] = f6;
		f_temp[7] = f7;
		f_temp[8] = f8;	
	for(i=0; i<9; i++)      f_before_collision[i] = f_temp[i];

	u_before_collision[0] = ( f8-f7-f6+f5-f3+f1 )/density;
	u_before_collision[1] = (-f8-f7+f6+f5-f4+f2 )/density;

	SetEquilibrium_f(density, u_before_collision);

	//after collision:
	real_t      f_unforced[9];

		f_temp[0] = f0;
		f_temp[1] = f1;
		f_temp[2] = f2;
		f_temp[3] = f3;
		f_temp[4] = f4;
		f_temp[5] = f5;
		f_temp[6] = f6;
		f_temp[7] = f7;
		f_temp[8] = f8;
	for(i=0; i<9; i++)      f_unforced[i]   = f_temp[i];
	for(i=0; i<2; i++)      u_temp[i]       = u_before_collision[i] + acceleration[i]/Omega ;


	acceleration[0]        +=  G_darcy(w, u_temp[0]);
	acceleration[1]        +=  G_darcy(w, u_temp[1]);
	u[0]                    = u_temp[0] + acceleration[0]/Omega ;
	u[1]                    = u_temp[1] + acceleration[1]/Omega ;

	SetEquilibrium_f(getRho(), u);

		f_temp[0] = f0;
		f_temp[1] = f1;
		f_temp[2] = f2;
		f_temp[3] = f3;
		f_temp[4] = f4;
		f_temp[5] = f5;
		f_temp[6] = f6;
		f_temp[7] = f7;
		f_temp[8] = f8;
	for(i=0; i<9; i++) {
		f_temp[i] = (f_before_collision[i] - f_unforced[i]) * (1 - Omega) + f_temp[i];
	}
	
	f_assign(f_temp);

	//=========== HEAT ===========
	//saving memory by using f-variables

	real_t  omega_T                 = 1.0/(3* AlfaT( w ) + 0.5),
			rhoT                    = g0 + g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8,
			g_before_collision[9],
			g_temp[9];

		g_temp[0] = g0;
		g_temp[1] = g1;
		g_temp[2] = g2;
		g_temp[3] = g3;
		g_temp[4] = g4;
		g_temp[5] = g5;
		g_temp[6] = g6;
		g_temp[7] = g7;
		g_temp[8] = g8;
	for(i=0; i<9; i++)      g_before_collision[i] = g_temp[i];
	SetEquilibrium_g(rhoT, u_before_collision);

	//after collision
		g_temp[0] = g0;
		g_temp[1] = g1;
		g_temp[2] = g2;
		g_temp[3] = g3;
		g_temp[4] = g4;
		g_temp[5] = g5;
		g_temp[6] = g6;
		g_temp[7] = g7;
		g_temp[8] = g8;
	real_t                      g_unforced[9];
	for(i=0; i<9; i++)      g_unforced[i]   = g_temp[i];
	rhoT = g0 + g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8;

	SetEquilibrium_g(rhoT, u);

	for(i=0; i<9; i++) {
		g_temp[i] = (g_before_collision[i] - g_unforced[i]) * (1 - omega_T) + g_temp[i];
	}

	g_assign(g_temp);

}

