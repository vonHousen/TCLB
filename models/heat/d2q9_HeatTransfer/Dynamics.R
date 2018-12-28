
AddDensity( name="f[0]", dx= 0, dy= 0, group="f")
AddDensity( name="f[1]", dx= 1, dy= 0, group="f")
AddDensity( name="f[2]", dx= 0, dy= 1, group="f")
AddDensity( name="f[3]", dx=-1, dy= 0, group="f")
AddDensity( name="f[4]", dx= 0, dy=-1, group="f")
AddDensity( name="f[5]", dx= 1, dy= 1, group="f")
AddDensity( name="f[6]", dx=-1, dy= 1, group="f")
AddDensity( name="f[7]", dx=-1, dy=-1, group="f")
AddDensity( name="f[8]", dx= 1, dy=-1, group="f")

AddDensity( name="g[0]", dx= 0, dy= 0, group="g")
AddDensity( name="g[1]", dx= 1, dy= 0, group="g")
AddDensity( name="g[2]", dx= 0, dy= 1, group="g")
AddDensity( name="g[3]", dx=-1, dy= 0, group="g")
AddDensity( name="g[4]", dx= 0, dy=-1, group="g")
AddDensity( name="g[5]", dx= 1, dy= 1, group="g")
AddDensity( name="g[6]", dx=-1, dy= 1, group="g")
AddDensity( name="g[7]", dx=-1, dy=-1, group="g")
AddDensity( name="g[8]", dx= 1, dy=-1, group="g")

AddDensity( name="w", group="w", parameter=T )

#=====================================================


AddQuantity(name="Rho",	unit="kg/m3")
AddQuantity(name="T",	unit="K")
AddQuantity(name="U",	unit="m/s",	vector=T)
AddQuantity(name="G",	unit="m/s2", vector=T)
AddQuantity(name="E",	unit="J/m3")
AddQuantity(name="Gr")
AddQuantity(name="Pr")
AddQuantity(name="Ra")
AddQuantity(name="W")
AddQuantity( name="RhoB", adjoint=T)


#=====================================================


AddSetting(	name="Omega",
			comment='One over relaxation time' )

AddSetting(	name="Nu",
			default=0.16666666,
			Omega='1.0/(3*Nu + 0.5)',
			comment='Viscosity' )

AddSetting(	name="Density",
			default=1000,
			comment='Density of the fluid' )

AddSetting(	name="InitVelocityX",
			default=0,
			comment='Initial velocity of the fluid - x direction',
			zonal=TRUE)

AddSetting(	name="InitVelocityY",
			default=0,
			comment='Initial velocity of the fluid - y direction',
			zonal=TRUE)

AddSetting( name="G_X",
			default=0,
			comment='Acceleration x')

AddSetting( name="G_Y",
			default=0,
			comment='Acceleration y')

AddSetting(	name="InitTemperature",
			default=1,
			comment='Init temperature')

AddSetting(	name="AlfaFluid",
			default=0.1,
			comment='Coefficent of heat transfer of the solid')

AddSetting(	name="AlfaSolid",
			default=0.05,
			comment='Coefficent of heat transfer of the solid')

AddSetting(	name="SourceTemperature1",
			default=300,
			comment='Temperature of the 1st source')

AddSetting(	name="SourceTemperature2",
			default=300,
			comment='Temperature of the 2nd source')

AddSetting(	name="InitSourceTemperature",
			default=300,
			comment='Temperature of the source at initialisation')

AddSetting(	name="Q",
			default=0,
			comment='Heat flux' )

AddSetting(	name="Beta",
			default=0,
			comment='Coefficent of heat expansion' )

AddSetting(	name="Tref",
			default=0,
			comment='Reference temperat. for Boussinesq aprox' )

AddSetting(	name="Ldim",
			default=0,
			comment='Character. dimens. for dimensionless numbers' )

AddSetting(	name="Adiffusity",
			default=0,
			comment='Thermal diffusivity' )

AddSetting(	name="G_Boussinesq_X",
			default=0,
			comment='Acceleration x for Boussinesq aproximation' )

AddSetting(	name="G_Boussinesq_Y",
			default=0,
			comment='Acceleration y for Boussinesq aproximation' )

AddSetting(	name="Nu_dup",
			default=0.1666666,
			comment='Viscosity - duplicate' )

AddSetting(	name="W_init",
			default=1,
			comment='Porosity factor' )

AddSetting(	name="InletVelocity",
			default=0,
			comment='Velocity on inlet' )



#=====================================================


AddNodeType(name="Heater",          group="TEMPBOUNDARY")
AddNodeType(name="InitHeater",      group="TEMPBOUNDARY")
AddNodeType(name="WallSouth",	    group="TEMPBOUNDARY" )
AddNodeType(name="WallNorth",	    group="TEMPBOUNDARY" )
AddNodeType(name="Wall",		    group="BOUNDARY" )
AddNodeType(name="BGK",			    group="COLLISION" )
AddNodeType(name="DefaultAlfa",     group="TEMPALFA")
AddNodeType(name="OtherAlfa",       group="TEMPALFA")
AddNodeType(name="PseudoWall",      group="POROUSWALL")
AddNodeType(name="InletW",          group="INLET")
AddNodeType(name="InletGauge",      group="GAUGE")
AddNodeType(name="OutletGauge",      group="GAUGE")

#TODO alfa not in the group "TEMPALFA"

#====================================================


AddGlobal(  name="TotalHeat",    comment='Energy',      unit="J" )
AddGlobal(  name="TotalMass",    comment='Total mass',  unit="kg" )
AddGlobal(  name="MassFlowOut",  comment='Mass flow out', unit="kg/s")
AddGlobal(  name="MassFlowIn",   comment='Mass flow in', unit="kg/s")
AddGlobal(  name="Penalty",      comment='Porosity penalty function')