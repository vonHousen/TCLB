
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


#=====================================================


AddQuantity(name="Rho",	unit="kg/m3")
AddQuantity(name="T",	unit="K")
AddQuantity(name="U",	unit="m/s",	vector=T)
AddQuantity(name="E",	unit="J/m3")

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

AddSetting(	name="FluidAlfa",
			default=1,
			comment='Coefficent of heat transfer')

AddSetting(	name="FluidAlfa2",
			default=1,
			comment='Coefficent of heat transfer - another one')

AddSetting(	name="SourceTemperature",
			default=300,
			comment='Temperature of the source')

AddSetting(	name="Q",
			default=0,
			comment='Heat flux' )

AddSetting(	name="Beta",
			default=0,
			comment='Coefficent of heat expansion' )

AddSetting(	name="Tref",
			default=0,
			comment='Reference temperature for Boussinesq aproximation' )

AddSetting(	name="G_Boussinesq_X",
			default=0,
			comment='Reference temperature for Boussinesq aproximation' )

AddSetting(	name="G_Boussinesq_Y",
			default=0,
			comment='Reference temperature for Boussinesq aproximation' )

#=====================================================


AddNodeType(name="Heater",      group="TEMPBOUNDARY")
AddNodeType(name="WallSouth",	group="BOUNDARY" )
AddNodeType(name="WallNorth",	group="BOUNDARY" )
AddNodeType(name="Wall",		group="BOUNDARY" )
AddNodeType(name="BGK",			group="COLLISION" )
AddNodeType(name="DefaultAlfa",          group="TEMPBOUNDARY")
AddNodeType(name="OtherAlfa",          group="TEMPBOUNDARY")

#====================================================


AddGlobal(  name="TotalHeat",    comment='Energy', unit="J" )