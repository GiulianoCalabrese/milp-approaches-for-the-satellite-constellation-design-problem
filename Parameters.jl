using SatelliteToolbox

#################### LIST of PARAMETERS ######################

# RAYON de la Terre
	global RAYON = 6378136.3 # m
	
# Circumeference of a sphere (Earth)
	global CircumferenceEarth = RAYON*2*pi # Km

# hai [1261.9898537417057, 893.6996722009753, 566.8052002788218] come differenti altitudini per 24h di simulazione
# ottieni questo con apertura 40 gradi come lato in KM del quadrato del satellite : 930.1114958405926 km   655.7710767836845 km    414.35089323796853 km

# Temps total de simulation
	global TEMPS_SIMU = 24.0*3600 # sec
	# sideral time in seconds 23*3600 + 56*60 + 4

# Temps max entre 2 revisite d'un point au sol par un satellite
	global DELTA_T = 12.0*3600 # sec

# Pas de temps de la simulation
	global dt = 180.0# sec

# Standard gravitational parameter of the central body
	const μ = 3.986004418e14 # [m^3/s^2]
	
# Nombre de pas de temps dans la simulation
	global NUM_TIME = convert(Int,round(TEMPS_SIMU/dt))

# Nombre de périodes où les points doivent être surveillés AU MOINS 1 fois
	global PERIOD = convert(Int,TEMPS_SIMU/DELTA_T)

# Nombre de pas de temps dans une période
	global NUM_PERIOD = convert(Int,round(DELTA_T/dt))
	
# Earth's rotation rate in rad/s
	global we= 7.2921e-5

# Nombre de points à observer sur la Terre (nous rappelons que la Terre est une sphère discrétisée)
	global RANGE_NUM_PIXEL = 1:1#1*fill(1, 10)#1:2
	
# Obtention ensemble de valeurs admis pour hauteur satellite en fonction du nombre de N revolutions du satellite
	N_test=1:24
	global altitudeSet = Float64[]
	global PERIOD_SAT = Float64[]
	for i in N_test
		altitudeSetVal = ( ( μ*((TEMPS_SIMU)/N_test[i])^2 ) / 4π^2 )^(1/3) - RAYON
		if altitudeSetVal >= 400000 && altitudeSetVal <= 1400000
			append!(altitudeSet, [altitudeSetVal])
			append!(PERIOD_SAT, [TEMPS_SIMU/N_test[i]])
		end
	end
	# println(altitudeSet)
	
	global PERIOD_SAT_MIN = PERIOD_SAT[end]

# Initialize variable in model
	global start_variable = 0.0
	global AltitudeVariable = altitudeSet[1]
	global PERIODVariable = PERIOD_SAT[1]

	global Modele = Val(:twobody)
	global color = ["brown","green","red","blue","violet","black"]
	global M_limit = 2*pi
	
	global time_zero_simulation = date_to_jd(2000,1,1,12,0,0)#1970 value: 2.451545e6

	# println(time_zero_simulation)
	# println(jd_to_date(time_zero_simulation))
	

