include("Parameters.jl")
include("Functions.jl")

using SatelliteToolbox, PyPlot
  
periodSat = 2*pi*sqrt(((RAYON+altitudeSet[end])^3)/μ)
global ThetaMin = ((2*pi*dt)/(2*PERIOD_SAT_MIN))*1.2
global numbOfDiscretize = Int(ceil(pi/ThetaMin))
alphaHalf = atan((sin(ThetaMin))/(((RAYON+altitudeSet[end])/RAYON)-cos(ThetaMin)))
theta = fill(0.0,length(altitudeSet))
for j in 1:length(altitudeSet)
	if ((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf) > 1
		# recalculer alpha limite
		alpha_lim = asin(RAYON / (RAYON + altitudeSet[j]))
		#calculer theta avec alpha limite
		global theta[j] = -alpha_lim + asin(1.0)#((RAYON+altitudeSet[j])/RAYON)*sin(alpha_lim ))
	else	
		global theta[j] = (-alphaHalf + asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf)))
	end
end

println("ThetaMin :"*string(ThetaMin*(180/pi)))
println("alpha "*string(alphaHalf*(180/pi)))
println("ThetaSet :"*string(theta*(180/pi)))
println("periode : "*string(periodSat))

lat = [-55.156239210990506]*(pi/180)
long = [116.51310287186844]*(pi/180)
RANGE_NUM_PIXEL = 1:length(lat)
global TimeSeen_track = zeros(length(RANGE_NUM_PIXEL),length(collect(0:dt:TEMPS_SIMU)))
CreateSfera(RAYON)
CreateEquatorialPlane()	
global legend_plot = []	
	
inclinaison = []
noeudAscendant = []
meanAnomaly = []

altitude = [7.640227192785597e6]

sin_inclinaison = [0.7108923944269709]
cos_inclinaison = [0.7033007916573732]
sin_noeudAscendant = [0.276814503693488]
cos_noeudAscendant = [-0.9609233739195484]
sin_meanAnomaly = [-0.9426898353797738]
cos_meanAnomaly = [-0.3336703077465163]

# sin_inclinaison = [sin(1.33)]
# cos_inclinaison = [cos(1.33)]
# sin_noeudAscendant = [sin(0.33)]
# cos_noeudAscendant = [cos(0.33)]
# sin_meanAnomaly = [sin(0.0)]
# cos_meanAnomaly = [cos(0.0)]

println(sin_inclinaison[1]^2 + cos_inclinaison[1]^2)

for i in 1:length(altitude)
	temp = asin(sin_inclinaison[i])
	append!(inclinaison, temp)
	temp = asin(sin_noeudAscendant[i])
	append!(noeudAscendant, temp)
	temp = asin(sin_meanAnomaly[i])
	append!(meanAnomaly, temp)
end
  
#creazione propagazione cartesiane dei satelliti nel tempo	
for i in 1:length(altitude)		

	local orbp = init_orbit_propagator(Modele, time_zero_simulation, altitude[i], 0.0, inclinaison[i], noeudAscendant[i],
				 0.0, meanAnomaly[i])
	local r,v = propagate!(orbp, collect(0:dt:TEMPS_SIMU)) 

	local x_sat=fill(0.0, length(r)) 
	local y_sat=fill(0.0, length(r)) 
	local z_sat=fill(0.0, length(r))
	local x_polar = []
	local y_polar = []
	local z_polar = []
	local LatitudeSat=fill(0.0, NUM_TIME+1) 
	local LongitudeSat=fill(0.0, NUM_TIME+1) 

	for j in 1:(length(r))

		# matrix rotation for ECI to ECEF
		local matrix_rotation = r_eci_to_ecef(TOD(), PEF(), time_zero_simulation+(dt*(j-1))/86400)
		local ECItoECEF = matrix_rotation*(r[j][:])
		
		local x_sat[j]=ECItoECEF[1]
		local y_sat[j]=ECItoECEF[2]
		local z_sat[j]=ECItoECEF[3]
		
		LatitudeSat[j] = asin(round(((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*cos(sqrt(μ/(altitude[i]^3))*((j-1)*dt))/altitude[i]) 
		+ ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))*((j-1)*dt))*sqrt(altitude[i]/μ)), digits=8))
		
		LongitudeSat[j] = ((-(calcul_angle_ECI_ECEF(time_zero_simulation) + (we*((j-1)*dt))) +
		atan((((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i]) + (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(sqrt(μ/(altitude[i]^3))*((j-1)*dt))/altitude[i])
		+ ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i]) + (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(sqrt(μ/(altitude[i]^3))*((j-1)*dt))
		*sqrt(altitude[i]/μ)),((((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i]) - (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
		cos(sqrt(μ/(altitude[i]^3))*((j-1)*dt))/altitude[i])+ ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
		cos_meanAnomaly[i]))*sin((sqrt(μ/(altitude[i]^3))*((j-1)*dt)))*sqrt(altitude[i]/μ))))) %(2*pi))# + (2*pi))%(2*pi))
		
		if LongitudeSat[j] <= 0 
			LongitudeSat[j] += 2*pi
		end
		
		x_temp,y_temp,z_temp = PolarToCartesian(altitude[i],pi/2-LatitudeSat[j],LongitudeSat[j])
		append!(x_polar,[x_temp])
		append!(y_polar,[y_temp])
		append!(z_polar,[z_temp])
		
		for indexTarget in RANGE_NUM_PIXEL
			# println(LongitudeSat[j])
			# println(long[indexTarget])
			if 	(theta[i] >= abs(LatitudeSat[j]-lat[indexTarget])) && (theta[i] >= abs(LongitudeSat[j]-long[indexTarget]))			
				# println("target  "*string(indexTarget)*"   "*string(j*dt))
				TimeSeen_track[indexTarget,j] = 1.0
			end
		end

	end

	# Position sat en ECEF
	# PyPlot.scatter3D(x_sat,y_sat,z_sat,color = "green")#, label = legend_plot[end-1])
	PyPlot.scatter3D(x_polar,y_polar,z_polar,color = "magenta")#, label = legend_plot[end])

end	
 
# append!(legend_plot,["Optimal Satellite_ECEF"])
# append!(legend_plot,["Optimal Satellite_ECI"])
legend()
	
CountRevisitTime(TimeSeen_track)

 