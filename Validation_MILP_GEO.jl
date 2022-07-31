include("Parameters.jl")
include("Functions.jl")

using SatelliteToolbox, PyPlot

lat = [14,-37, -57]*(pi/180)
long = [206, 18, 56]*(pi/180)

global TimeSeen_track = zeros(length(lat),length(collect(0:dt:TEMPS_SIMU)))  
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

global legend_plot = []	

CreateSfera(RAYON)
CreateEquatorialPlane()

for i in 1: length(lat)
	x_temp,y_temp,z_temp = PolarToCartesian(RAYON,pi/2-lat[i],long[i])
	PyPlot.scatter3D(x_temp,y_temp,z_temp, label = "Target "*string(i))
end

# Discretization Keplerian parameters
inclinaisonSet = Discretization(numbOfDiscretize,1*pi)
noeudAscendantSet = Discretization(numbOfDiscretize,2*pi)
meanAnomalySet = copy(noeudAscendantSet)

index_semiAxis = [1]
k =  [6]
l =  [23]
s =  [13]
  
#creazione propagazione cartesiane dei satelliti nel tempo	
for i in 1:length(index_semiAxis)	
	local semiAxis = RAYON+altitudeSet[index_semiAxis[i]]
	local LatitudeSat,LongitudeSat = ProjSatPlotPOLARDominique(semiAxis, inclinaisonSet[s[i]], noeudAscendantSet[l[i]], meanAnomalySet[k[i]])
	println("theta :"*string(theta[i]))
	println("a,i,nA,mA :"*string([semiAxis,inclinaisonSet[s[i]],noeudAscendantSet[l[i]],meanAnomalySet[k[i]]]))

	local orbp = init_orbit_propagator(Modele, time_zero_simulation, semiAxis, 0.0, inclinaisonSet[s[i]], noeudAscendantSet[l[i]],0.0, meanAnomalySet[k[i]])
	local r,v = propagate!(orbp, collect(0:dt:TEMPS_SIMU)) 

	local x_sat=fill(0.0, length(r)) 
	local y_sat=fill(0.0, length(r)) 
	local z_sat=fill(0.0, length(r))
	local x_polar = []
	local y_polar = []
	local z_polar = [] 
	
	for j in 1:(length(r))
		# matrix rotation for ECI to ECEF
		local matrix_rotation = r_eci_to_ecef(TOD(), PEF(), time_zero_simulation+(dt*(j-1))/86400)
		local ECItoECEF = matrix_rotation*(r[j][:])
		
		local x_sat[j]=ECItoECEF[1]
		local y_sat[j]=ECItoECEF[2]
		local z_sat[j]=ECItoECEF[3]
		
		x_temp,y_temp,z_temp = PolarToCartesian(semiAxis,pi/2-LatitudeSat[j],LongitudeSat[j])
		append!(x_polar,[x_temp])
		append!(y_polar,[y_temp])
		append!(z_polar,[z_temp])
		
		# println("x :"*string([x_polar[j],x_sat[j]]))
		# println("y :"*string([y_polar[j],y_sat[j]]))
		# println("z :"*string([z_polar[j],z_sat[j]]))
		
		for indexTarget in 1:length(lat)
			# println("diff lat:"*string(abs(LatitudeSat[j]-lat[indexTarget])))
			# println("diff long:"*string(abs(LongitudeSat[j]-long[indexTarget])))
			if 	((theta[i] >= abs(LatitudeSat[j]-lat[indexTarget])) && (theta[i] >= abs(LongitudeSat[j]-long[indexTarget])))			
				TimeSeen_track[indexTarget,j] = 1.0
			end
		end

	end

	# Position sat en ECEF
	PyPlot.plot(x_sat,y_sat,z_sat,color = "magenta")#, label = legend_plot[end-1])
	# PyPlot.scatter3D(x_polar,y_polar,z_polar,color = "magenta")#, label = legend_plot[end])

end	
 
# append!(legend_plot,["Optimal Satellite_ECEF"])
# append!(legend_plot,["Optimal Satellite_ECI"])
legend()	
CountRevisitTime(TimeSeen_track)