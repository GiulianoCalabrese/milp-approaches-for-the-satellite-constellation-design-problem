include("Parameters.jl")
include("Functions.jl")

using SatelliteToolbox, PyPlot

global TimeSeen_track = zeros(length(RANGE_NUM_PIXEL),length(collect(0:dt:TEMPS_SIMU)))
  
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
		global theta[j] = -alpha_lim + asin(((RAYON+altitudeSet[j])/RAYON)*sin(alpha_lim ))
	else	
		global theta[j] = (-alphaHalf + asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf)))
	end
end

println("ThetaMin :"*string(ThetaMin*(180/pi)))
println("alpha "*string(alphaHalf*(180/pi)))
println("ThetaSet :"*string(theta*(180/pi)))
println("periode : "*string(periodSat))

number_targets = RANGE_NUM_PIXEL[end]
global x_target,y_target,z_target = generateTargetsCartesians(number_targets)
ρ = Float64[]
lat = Float64[]
long = Float64[]

global legend_plot = []	

CreateSfera(RAYON)
CreateEquatorialPlane()
	
# Creazione targets
for NUM_PIXEL in RANGE_NUM_PIXEL
	local Results=[]
	ρ_temp,colat_temp,long_temp = CartesianToPolar(x_target[NUM_PIXEL],y_target[NUM_PIXEL],z_target[NUM_PIXEL])
	append!(ρ,ρ_temp)
	append!(lat,(pi/2-colat_temp))
	append!(long,long_temp)
	# PyPlot.title("Reference system : ECEF")
	append!(legend_plot, ["Target "*string(NUM_PIXEL)])
	PyPlot.scatter3D(x_target[NUM_PIXEL],y_target[NUM_PIXEL],z_target[NUM_PIXEL], color = color[2+NUM_PIXEL], label = legend_plot[NUM_PIXEL], marker = "o")#, markersize = 5)
	println("lat,long :"*string([lat[NUM_PIXEL],long[NUM_PIXEL]]))
end

# Discretization Keplerian parameters
inclinaisonSet = Discretization(numbOfDiscretize,1*pi)
noeudAscendantSet = Discretization(numbOfDiscretize,2*pi)
meanAnomalySet = copy(noeudAscendantSet)
		
index_semiAxis=1
k=6
l=1
s=1
index_semiAxis =   1
k =  1
l =  2
s =  1
semiAxis = RAYON+altitudeSet[index_semiAxis]

println(semiAxis)
println(meanAnomalySet[k])
println(noeudAscendantSet[l])
println(inclinaisonSet[s])

#creazione propagazione cartesiane dei satelliti nel tempo	
for i in 1:length(index_semiAxis)	

	local LatitudeSat,LongitudeSat = ProjSatPlotPOLARDominique(semiAxis, inclinaisonSet[s], noeudAscendantSet[l], meanAnomalySet[k])
	println("theta :"*string(theta[index_semiAxis]))
	println("a,i,nA,mA :"*string([semiAxis,inclinaisonSet[s],noeudAscendantSet[l],meanAnomalySet[k]]))

	local orbp = init_orbit_propagator(Modele, time_zero_simulation, semiAxis, 0.0, inclinaisonSet[s], noeudAscendantSet[l],0.0, meanAnomalySet[k])
	local r,v = propagate!(orbp, collect(0:dt:TEMPS_SIMU)) 

	local x_sat=fill(0.0, length(r)) 
	local y_sat=fill(0.0, length(r)) 
	local z_sat=fill(0.0, length(r))
	local x_polar = []
	local y_polar = []
	local z_polar = [] 
	
	# for indexTarget in RANGE_NUM_PIXEL
		# println("lat:"*string(lat[indexTarget]))
		# println("long:"*string(long[indexTarget]))
	# end
	# println("theta:"*string(theta[index_semiAxis]))
	for j in 1:(length(r)-1)
		# matrix rotation for ECI to ECEF
		local matrix_rotation = transpose(r_eci_to_ecef(TOD(), PEF(), (dt*(j-1))/86400))

		local ECItoECEF = matrix_rotation*(r[j][:])
		
		local x_sat[j]=ECItoECEF[1]
		local y_sat[j]=ECItoECEF[2]
		local z_sat[j]=ECItoECEF[3]
		
		x_temp,y_temp,z_temp = PolarToCartesian(semiAxis,pi/2-LatitudeSat[j],LongitudeSat[j])

		result_temp = [x_temp,y_temp,z_temp]
		append!(x_polar,[result_temp[1]])
		append!(y_polar,[result_temp[2]])
		append!(z_polar,[result_temp[3]])
		
		for indexTarget in RANGE_NUM_PIXEL
			# println("diff lat:"*string(abs(LatitudeSat[j]-lat[indexTarget])))
			# println("diff long:"*string(abs(LongitudeSat[j]-long[indexTarget])))
			if 	(((theta[index_semiAxis] >= abs(LatitudeSat[j]-lat[indexTarget])))
			&& 
				((theta[index_semiAxis] >= abs(LongitudeSat[j]-long[indexTarget])) || (theta[index_semiAxis] >= abs(LongitudeSat[j]-long[indexTarget]-2*pi))))			
				# if indexTarget == 1
				# println("target  "*string(indexTarget)*"   "*string(j*dt))	
				# end
				TimeSeen_track[indexTarget,j] = 1.0
			end
		end

	end

	# Position sat en ECEF
	#PyPlot.scatter3D(x_sat,y_sat,z_sat,color = "green")#, label = legend_plot[end-1])
	
	#PyPlot.scatter3D(x_polar,y_polar,z_polar,color = "magenta")#, label = legend_plot[end])

end	
 
# append!(legend_plot,["Optimal Satellite_ECEF"])
# append!(legend_plot,["Optimal Satellite_ECI"])
legend()
	
CountRevisitTime(TimeSeen_track)

 