using JuMP, CPLEX, Random
include("Parameters.jl")
include("Functions.jl")
# Random.seed!(1)

global theta_min = ((2*pi*dt)/(2*PERIOD_SAT_MIN))*1.2
global numbOfDiscretize = Int(ceil(pi/theta_min))
global alphaHalf = atan((sin(theta_min))/(((RAYON+altitudeSet[end])/RAYON)-cos(theta_min)))

Theta = fill(0.0,length(altitudeSet))
for j in 1:length(altitudeSet)
	if ((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf) > 1
		alpha_lim = asin(RAYON / (RAYON + altitudeSet[j]))
		global Theta[j] = (-alpha_lim + asin(round(((RAYON+altitudeSet[j])/RAYON)*sin(alpha_lim),digits=6)))
	else	
		global Theta[j] = (-alphaHalf + asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf)))
	end
end

# File name of output file	
global name_file = "MILP-Routine_GEO_MULTITARGET_"*string(numbOfDiscretize)*"_theta_"*string(round(theta_min*180/pi))*"_Tsimu"*string(Int(TEMPS_SIMU/3600))*"H_DeltaT"*string(Int(DELTA_T/3600))*"H_dt"*string(Int(dt/60))*"Mn.txt" 

# Les satellites couvrent une zone de la Terre qui corresponds à un carré, RAYON_SAT c'est la longueur du DEMI-coté
println([numbOfDiscretize, PERIOD_SAT_MIN])
println(Theta*(180/pi))
println(alphaHalf*(180/pi))
println(name_file)

global lat_temp = [80, 61, 33, 22, 16, 8, 1, -13, -15, -24, -31, -38, -53]*(pi/180)
global long_temp = [358, 116, 139, 343, 340, 55, 220, 24, 100, 9, 152, 144, 46]*(pi/180)
global ρ_temp = [6.3781363e6, 6.3781363e6, 6.3781363e6, 6.3781363e6, 6.3781363e6, 6.3781363e6, 6.3781363e6, 6.3781363e6, 6.3781363e6, 6.3781363e6]
global instance = 0

open(name_file, "w") do Output
	for NUM_PIXEL in RANGE_NUM_PIXEL
		println("Nombre cibles : ",NUM_PIXEL)
	
		# Obtenir coordonnées polaires des cartésiennes		
		# x,y,z = generateTargetsCartesians(NUM_PIXEL)
		
		# Position pixel statique, référentiel ECEF
		write(Output," coordonnées cibles : ")
		ρ = Float64[]
		lat = Float64[]
		long = Float64[]
		# for i in 1:NUM_PIXEL
			# ρ_temp,colat_temp,long_temp = CartesianToPolar(x[i],y[i],z[i])
			# append!(ρ,ρ_temp)
			# append!(lat,pi/2-colat_temp)
			# # longitude doit être comprise entre 0:2pi sinon incohérence avec formules Dominique
			# append!(long,long_temp+pi)
		# end

		for i in 1:NUM_PIXEL
			append!(ρ,ρ_temp[i])
			append!(lat,lat_temp[instance+i])
			append!(long,long_temp[instance+i])
		end	
		
		global instance += 1
		
		println([lat*(180/pi),long*(180/pi)])
		write(Output," ρ : ",  string(ρ))
		write(Output," lat : ",  string(lat*(180/pi)))
		write(Output," long : ",  string(long*(180/pi)))
		write(Output,"\n"," ")

		# Initialisation du nombre de satellites qui peuvent observer la Terre
		global NUM_SATELLITE = 0
		global solvedTime = 0.0
		
		# Discretization Keplerian parameters
		inclinaisonSet = Discretization(numbOfDiscretize,1*pi)
		noeudAscendantSet = Discretization(numbOfDiscretize,2*pi)
		meanAnomalySet = copy(noeudAscendantSet)
		
		n_inclinaison = numbOfDiscretize
		n_noeudAscendant = numbOfDiscretize
		n_meanAnomaly = numbOfDiscretize
		
		@time begin			
			CoverageSatLat = fill(0.0,NUM_PIXEL,NUM_TIME,numbOfDiscretize,numbOfDiscretize,numbOfDiscretize,length(altitudeSet))
			CoverageSatLong = fill(0.0,NUM_PIXEL,NUM_TIME,numbOfDiscretize,numbOfDiscretize,numbOfDiscretize,length(altitudeSet))
			
			for j in 1:length(altitudeSet)
				for k in 1:numbOfDiscretize		
					for l in 1:numbOfDiscretize	
						for s in 1:numbOfDiscretize
							global lat_Sat,long_Sat = ProjSatPlotPOLARDominique(RAYON+altitudeSet[j], inclinaisonSet[k],
							noeudAscendantSet[l], meanAnomalySet[s])
							
							for p in 1:NUM_TIME
								for t in 1:NUM_PIXEL
									CoverageSatLat[t,p,s,l,k,j] = abs(lat[t] - lat_Sat[p])
									CoverageSatLong[t,p,s,l,k,j] = abs(long[t] - long_Sat[p])									
								end
							end
						end
					end	
				end
			end
		end		
		
		status = "nothing"
		
		while (string(status) != "LOCALLY_SOLVED") & (string(status) != "OPTIMAL") & (string(status) != "TIME_LIMIT") & (string(status) != "MEMORY_LIMIT") & (NUM_SATELLITE < NUM_PIXEL*PERIOD)
			println(status)
			global NUM_SATELLITE = NUM_SATELLITE+1 
			println("Nombre Satellites : ",NUM_SATELLITE)
			write(Output,"NUM_SATELLITE : ", string(NUM_SATELLITE))
			write(Output,"\n"," ")
			model = Model(CPLEX.Optimizer)
			# set_optimizer_attribute(model, "CPX_PARAM_TILIM", 3600)

			# Altitude de l'orbite [Km]
			@variable(model, altitude[1:NUM_SATELLITE], start = AltitudeVariable)
			@variable(model, activation_altitude[1:NUM_SATELLITE,1:length(altitudeSet)], Bin, start = start_variable)
			
			# Inclinaison de l'orbite
			@variable(model, inclinaison[1:NUM_SATELLITE], start = 0)
			@variable(model, activation_inclinaison[1:NUM_SATELLITE,1:n_inclinaison], Bin, start = start_variable)	
			
			# Noeud ascendant
			@variable(model, noeudAscendant[1:NUM_SATELLITE], start = pi)
			@variable(model, activation_noeudAscendant[1:NUM_SATELLITE,1:n_noeudAscendant], Bin, start = start_variable)

			# Anomalie moyenne du plan orbit
			@variable(model, meanAnomaly[1:NUM_SATELLITE], start = 0)
			@variable(model, activation_meanAnomaly[1:NUM_SATELLITE,1:n_meanAnomaly], Bin, start = start_variable)			 
			@variable(model, η==0)
		
			# Variables de linéarisation 
			@variable(model, Xi[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL], Bin, start = start_variable)
		
			# Fonction objectif : Minimize the number of active satellite (choice between the next 2:)
			# minimiser sur nombre de sat fixés
			@objective(model, Min, η) 

			# Compute the geocentric distance.
			@constraint(model, [i=1:NUM_SATELLITE], sum(activation_altitude[i,j] for j = 1:length(altitudeSet))==1)
			@constraint(model, [i=1:NUM_SATELLITE], altitude[i] == RAYON + sum(altitudeSet[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))
			# Verification sinus/ cosinus relationship for angle inclinaison
			@constraint(model, [i=1:NUM_SATELLITE], sum(activation_inclinaison[i,j] for j = 1:n_inclinaison)==1)
			@constraint(model, [i=1:NUM_SATELLITE], inclinaison[i] == sum(inclinaisonSet[j]*activation_inclinaison[i,j] for j=1:n_inclinaison))
			# Verification sinus/ cosinus relationship for angle RAAN
			@constraint(model, [i=1:NUM_SATELLITE], sum(activation_noeudAscendant[i,j] for j = 1:n_noeudAscendant)==1)
			@constraint(model, [i=1:NUM_SATELLITE], noeudAscendant[i] == sum(noeudAscendantSet[j]*activation_noeudAscendant[i,j] for j=1:n_noeudAscendant))
			# Verification sinus/ cosinus relationship for angle Mean Anomaly
			@constraint(model, [i=1:NUM_SATELLITE], sum(activation_meanAnomaly[i,j] for j = 1:n_meanAnomaly)==1)
			@constraint(model, [i=1:NUM_SATELLITE], meanAnomaly[i] == sum(meanAnomalySet[j]*activation_meanAnomaly[i,j] for j=1:n_meanAnomaly))

			# Pour la linéarisation des variables en utilisant des variables binaires
			@variable(model, act4[i=1:NUM_SATELLITE, j=1:length(altitudeSet), k=1:n_meanAnomaly, l=1:n_noeudAscendant, s=1:n_inclinaison], Bin, start = start_variable)
			@constraint(model, [i=1:NUM_SATELLITE, j=1:length(altitudeSet), k=1:n_meanAnomaly, l=1:n_noeudAscendant, s=1:n_inclinaison], act4[i,j,k,l,s] <= activation_altitude[i,j])
			@constraint(model, [i=1:NUM_SATELLITE, j=1:length(altitudeSet), k=1:n_meanAnomaly, l=1:n_noeudAscendant, s=1:n_inclinaison], act4[i,j,k,l,s] <= activation_meanAnomaly[i,k])
			@constraint(model, [i=1:NUM_SATELLITE, j=1:length(altitudeSet), k=1:n_meanAnomaly, l=1:n_noeudAscendant, s=1:n_inclinaison], act4[i,j,k,l,s] <= activation_inclinaison[i,s])
			@constraint(model, [i=1:NUM_SATELLITE, j=1:length(altitudeSet), k=1:n_meanAnomaly, l=1:n_noeudAscendant, s=1:n_inclinaison], act4[i,j,k,l,s] <= activation_noeudAscendant[i,l])
			@constraint(model, [i=1:NUM_SATELLITE, j=1:length(altitudeSet), k=1:n_meanAnomaly, l=1:n_noeudAscendant, s=1:n_inclinaison], act4[i,j,k,l,s] >= 
			activation_altitude[i,j] + activation_meanAnomaly[i,k] + activation_inclinaison[i,s] + activation_noeudAscendant[i,l] - 3)

			# Thetamax permet de déterminer la surface couverte par le satellite à partir du demi angle d'ouverture des capteurs, de l'altitude et de l'excentricité
			@variable(model, Thetamax[1:NUM_SATELLITE])
			@constraint(model, [i=1:NUM_SATELLITE], Thetamax[i] == sum(Theta[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))
			# @constraint(model, [i=1:NUM_SATELLITE], Thetamax[i] == (-alphaHalf + sum(activation_altitude[i,j]*asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf)) for j=1:length(altitudeSet))))
					
			@variable(model, ThetaTargetLat[I=1:NUM_SATELLITE,j=1:NUM_PIXEL, p=1:NUM_TIME])
			@constraint(model, [i=1:NUM_SATELLITE,j=1:NUM_PIXEL, p=1:NUM_TIME], ThetaTargetLat[i,j,p] == 
			sum(act4[i,a,k,l,s]*CoverageSatLat[j,p,k,l,s,a]
			for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))) 
			
			@variable(model, ThetaTargetLong[i=1:NUM_SATELLITE,j=1:NUM_PIXEL, p=1:NUM_TIME])
			@constraint(model, [i=1:NUM_SATELLITE,j=1:NUM_PIXEL, p=1:NUM_TIME], ThetaTargetLong[i,j,p] == 
			sum(act4[i,a,k,l,s]*CoverageSatLong[j,p,k,l,s,a]
			for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))) 

			@constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi[i,p,j] => {-ThetaTargetLat[i,j,p] <= Thetamax[i]})
			@constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi[i,p,j] => {ThetaTargetLat[i,j,p] <= Thetamax[i]})
			@constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi[i,p,j] => {-ThetaTargetLong[i,j,p] <= Thetamax[i]})
			@constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi[i,p,j] => {ThetaTargetLong[i,j,p] <= Thetamax[i]})
			
			# Linearisation du modèle 
			@constraint(model, [j=1:NUM_PIXEL,k=1:PERIOD], sum(Xi[i,p,j] for p in (((k-1)*NUM_PERIOD)+1):(NUM_PERIOD*k) for i in 1:NUM_SATELLITE) >= 1)	 
						
			optimize!(model)
			
			status = termination_status(model)
			println("status : ", status)
			write(Output," status : ", string(status)) 
			write(Output,";  time CPU :", string(solve_time(model))) 
			write(Output,"\n  dt : ",string(dt))
			write(Output,";  alpha : ",string(alphaHalf*180/pi))
			write(Output,";  theta_min : ",string(theta_min*180/pi))
			if string(status) == "LOCALLY_SOLVED" || string(status) == "OPTIMAL"
				write(Output,"\n  demi-axe = ", string((value.(model[:altitude]))))
				write(Output,"\n  meanAnomaly = ",string(value.(model[:meanAnomaly])))
				write(Output,"\n  noeudAscendant = ",string(value.(model[:noeudAscendant])))
				write(Output,"\n  inclinaison = ",string(value.(model[:inclinaison])))
				write(Output,"\n  theta_max_i = ",string((value.(model[:Thetamax])*180)/pi))
				write(Output,"\n  index_semiAxis =  ",string(findall(a->a==1,value.(model[:activation_altitude]))))
				write(Output,"\n  k =  ",string(findall(a->a==1,value.(model[:activation_meanAnomaly])[:,:])))
				write(Output,"\n  l =  ",string(findall(a->a==1,value.(model[:activation_noeudAscendant])[:,:])))
				write(Output,"\n  s =  ",string(findall(a->a==1,value.(model[:activation_inclinaison])[:,:])))
				write(Output,"\n"," ")
				for i=1:NUM_SATELLITE
					for j=1:NUM_PIXEL
						for p=1:NUM_TIME
							if (((((value.(model[:Thetamax])[i][1]) >= (value.(model[:ThetaTargetLat])[i,j,p][1]))))
							&& 
							   ((value.(model[:Thetamax])[i][1]) >= (value.(model[:ThetaTargetLong])[i,j,p][1])) )			
								println("observing time (Theta): ",p*dt)
							end
							if(value.(model[:Xi])[i,p,j][1] > 0)
								println("observing time (Xi): ",p*dt)
							end
						end
					end
				end
			end
			println(status)
			write(Output,"\n"," ")
		end
	end			
end