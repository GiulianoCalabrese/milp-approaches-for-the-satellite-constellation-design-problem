import Pkg
using JuMP, AmplNLWriter, Ipopt
include("Parameters.jl")
include("Functions.jl")

theta_min = ((2*pi*dt)/(2*PERIOD_SAT_MIN))*1.2
Theta = fill(0.0,length(altitudeSet))
# demi-angle d'ouverture du satellite choisit arbitrairement en degrés
global alphaHalf = atan((sin(theta_min))/(((RAYON+altitudeSet[end])/RAYON)-cos(theta_min)))
for j in 1:length(altitudeSet)
	if ((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf) > 1
		global Theta[j] = pi/2 - alphaHalf
	else	
		global Theta[j] = (-alphaHalf + asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf)))
	end
end
println(Theta)

global name_file = "MINLP-routine_GEO_theta_"*string(round(theta_min*180/pi))*"_Tsimu"*string(Int(TEMPS_SIMU/3600))*"H_DeltaT"*string(Int(DELTA_T/3600))*"H_dt"*string(Int(dt/60))*"Mn.txt"
println(name_file)

global theta0 = calcul_angle_ECI_ECEF(time_zero_simulation)

open(name_file, "w") do Output
	for NUM_PIXEL in RANGE_NUM_PIXEL

		status = "nothing"
		println("nombre target"*string(NUM_PIXEL))

		# Obtenir coordonnées polaires des cartésiennes
		x,y,z = generateTargetsCartesians(NUM_PIXEL)

		# Position pixel statique, référentiel ECEF
		write(Output," coordonnées cibles : ")
		ρ = Float64[]
		lat = Float64[]
		long = Float64[]
		for i in 1:NUM_PIXEL
			ρ_temp,colat_temp,long_temp = CartesianToPolar(x[i],y[i],z[i])
			append!(ρ,ρ_temp)
			append!(lat,pi/2-colat_temp)
			append!(long,long_temp+pi)
		end
		write(Output," ρ : ",  string(ρ))
		write(Output," lat : ",  string(lat*(180/pi)))
		write(Output," long : ",  string(long*(180/pi)))
		write(Output,"\n"," ")

		# Initialisation du nombre de satellites qui peuvent observer la Terre
		global NUM_SATELLITE = 0
		global solvedTime = 0.0

		while (string(status) != "LOCALLY_SOLVED") & (string(status) != "OPTIMAL") & (string(status) != "TIME_LIMIT") & (string(status) != "OTHER_LIMIT") & (NUM_SATELLITE < NUM_PIXEL*PERIOD)

			global NUM_SATELLITE = NUM_SATELLITE+1
			global M_limit_1 = NUM_PIXEL * NUM_SATELLITE * NUM_TIME
			println("Nombre Satellites : ",NUM_SATELLITE)
			write(Output,"NUM_SATELLITE : ", string(NUM_SATELLITE))
			write(Output,"\n"," ")

			# model = Model(() -> AmplNLWriter.Optimizer("./mac-osx/couenne"))
			#model = Model(() -> AmplNLWriter.Optimizer("./mac-osx/bonmin"))
			#model = Model(Ipopt.Optimizer)
			# model = Model(() -> AmplNLWriter.Optimizer("./Windows/couenne/couenne"))
			model = Model(() -> AmplNLWriter.Optimizer("./Windows/bonmin/bonmin"))
			#set_optimizer_attribute(model, "max_cpu_time", 3600.0)
			#set_optimizer_attribute(model, "max_iter", 50000)
			#set_optimizer_attribute(model, "acceptable_dual_inf_tol", 10e-3)
			#set_optimizer_attribute(model, "print_level", 0)
			JuMP.register(model, :calcul_angle_ECI_ECEF, 1, calcul_angle_ECI_ECEF; autodiff = true)
			JuMP.register(model, :atan, 2, atan; autodiff = true)
			JuMP.register(model, :tan, 1, atan, autodiff=true)
			JuMP.register(model, :sin, 1, sin, autodiff=true)
			JuMP.register(model, :asin, 1, asin, autodiff=true)
			JuMP.register(model, :cos, 1, cos, autodiff=true)
			JuMP.register(model, :acos, 1, acos, autodiff=true)
			JuMP.register(model, :%, 2, %; autodiff = true)

			# Altitude de l'orbite [Km]
			@variable(model, altitude[1:NUM_SATELLITE])
			@variable(model, activation_altitude[1:NUM_SATELLITE,1:length(altitudeSet)], Bin)

			# inclinaison de l'orbite
			@variable(model, 0 <= sin_inclinaison[1:NUM_SATELLITE] <= 1)
			@variable(model, -1 <= cos_inclinaison[1:NUM_SATELLITE] <= 1)

			# Noeud ascendant
			@variable(model, -1 <= sin_noeudAscendant[1:NUM_SATELLITE] <= 1)
			@variable(model, -1 <= cos_noeudAscendant[1:NUM_SATELLITE] <= 1)

			# Anomalie moyenne du plan orbit
			@variable(model, -1 <= sin_meanAnomaly[1:NUM_SATELLITE] <= 1)
			@variable(model, -1 <= cos_meanAnomaly[1:NUM_SATELLITE] <= 1)
			
			#Fonction objectif : Minimize the number of active satellite
			@objective(model, Min, 0)

			# Identité trigonométrique
			@constraint(model, [i=1:NUM_SATELLITE], ((sin_inclinaison[i]^2) + (cos_inclinaison[i]^2)) == 1)
			@constraint(model, [i=1:NUM_SATELLITE], ((sin_noeudAscendant[i]^2) + (cos_noeudAscendant[i]^2)) == 1)
			@constraint(model, [i=1:NUM_SATELLITE], ((sin_meanAnomaly[i]^2) + (cos_meanAnomaly[i]^2)) == 1)

			# Compute the geocentric distance.
			@constraint(model, [i=1:NUM_SATELLITE], sum(activation_altitude[i,j] for j = 1:length(altitudeSet))==1)
			@constraint(model, [i=1:NUM_SATELLITE], altitude[i] == RAYON + sum(altitudeSet[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))

			# Période de révolution de l'orbite [Km]
			@variable(model, periodSat[1:NUM_SATELLITE], start = PERIODVariable)
			@constraint(model, [i=1:NUM_SATELLITE], periodSat[i] == sum(PERIOD_SAT[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))

			# Thetamax permet de déterminer la surface couverte par le satellite à partir du demi angle d'ouverture des capteurs, de l'altitude et de l'excentricité
			@variable(model, Thetamax[1:NUM_SATELLITE])
			@constraint(model, [i=1:NUM_SATELLITE], Thetamax[i] == sum(Theta[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))
			
			# Variables de linéarisation (MINLP)
			@variable(model, Xi[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL], Bin)
			
			@NLexpression(model, LatitudeSat[i=1:NUM_SATELLITE, p=1:NUM_TIME],asin(((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*
			cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i]) + ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))*((p-1)*dt))*
			sqrt(altitude[i]/μ))))

			@variable(model, 0 <= LongitudeSat[1:NUM_SATELLITE,1:NUM_TIME] <= 2*pi)				
			@NLconstraint(model, [i=1:NUM_SATELLITE,p=1:NUM_TIME],LongitudeSat[i,p] == (-(theta0 + (we*((p-1)*dt))) + 
			atan((((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i]) + (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i])
			+ ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i]) + (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(sqrt(μ/(altitude[i]^3))*((p-1)*dt))
			*sqrt(altitude[i]/μ)),((((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i]) - (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
			cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i])+ ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
			cos_meanAnomaly[i]))*sin((sqrt(μ/(altitude[i]^3))*((p-1)*dt)))*sqrt(altitude[i]/μ))))))	
					
			@NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], LatitudeSat[i,p] - lat[j] <= Thetamax[i]+(1-Xi[i,p,j])*M_limit)
			@NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], - LatitudeSat[i,p] + lat[j] <= Thetamax[i]+(1-Xi[i,p,j])*M_limit)
			@NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], LongitudeSat[i,p] - long[j] <= Thetamax[i]+(1-Xi[i,p,j])*M_limit)
			@NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], - LongitudeSat[i,p] + long[j] <= Thetamax[i]+(1-Xi[i,p,j])*M_limit)

			# Linearisation du modèle (approche MINLP)
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
				write(Output,"\n  sin_inclinaison = ",string(value.(model[:sin_inclinaison])))
				write(Output,"\n  cos_inclinaison = ",string(value.(model[:cos_inclinaison])))
				write(Output,"\n  sin_noeudAscendant = ",string(value.(model[:sin_noeudAscendant])))
				write(Output,"\n  cos_noeudAscendant = ",string(value.(model[:cos_noeudAscendant])))
				write(Output,"\n  sin_meanAnomaly = ",string(value.(model[:sin_meanAnomaly])))
				write(Output,"\n  cos_meanAnomaly = ",string(value.(model[:cos_meanAnomaly])))
				write(Output,"\n  theta_max_i = ",string((value.(model[:Thetamax])*180)/pi))
				write(Output,"\n"," ")
				for i=1:NUM_SATELLITE
					println(value.(model[:Thetamax])[i][1]*(180/pi))
					for j=1:NUM_PIXEL
						for p=1:NUM_TIME
							if(((abs(value.(model[:LatitudeSat])[i,p][1]-lat[j])) <= (value.(model[:Thetamax])[i][1])) && ((abs(value.(model[:LongitudeSat])[i,p][1]-long[j])) <= (value.(model[:Thetamax])[i][1])))
								# println(value.(model[:LatitudeSat])[i,p])							
								# println(((value.(model[:Thetamax])[i][1])>= abs(value.(model[:LatitudeSat])[i,p]-lat[j])))
								println("observing time (Theta): ",p*dt)
								#println("target  "*string(j)*"   "*string((value.(model[:LatitudeSat])[i,p][1])))
							end
							if(value.(model[:Xi])[i,p,j][1] > 0)
								println("observing time (Xi): ",p*dt)
							end
						end
					end
				end
			end
			write(Output,"\n"," ")
		end
	end
end
