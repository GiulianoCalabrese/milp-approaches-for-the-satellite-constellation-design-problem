import Pkg
using JuMP, AmplNLWriter, Ipopt, CPLEX, Juniper, Alpine, Pavito, SCIP, BARON, NEOSServer
using Distributions, Random

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

global name_file = "MINLP_GEO_theta_"*string(round(theta_min*180/pi))*"_Tsimu"*string(Int(TEMPS_SIMU/3600))*"H_DeltaT"*string(Int(DELTA_T/3600))*"H_dt"*string(Int(dt/60))*"Mn.txt" 
println(name_file)

global angle0 = calcul_angle_ECI_ECEF(time_zero_simulation)

open(name_file, "w") do Output
	for NUM_PIXEL in RANGE_NUM_PIXEL

		#Initialisation du nombre de satellites qui peuvent observer la Terre
		global NUM_SATELLITE = NUM_PIXEL*PERIOD
		println("nombre target"*string(NUM_PIXEL))
		global M_limit_1 = NUM_PIXEL * NUM_SATELLITE * NUM_TIME

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

		# MIP optimizer
			cplex = optimizer_with_attributes(CPLEX.Optimizer) 
			# NLP optimizer
			ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "max_iter"=>1000000, "max_cpu_time"=>3600.0, "tol"=>0.9)
			# MINLP optimizer
			juniper = optimizer_with_attributes(Juniper.Optimizer, 
													 "nl_solver" => ipopt) 	
			
			# Global optimizer
			alpine = optimizer_with_attributes(Alpine.Optimizer, 
													 "nlp_solver" => ipopt,
													 "mip_solver" => cplex,
													 "minlp_solver" => juniper)
			pavito = optimizer_with_attributes(Pavito.Optimizer, 
													 "cont_solver" => ipopt,
													 "mip_solver" => cplex)
			# model = Model(alpine)
			# model = Model(juniper)
			# model = Model(pavito)
			# model = Model(SCIP.Optimizer) 
			# model = Model(BARON.Optimizer) 

			# model = Model() do 
				# NEOSServer.Optimizer(email="luca.mencarelli.university@gmail.com", solver="OCTERACT")
			# end
			# model = Model(() -> AmplNLWriter.Optimizer("D:/jufloquet/Desktop/RAPOSa_3.0.0_Windows"))
			# model = Model(() -> AmplNLWriter.Optimizer("./Windows/couenne/couenne"))
			model = Model(() -> AmplNLWriter.Optimizer("./Windows/bonmin/bonmin"))
			# model = Model(() -> AmplNLWriter.Optimizer("D:/jufloquet/Octeract/bin/octeract-engine"))
			
			#set_optimizer_attribute(model, "max_cpu_time", 3600.0)
			#set_optimizer_attribute(model, "max_iter", 50000)
			#set_optimizer_attribute(model, "acceptable_dual_inf_tol", 10e-3)
			#set_optimizer_attribute(model, "print_level", 1)

			# Altitude de l'orbite [Km]
			@variable(model, altitude[1:NUM_SATELLITE] >= RAYON, start=RAYON)
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

			# Thetamax permet de déterminer la surface couverte par le satellite à partir du demi angle d'ouverture des capteurs, de l'altitude et de l'excentricité
			@variable(model, Thetamax[1:NUM_SATELLITE])
			@constraint(model, [i=1:NUM_SATELLITE], Thetamax[i] == sum(Theta[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))
			
			# Variables de linéarisation (MINLP)
			@variable(model, Xi[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL], Bin)

			@variable(model, -pi/2 <= LatitudeSat[1:NUM_SATELLITE, 1:NUM_TIME] <= pi/2)
			@NLconstraint(model, [i=1:NUM_SATELLITE, p=1:NUM_TIME], sin(LatitudeSat[i,p]) == (((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*
			cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i]) + ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))
			*((p-1)*dt))*sqrt(altitude[i]/μ))))
			@NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], abs(- LatitudeSat[i,p] + lat[j]) <= Thetamax[i] + (1-Xi[i,p,j])*M_limit)
			
			@variable(model, LongitudeSat[1:NUM_SATELLITE, 1:NUM_TIME])
			@NLconstraint(model, [i=1:NUM_SATELLITE, p=1:NUM_TIME], (((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
			- (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i]) +
			((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
			cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ)))) * tan(LongitudeSat[i,p]) == ((((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])+
			(cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])+
			((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])	+ (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
			cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))))
			@variable(model, k[i=1:NUM_SATELLITE,p=1:NUM_TIME], Int)
			@variable(model, 0 <= y[i=1:NUM_SATELLITE,p=1:NUM_TIME] <= 2*pi-0.0001)
			@variable(model, k1[i=1:NUM_SATELLITE,p=1:NUM_TIME], Int)
			@variable(model, 0 <= y1[i=1:NUM_SATELLITE,p=1:NUM_TIME] <= 2*pi-0.0001)

			#https://math.stackexchange.com/questions/2513559/how-to-model-a-modulo-with-linear-constraints
			@constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], y1[i,p] - long[j] <= Thetamax[i] + (1-Xi[i,p,j])*M_limit)
			@constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], - y1[i,p] + long[j] <= Thetamax[i] + (1-Xi[i,p,j])*M_limit)
			@constraint(model, [p=1:NUM_TIME,i=1:NUM_SATELLITE], -(angle0 + (we*((p-1)*dt))) + LongitudeSat[i,p] == k[i,p] + (2*pi) * y[i,p])
			@constraint(model, [p=1:NUM_TIME,i=1:NUM_SATELLITE], y[i,p] + (2*pi) == k1[i,p] + (2*pi) * y1[i,p])

			# Linearisation du modèle (approche MINLP)
			@constraint(model, [j=1:NUM_PIXEL,k=1:PERIOD], sum(Xi[i,p,j] for p in (((k-1)*NUM_PERIOD)+1):(NUM_PERIOD*k) for i in 1:NUM_SATELLITE) >= 1)

			for index=1:20
				for i=1:NUM_SATELLITE
					start_variable_inc = rand(Uniform(0,1))
					set_start_value(sin_inclinaison[i], start_variable_inc)
					set_start_value(cos_inclinaison[i], sqrt(1-start_variable_inc^2))
					start_variable_noeud = rand(Uniform(-1,1))
					set_start_value(sin_noeudAscendant[i], start_variable_noeud)
					set_start_value(cos_noeudAscendant[i], sqrt(1-start_variable_noeud^2))
					start_variable_mean = rand(Uniform(-1,1))
					set_start_value(sin_meanAnomaly[i], start_variable_mean)
					set_start_value(cos_meanAnomaly[i], sqrt(1-start_variable_mean^2))
					set_start_value(altitude[i], RAYON+altitudeSet[1])
				end

				optimize!(model)
				status = string(termination_status(model))
				if (status != "INFEASIBLE") & (status != "LOCALLY_INFEASIBLE")
					break
				end
			end
		
		println("status : ", status)
		write(Output," status : ", string(status)) 
		write(Output,";  time CPU :", string(solve_time(model))) 
		write(Output,"\n  dt : ",string(dt))
		write(Output,";  alpha : ",string(alphaHalf*180/pi))
		write(Output,";  theta_min : ",string(theta_min*180/pi))
		if string(status) == "LOCALLY_SOLVED" || string(status) == "OPTIMAL"
			write(Output,"\n NUM SATELLITES : ", string(NUM_SATELLITE))
			write(Output,"\n  demi-axe : ", string((value.(model[:altitude]))))
			write(Output,"\n  sin_inclinaison : ",string(value.(model[:sin_inclinaison])))
			write(Output,"\n  sin_noeudAscendant : ",string(value.(model[:sin_noeudAscendant])))
			write(Output,"\n  meanAnomaly : ",string(value.(model[:meanAnomaly])))
			write(Output,"\n  theta_max_i : ",string((value.(model[:Thetamax])*180)/pi))
			write(Output,"\n"," ")
		end
		write(Output,"\n"," ")
	end			
end

		
			
		
