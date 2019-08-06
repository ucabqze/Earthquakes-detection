using Random
using PyPlot
using DelimitedFiles

function oneOrTwo_122()
    S=2000 # number of points on the spiral
    rate=25 # angular rate of spiral
    sd=0.2 # standard deviation of the sensor Gaussian noise

    #Establishment of coordinates for points on the spiral
    x=zeros(S); y=zeros(S)
    for s=1:S
        theta=rate*2*pi*s/S;  r=s/S
        x[s]=r*cos(theta); y[s]=r*sin(theta)
    end

    #the observed sensor values
    val = readdlm("EarthquakeExerciseData.txt", '\t', Float64, '\n') #using DelimitedFiles

    #Calcualte the case of two explosions
    N=30
    x_sensor_2=zeros(N); y_sensor_2=zeros(N)
    v_2=zeros(S,S,N)
    for sensor=1:N
        theta_sensor=2*pi*sensor/N
        x_sensor_2[sensor]=cos(theta_sensor); y_sensor_2[sensor]=sin(theta_sensor)
        for s1=1:S
            for s2=1:S
                v_2[s1,s2,sensor]=value_H2(x[s1],y[s1],x[s2], y[s2], x_sensor_2[sensor],y_sensor_2[sensor])
            end
        end
    end

    # logp(vi│s1,s2,H2)
    logpv_s1s2H2=zeros(S,S,N)
    for sensor=1:N
        for s1=1:S
            for s2=1:S
                logpv_s1s2H2[s1,s2,sensor] = -0.5*(val[sensor]-v_2[s1,s2,sensor])^2/(sd^2) # Gaussian distribution
            end
        end
    end

    #p(v_i│s1,s2,H2)
    pv_s1s2H2 = exp.(logpv_s1s2H2)###

    #logp(vi│H2)
    pv_H2 = zeros(N)###
    for sensor=1:N
        for s1=1:S
            for s2=1:S
                pv_H2[sensor] += pv_s1s2H2[s1,s2,sensor]/(S*S)###
            end
        end
    end

    logpv_H2 = log.(pv_H2)###


    #Calcualte the case of one explosions
    x_sensor_1=zeros(N); y_sensor_1=zeros(N)
    v_1=zeros(S,N)
    for sensor=1:N
        theta_sensor=2*pi*sensor/N
        x_sensor_1[sensor]=cos(theta_sensor); y_sensor_1[sensor]=sin(theta_sensor)
        for s=1:S
            v_1[s,sensor]=value_H1(x[s],y[s],x_sensor_1[sensor],y_sensor_1[sensor]) # explosion value
        end
    end

    #logp(vi│s,H1)
    logpv_sH1=zeros(S,N)
    for sensor=1:N
        for s=1:S
            logpv_sH1[s,sensor] = -0.5*(val[sensor]-v_1[s,sensor])^2/(sd^2) # Gaussian distribution
        end
    end

    #p(v_i│s,H1)
    pv_sH1 = exp.(logpv_sH1)###

    #logp(vi│H1)
    pv_H1 = zeros(N)###
    for sensor=1:N
        for s=1:S
            pv_H1[sensor] += pv_sH1[s,sensor]/S ###
        end
    end

    logpv_H1 = log.(pv_H1)###

    #logp(v│H2) -logp(v│H1)
    dif = sum(logpv_H2)-sum(logpv_H1)

    print("log p(v|H2) − log p(v|H1) is: ")
    print(dif)



end
#Function to calculate value
function value_H1(x_true,y_true,x_sensor,y_sensor)
   return 1/(0.1+ (x_true-x_sensor)^2 + (y_true-y_sensor)^2)
end

function value_H2(x_1,y_1,x_2, y_2, x_sensor,y_sensor)
    d_1_square = (x_1-x_sensor)^2 + (y_1-y_sensor)^2
    d_2_square = (x_2-x_sensor)^2 + (y_2-y_sensor)^2
   return 1/(0.1+ d_1_square) + 1/(0.1+ d_2_square)#+sd*randn() #noise is drawn from a zero mean
end
