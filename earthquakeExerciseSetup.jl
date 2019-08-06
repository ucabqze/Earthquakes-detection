using Random
using PyPlot
using DelimitedFiles

function earthquake_exercise_setup()
    # explosion detector (using spiral coordinate system)
    # define the coordinate system:
    S=2000 # number of points on the spiral
    rate=25 # angular rate of spiral
    sd=0.2 # standard deviation of the sensor Gaussian noise

    #Establishment of coordinates for points on the spiral
    x=zeros(S); y=zeros(S)
    for s=1:S
        theta=rate*2*pi*s/S;  r=s/S
        x[s]=r*cos(theta); y[s]=r*sin(theta)
    end

    # define the locations of the detection stations on the surface
    # Also define what value on each sensor would be generated by explostions at s1 and s2
    N=30 # number of stations
    x_sensor=zeros(N); y_sensor=zeros(N)
    v=zeros(S,S,N)
    for sensor=1:N
        theta_sensor=2*pi*sensor/N
        x_sensor[sensor]=cos(theta_sensor); y_sensor[sensor]=sin(theta_sensor)
        for s1=1:S
            for s2=1:S
                v[s1,s2,sensor]=value(x[s1],y[s1],x[s2], y[s2], x_sensor[sensor],y_sensor[sensor])
            end
        end
    end

    #the observed sensor values
    val = readdlm("EarthquakeExerciseData.txt", '\t', Float64, '\n') #using DelimitedFiles

    # Perform inference p(location|observed sensor values) given these sensor value
    logp=zeros(S,S)
    for s1=1:S
        for s2 = 1:S
            for sensor=1:N
                logp[s1,s2] += -0.5*(val[sensor]-v[s1,s2,sensor])^2/(sd^2) # Gaussian distribution
            end
        end
    end

    #get p(s1,s2│v)
    max_logp = maximum(logp)*ones(S,S)
    p = exp.(logp-max_logp) #do exponential
    p = p/sum(p) # normalise

    #get p(s1│v)
    p_s1 = zeros(S)
    for s1=1:S
        for s2 = 1:S
            p_s1[s1] += p[s1,s2]
        end
    end

    p_s1=p_s1/sum(p_s1) # normalise


    # plot the posterior and most likely location of the explosion:
    maxp,maxind = findmax(p_s1) #find posterior and most likely location for s1

    maxp_s2,maxind_s2 = findmax(p[maxind,:]) #find the most likely location for s2

    #val_clean
    val_clean = val .- sd*randn()

    figure()
    for s=1:S
        plot(x[s],y[s],".",color=(1-(p_s1[s]/maxp))*[1,1,1])
    end
    for theta=0:0.01:2*pi
       plot(cos(theta),sin(theta),".",color=[0,0,0])
    end
        sf=0.2
    for sensor=1:N
        plot(x_sensor[sensor],y_sensor[sensor],"o",color=[1,0,0])
        theta_sensor=2*pi*sensor/N
        x_sensor_tm=(1+sf*val_clean[sensor])*cos(theta_sensor)
        y_sensor_tm=(1+sf*val_clean[sensor])*sin(theta_sensor)
        plot([x_sensor[sensor], x_sensor_tm],[y_sensor[sensor], y_sensor_tm],"-r")
        x_sensor_nm=(1+sf*val[sensor])*cos(theta_sensor+0.05)
        y_sensor_nm=(1+sf*val[sensor])*sin(theta_sensor+0.05)
        plot([x_sensor[sensor], x_sensor_nm],[y_sensor[sensor], y_sensor_nm],"-m")
    end

    plot(x[maxind],y[maxind],"m+",markersize=20,label="estimated (most likely s1)")
    plot(x[maxind_s2], y[maxind_s2],"rx",markersize=20,label="estimated (most likely s2)")
    print("the posterior p(s1|v) is ", maxp, ", and the most likely location is (", x[maxind],", ",y[maxind],")")
    legend()

end



function value(x_1,y_1,x_2, y_2, x_sensor,y_sensor)
    d_1_square = (x_1-x_sensor)^2 + (y_1-y_sensor)^2
    d_2_square = (x_2-x_sensor)^2 + (y_2-y_sensor)^2
   return 1/(0.1+ d_1_square) + 1/(0.1+ d_2_square)#+sd*randn() #noise is drawn from a zero mean
end
