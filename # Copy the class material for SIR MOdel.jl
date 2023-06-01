# Copy the class material for SIR MOdel


# tips: function: return -priority
#f(x,y)=x+y
#f(2,3)
#1. Initialize
function updateSIR(popnvector)
    susceptibles = popnvector[1];
    infecteds    = popnvector[2]; 
    removeds     = popnvector[3];
    newS= susceptibles-lambda*susceptibles*infecteds*dt;
    newI = infecteds + lambda*susceptibles*infecteds*dt - gam*infecteds*dt  
    newR = removeds + gam*infecteds*dt
    return [newS newI newR]
end

# likewise, a run of the model uses exactly the same code ... but we'll play a bit with the values that determine a run

gam = 1/20.         # recovery rate parameter  (ditto)
lambda = 0.0005     # infection rate parameter
dt = 0.5            # length of time step in days
tfinal = 610.;      # respecting community values: lowercase only in the names 
s0 = 2000.          # initial susceptibles, note that we use the  type Float64 from the start
i0 = 4.             # initial infecteds; set this to 1. to  mimic an epidemic with an index case
r0 = 0.             # not always the case, of course

# Initialize 
nsteps=round(Int64,tfinal/dt)
resultvals=Array{Float64}(undef,nsteps+1,3)  
timevec=Array{Float64}(undef,nsteps+1)
resultvals[1,:]=[s0,i0,r0]
timevec[1]=0

for step=1:nsteps
    resultvals[step+1,:]=updateSIR(resultvals[step,:])
    timevec[step+1]=timevec[step]+dt
end
svals = resultvals[:,1];  # get the results
ivals = resultvals[:,2];

plot(svals, ivals, 
    title = "First look at I vs S plot",
    xlabel = "Susceptibles",
    ylabel = "Infecteds")        # and plot them

# we first vary S(0); the initial threshold is (0.1)/(1/20000) = 2000
# one can also vary gam and lambda and dt



