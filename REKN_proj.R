
#Geolocation analysis with open source tools 2016: NAOC, WA DC
#E. Bridge, M.T. Hallworth, S. Livoski, and E. Rakhimberdiev

# Check to make sure the required packages are installed on your machine
# If not, they will be installed

reqPackages <- c("devtools","digest","GeoLight","geosphere","raster","fields","forecast",
                 "circular","truncnorm","parallel","bit","rgdal","CircStats","Rcpp",
                 "RcppArmadillo","ggmap","ggsn","sp","maptools","rgeos","MASS")


get.packages <- reqPackages[!(reqPackages %in% installed.packages()[,"Package"])]

if(length(get.packages)>0){
  install.packages(get.packages,repos = "https://cloud.r-project.org/", dependencies = TRUE)
}

# Install necessary packages from Github using the devtools library #

library(devtools)
install_github("SWotherspoon/SGAT")
install_github("SLisovski/TwGeos") 
install_github("SLisovski/GeoLight", ref = "Update_2.01", force = T)
install_github("eldarrak/FLightR")

library(GeoLight)
library(TwGeos)
library(SGAT)
library(FLightR)
library(MASS)  #needed for fitting distributions
library(maptools)


#load the TwGeos package
library("TwGeos")                                 


##wd is: T:\redknot_geo

#For BAS geos: read the data into a dataframe called d.lig
d.lig <- readLig("data/TRES.lig", skip = 0)

#reduce the dataframe to just Date and Light
d.lig <- subset(d.lig,select=c("Date","Light"))   

#For Intigeo geos: read the data into a dataframe called d.lux
d.lux <- readMTlux("data/Godwit.lux")     
#log transformation to reduce range of light values (max value = 70,000)
d.lux$Light<- log(d.lux$Light)
head(d.lux)

str(d.lux)

#Plot the light data. 
#Each day is represented by a thin horizontal line that plots the light values as 
#grayscale (dark = low light, white = max light)in order from bottom to top. 
#A light image allows you to visualize an entire dataset at once, and easily spot discrepancies 
#in light to dark transitions. Notive dark areas either in the beginning (pre-deployment) and 
#the end (post-retrieval)

lightImage( tagdata = d.lig, # light data
            offset = 18,     # adjusts the y-axis to put night (dark shades) in the middle
            zlim = c(0, 64), # y axis
            dt = 120)        # miniumn dark period

#Defining Transitions by setting a threshold (sets time of sunrise and sunset)
#Choose the lowest threshold value that is consistently above any noise in the nighttime
#light levels. For many .lig datasets 2.5 is above any nighttime noise. 

threshold <- 2.5

library(SGAT)
?preprocessLight()

# For PC 
twl <- preprocessLight(tagdata = d.lig,
                       threshold = threshold,
                       offset = 18, 
                       lmax = 12, 
                       gr.Device = "default")

#Another (better?) way to find the twilights without using the interactive process is the findTwilights
#function in TWGeos package. This function finds the twilight times without individual inspection and
#without editing. 
#First we need to know one data and time when it is night within the dataset. This date and time 
#where night is known is the seed. You can specify the seed in as.POSIXct format or you can specify
#the seed date by clicking on the light image. 

#Specify seed manually
seed <- as.POSIXct("2011-11-01 04:00", origin  = "1970-01-01", tz = "GMT")

twl  <- findTwilights(tagdata = d.lig, 
                      threshold = threshold, 
                      include = seed,
                      dark.min = 0) # minimum dark period in minutes

#Specify seed interactively
# Plot the data 
plot(d.lig$Date[3000:5000], d.lig$Light[3000:5000], type = "o", pch = 16, cex = 0.5)

seed <- as.POSIXct(locator(n=1)$x, origin  = "1970-01-01", tz = "GMT") # click at any time during the night

twl  <- findTwilights(tagdata = d.lig, 
                      threshold = threshold, 
                      include = seed)
#If the twl object is empty check the following things: 1) If you used locator() function make sure you 
#clicked on a date within the time frame that you may have subsetted above. 2) Double check that the light 
#levels fall below the threshold level.

#Plot twilights overlaid on the light image
lightImage(tagdata = d.lig, 
           offset = 18, 
           zlim = c(0, 12),
           dt = 120)
tsimagePoints(twl$Twilight, 
              offset = 12, 
              pch = 16, 
              cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

#Removing/Fixing Twilight Outliers
#Rule: if a twilight time is outlier.mins (e.g. 45) minutes different to its surrounding twilight times, 
#defined by window - the number of surrounding twilight times (e.g. 4), and the surrounding twiligth times 
#are within stationary.mins (e.g. 15) minutes, the outlayer will be moved in between the two neighboring twilights 
#or deleted if the sourrounding twilights are > stationary.mins.

twl <- twilightEdit(twilights = twl, 
                    window = 4,           # two days before and two days after
                    outlier.mins = 45,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)

#This function provides a dataframe with both the edited twlight times (in column “Twilight”) and the original 
#twilight times (in column “Twilight0”). Now that we have adjusted the twilights - we have new information in out twl file.

#It’s a good idea to save this dataframe for future reference so that the edits made are documented and your analysis repeatable
write.csv(twl, file = "Edited_twilights", quote = FALSE, row.names = FALSE)

#In order to take advantage of the edits - you can subset the data using when Deleted == FALSE this removes 
#any of the twilights which were deleted.
twl <- subset(twl, Deleted == FALSE)

#Check the sampling interval of the geolocator
head(d.lig)
d.lig$Date[3]-d.lig$Date[2]

#Here we adjust the twilight time back to account for the max light level readings every # minutes
twl <- twilightAdjust(twilights = twl, 
                      interval = 120) # The unit here is seconds

#Truncate data to only include deployment period
On_Bird <- as.POSIXct(c("2011-06-26", "2012-05-11"), tz = "GMT")

#Make a light image and use abline() to draw some orange lines that enclose the time period when the geolocator was on the bird.

lightImage(d.lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

abline(v  = On_Bird, col  ="orange", lwd = 3)

#You could also specify the dates when the tag was on the bird by using the locator() function and 
#clicking on the light image to specify the dates.

lightImage(d.lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

On_Bird <- as.POSIXct(locator(n=2)$x, origin = "1970-01-01", tz = "GMT")

#The following code subsets the twilight data to include only those data between the capture dates.
twl <- twl[twl$Twilight > On_Bird[1] & twl$Twilight < On_Bird[2],]
##plot just the on-bird twilights on top of a light image.
lightImage(d.lig, offset = 19)
tsimagePoints(twl$Twilight, offset = 19, pch = 16, cex = 0.5,
              col = ifelse(twl$Rise, "dodgerblue", "firebrick"))

#SGAT process

#Calibration
#To perform calibration we must have data from a known location (usually the tag deployment or animal recapture site). 
#It is best if the tag is on the animal during the calibration period. Let’s create an object to store the capture location.
# Calibration Coordinates

CapLocs <- cbind(-80.46,42.62) # Longitude, latitude

#Make a subset of data that includes known calibration period. 

cal <- subset(twl, twl$Twilight < as.POSIXct("2011-07-06", "GMT"))

#The alternative is to use the the light-level information direclty off the geolocator to determine when the animal was 
#at the capture location. In order to use the second approach we need to know what the actual twilight events are and how 
#they match up to the observed twilights derived from the geolocator data.

Times <- seq(from = d.lig$Date[1], 
             to = d.lig$Date[length(d.lig$Date)], 
             by = "day")

rise <- rep(c(TRUE, FALSE), length(Times))


# making predicted twilight times given location and zenith #

KnownTwl <- data.frame(Twilight = twilight(rep(Times, each = 2),
                                           lon = CapLocs[1], 
                                           lat = CapLocs[2], 
                                           rise = rise,
                                           zenith = 94),
                       Rise = rise) 

#Once we have the known twilight times from the capture location you can plot the light data from the geolocator (d.lig) 
#and overlay the known twilight times from the capture location. The sun rise at the known location is shown in blue and the 
#sunset is shown in red. When the light data from the tag corresponds with the known twilights at the capture location the 
#animal was likely stationary at/around the capture site. When the measured light levels deviate from the known twilights 
#the bird was somewhere else. You can use the locator() function to click on the image to determine dates, or you can create 
#an R object that specifies the dates.

lightImage(d.lig,
           offset = 19, 
           zlim = c(0,64), 
           main = "Light Image") 

tsimagePoints(KnownTwl$Twilight, 
              offset = 19, 
              pch = 16, cex = 0.5, 
              col = ifelse(KnownTwl$Rise, "blue", "red"))

# adds line at two equinoxes for reference. Change the dates if necessary (can vary by year) #
eqnx<-as.POSIXct(c("2011-09-23", "2012-03-20"), tz = "GMT") 
abline(v = eqnx, lwd=3, lty=3, col="purple")

#Look for period when deviations between measured light levels from known twilights at capture location

calibrationPeriod <- subset(twl, Twilight0 < as.POSIXct("2011-07-15", "GMT"))
str(calibrationPeriod)

#Get calibration data into GeoLight format and calculate Sun-elevation level
##change data format
cal2 <- data.frame(tFirst = calibrationPeriod$Twilight[1:(nrow(calibrationPeriod)-1)],
                   tSecond =  calibrationPeriod$Twilight[2:nrow( calibrationPeriod)],
                   type = abs( calibrationPeriod$Rise[1:(nrow(calibrationPeriod)-1)]-2))

##apply getElevation"
sunAngle <- GeoLight::getElevation(twl = cal2,
                                   known.coord = CapLocs,
                                   plot = FALSE, 
                                   lnorm.pars = FALSE)

##sunAngle is the angle of the sun relative to the horizon where the horizon is zero. 
##Positive values are above the horizon and negative values are below the horizon.
sunAngle

#Determine Zenith angle

library(MASS)

# Calculate solar time from calibration data 
sun <- solar(calibrationPeriod[,1])

# Adjust the solar zenith angle for atmospheric refraction
zenithAngle <- refracted( zenith(sun = sun,
                                 lon = CapLocs[1], 
                                 lat = CapLocs[2]))

twilight_time <- twilight(tm = calibrationPeriod[,1],
                          lon = CapLocs[1], 
                          lat = CapLocs[2], 
                          rise = calibrationPeriod[,2],
                          zenith = quantile(zenithAngle,probs=0.95))

# Determine the difference in minutes from when the sun rose and the geolocator said it rose 
twl_deviation <- ifelse(calibrationPeriod$Rise, 
                        as.numeric(difftime(calibrationPeriod[,1], twilight_time, units = "mins")),
                        as.numeric(difftime(twilight_time, calibrationPeriod[,1], units = "mins")))

# Describe the distribution of the error 
twl.dist <- fitdistr(abs(twl_deviation), "log-Normal")

# save the Twilight model parameters
alpha <- c(twl.dist$estimate[1], twl.dist$estimate[2]) 

# Make some plots to visualize the data 
par(mfrow=c(1,2))  
hist(abs(twl_deviation), freq = F,
     yaxt="n",
     xlim = c(0, 60),
     breaks=10,
     col="gray",
     main = "",
     xlab = "Twilight error (mins)")
axis(2,las=2)
lines(seq(0,60, length = 100),
      dlnorm(seq(0,60, length = 100), alpha[1], alpha[2]), col ="red",lwd = 3, lty = 2)

#Zenith angle plot
par(bty="l")
plot(median(zenithAngle,na.rm=TRUE),ylim = c(80,120),pch=19,cex=1.25,ylab="Zenith Angle")
segments(1,quantile(zenithAngle,probs=0.025),1,quantile(zenithAngle,probs=0.975))
points(1,max(zenithAngle,na.rm=TRUE),col="red",pch=20)


#Assign the zenith that we will use in the analyses. If using a simple threshold approach (GeoLight) 
#the median zenith Angle is probably good enough.

Zenith <- quantile(zenithAngle, probs = 0.95)

#Initial Path
#Look at a subset o the data which includes twilights that 1) occurred after the first calibrartion date 
#(presumably the deployment date) and 2) were not deleted according to the output of twilightEdit().
#Important - the tol setting in the thresholdPath model sets the tolerance around equinox (i.e. filter out those points). 
#Latitudinal estimates are unrelibale around the equinox periods because the change in day length is similar everywhere. 
#tol values = 0 indicates save all points - larger tol values (0.2) filter out quite a few points. 
#For this example, I set the tolerance to 0 (not the default) so I can look at all the data - even around equinox periods.

#Geolight has simple functions for generating coordinates from twlight data (coord()) and for plotting location data (tripMap()). 
#But first you have to get the data into the tFirst-tSecond format.

library(GeoLight)

##change data format
twlgl <- data.frame(tFirst = twl$Twilight[1:(nrow(twl)-1)],
                    tSecond =  twl$Twilight[2:nrow( twl)],
                    type = abs( twl$Rise[1:(nrow(twl)-1)]-2))

##generate coordinates

InitialPath <- coord(twl = twlgl, degElevation = sunAngle)

##plot the data on a map
tripMap(crds = InitialPath, equinox = TRUE, xlim = c(-90,-70), ylim = c(10,50), legend = TRUE)

#The SGAT package offers a slightly more sophisticated function to generate an initial path.

InitialPath <- thresholdPath(twilight = twl$Twilight,
                             rise = twl$Rise,
                             zenith = Zenith,
                             tol = 0)

#You can use tripmap to show the path, as we did above…

tripMap(crds = InitialPath$x,
        equinox = TRUE, 
        xlim = c(-90,-70),
        ylim = c(10,50), 
        legend = TRUE)

#…or here is some alternative code for making figures to show the inital path.

year<-c("2011-01-01","2012-01-01")

layout(matrix(c(1,3,
                2,3), 2, 2, byrow = TRUE))

par(mar=c(2,4,2,0))

plot(InitialPath$time, InitialPath$x[, 2], 
     type = "b", 
     pch = 16, 
     cex = 0.5, 
     ylab = "Latitude", 
     xlab = '', 
     xaxt="n",
     col=ifelse(InitialPath$time<as.POSIXct(year[2],format="%Y-%m-%d"),"blue","green"))

abline(h = CapLocs[2])
abline(v = as.POSIXct("2011-09-23"),col="red",lty=2,lwd=1.5)
abline(v = as.POSIXct("2012-03-21"),col="red",lty=2,lwd=1.5)

par(mar=c(2,4,2,0))
plot(InitialPath$time, InitialPath$x[, 1],
     type = "b",
     pch = 16,
     cex = 0.5,
     ylab = "Longitude",
     xlab = '',
     col=ifelse(InitialPath$time<as.POSIXct(year[2],format="%Y-%m-%d"),"blue","green"))

abline(h = CapLocs[1])
abline(v = as.POSIXct("2011-09-23"),col="red",lty=2,lwd=1.5)
abline(v = as.POSIXct("2012-03-20"),col="red",lty=2,lwd=1.5)

library(maptools)

#Checking availability of rgeos

data("wrld_simpl")
plot(wrld_simpl,
     col = "grey95",
     xlim = c(-120,-60),
     ylim=c(0,40))
box(which="plot")
lines(InitialPath$x,
      col = ifelse(InitialPath$time<as.POSIXct(year[2],format="%Y-%m-%d"),"blue","green"))

points(InitialPath$x, 
       pch = 16, 
       cex = 0.5, 
       col=ifelse(InitialPath$time<as.POSIXct(year[2],format="%Y-%m-%d"),"blue","green"))
points(CapLocs[1],CapLocs[2],pch=19,cex=2)


#Creating path with SGAT

#The SGAT package contains functions that let you apply a movement model, landcover masking, and estimates 
#of light detection error to generate a Markov Chain Monte Carlo algorithm that generates many likely movement paths 
#from an individual data set. To run this sort of analysis we have to provide several parameters based on what we know 
#about our data and the animal we are studying. We also have to provide an initial estimate of the movement path, 
#which in this example will be the simple threshold path we generated above.

#The workflow for creating the final path is:
#1.Give the model an inital path - just created above
#2.Define the mid-points between locations (needed to generate path)
#3.Define a movement model (usually a distribution of flight speeds or daily flight distances wherein long flights are rare and short flights are common.)
#4.Establish the error distribution for twilight detection based on your calibration data
#5.(Optional) Establish a mask or grid that delineates location that are extremely unlikely to be part of the movement path (e.g. small landbirds will not reside on large water bodies).
#6.Run model multiple times until model converges and throw out the first interations as burn-in
#7.Refine the model further from previous runs
#8.Run the model several hundred times to create the final path

#First, we rename the initial path from our results above (the mapped route), and then find the midpoints 
#between each pair of consecutive locations.

#Initialize SGAT Estella model

#Next, we create an object called fixedx that specifies the known locations for some or all of the calibration 
#data and categorizes those locations as fixed so that the model does not try to estimate them.

# Specify the dates when locations are fixed # 

fixedx <- rep(FALSE, nrow(x0))

fixedx[c(1:10,(nrow(x0)-3):nrow(x0))] <- TRUE

x0[fixedx, 1] <- CapLocs[1]
x0[fixedx, 2] <- CapLocs[2]

z0 <- trackMidpts(x0) # update z0 positions

#Establish a movement model

#We need to provide a mean and standard deviation for a gamma distribution of flight speeds that get applied to 
#each day of the analysis period. We typically want short (near zero) distance flights to be common and long distance 
#flights to be relatively rare. So both mean and distribution should be small.

# Flight model
beta <- c(0.7, 0.08) #mean and sd 

#Restrict path to land (Not suitable for Red Knots!)

#You can restrict the path to locations on land - birds are still able to fly over water but have stationary locations on land. 
#This makes sense for a terrestrial bird like warblers. It may not for other species, as always, consider the ecology of the 
#species when conducting the analysis. Here, we will restrict the paths to the wrld_simpl.

#We include a prior distribution so only locations on land within the Americas are used. We first create a function to covert 
#the wrld_smpl shapefile to a binary surface.

## Function to construct a land/sea mask
distribution.mask <- function(xlim, ylim, n = 4, land = TRUE, shape) {
  
  r <- raster(nrows = n * diff(ylim),
              ncols = n * diff(xlim),
              xmn = xlim[1], 
              xmx = xlim[2], 
              ymn = ylim[1], 
              ymx = ylim[2], 
              crs = proj4string(shape))
  
  r <- cover(rasterize(shape, shift = c(-360, 0), r,1, silent = TRUE), 
             rasterize(shape, r, 1, silent = TRUE), rasterize(shape, r, 1, silent = TRUE))
  
  r <- as.matrix(is.na(r))[nrow(r):1, ]
  
  if (land) 
    r <- !r
  xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
  ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)
  
  function(p) {
    r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
  }
}

#Set the x and y limits for the raster surface. It’s important that you set these limits to include all longitude and latitude 
#values where the individual may occur throughout the year. If your xlim and ylim don’t include all areas you won’t get results 
#for that portion of the range.

xlim <- c(-120,-60) 
ylim <- c(0,70)

#Here we define the distribution mask to only include stationary locations that occur on land. This makes sense for something 
#like a Tree Swallow but may not for other species - always consider the ecology of the species while conducting the analysis.

## Define the distribution mask
is.dist <- distribution.mask(shape=wrld_simpl,
                             xlim = xlim,
                             ylim = ylim,
                             n = 4,
                             land = TRUE)

# Define the log prior for x and z
log.prior <- function(p) {
  f <- is.dist(p)
  ifelse(f | is.na(f), 0, -10)
}

#The log prior function produces a probability value for locations on land and locations at sea. If you feed the function a 
#set of coordinates that is at sea the prior probilbity is -10

## try a location in the Gulf of Mexico
log.prior(cbind(-90, 25))

#If you feed the function a set of coordinates that is at sea the prior probilbity is 0

## try a location on land
log.prior(cbind(-80, 45))

#In the SGAT package we need to specify a model for the analysis. Below we specify a few key parameters.
#1.twilight = twilight times that we determined above.
#2.rise = a logical vector sunrise = TRUE - this is calculated at the same time when you define twilights.
#3.twilight.model = the distribution type for the difference between observed twilight and expected twilight.
#4.alpha = the shape of the twilight.model distribution
#5.beta = the movement model parameter
#6.logp.x and logp.z = constraints set on the x and z (intermediate) positions. This is where you set the constraints for land
#7.x0 = initial values for the birds path (x positions)
#8.z0 = initial values for the birds path (z positions)
#9.zenith = the zenith angle to be used. This can take a single value (no change in zenith throughout the year) or a vector of nrow(twl) if you want to use different zenith angles.
#10.fixedx = a vector telling the model which locations need to be estimated because positions are unknown.

# Define the threshold model - slimilar to above #

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "ModifiedLogNormal",
                        alpha = alpha,
                        beta = beta,
                        # Here is where we set the constraints for land
                        logp.x = log.prior, logp.z = log.prior, 
                        x0 = x0,
                        z0 = z0,
                        zenith = Zenith,
                        fixedx = fixedx)

#We also need to define the error distribution around each location. We set that using a multivariate normal distribution

# This defines the error distribution around each location #

proposal.x <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0))

#We then fit the model using the estelleMetropolis sampler. Here you can set the number of iterations (iters), 
#the thinning rate (thin) and the number of chains to run (chains).

fit <- estelleMetropolis(model = model,
                         proposal.x = proposal.x,
                         proposal.z = proposal.z,
                         iters = 1000, # This value sets the number of iterations to run
                         thin = 10,
                         chains = 3)

### Fine Tuning 

proposal.x <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(chainLast(fit$x)))
proposal.z <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(chainLast(fit$z)))

# Summarize the results of the previous chains to initialize the next model. 
# This makes full use of all the data

xsum <- locationSummary(fit$x,
                        time = fit$model$time,
                        collapse = TRUE) 

zsum <- locationSummary(fit$z,
                        time = fit$model$time,
                        collapse = TRUE)

fit <- estelleMetropolis(model = model,
                         proposal.x = proposal.x,
                         proposal.z = proposal.z,
                         # initialize the model using the median path from pervious runs
                         x0 = cbind(xsum$'Lon.50%',xsum$'Lat.50%'), 
                         z0 = cbind(zsum$'Lon.50%',zsum$'Lat.50%'),
                         iters = 1000, # This value sets the number of iterations to run
                         thin = 10,
                         chains = 3)

# Final Run

proposal.x <- mvnorm(chainCov(fit$x),s=0.1)
proposal.z <- mvnorm(chainCov(fit$z),s=0.1)

# Note the increase in number of interations - this takes a bit longer to run
xsum <- locationSummary(fit$x,
                        time = fit$model$time,
                        collapse = TRUE) 

zsum <- locationSummary(fit$z,
                        time = fit$model$time,
                        collapse = TRUE)

fit <- estelleMetropolis(model = model,
                         proposal.x = proposal.x,
                         proposal.z = proposal.z,
                         # initialize the model using the median path from pervious runs
                         x0 = cbind(xsum$'Lon.50%',xsum$'Lat.50%'),
                         z0 = cbind(zsum$'Lon.50%',zsum$'Lat.50%'),
                         iters = 5000, # This value sets the number of iterations to run
                         thin = 10,
                         chains = 3)

#Initial results
#To view the results of the MCMC process we first need to create an empty raster to fill with the results.

# This step makes an empty raster #
r <- raster(nrows=4*diff(ylim), # this sets the spatial resolution 
            ncols=4*diff(xlim), # ditto
            xmn=xlim[1],
            xmx=xlim[2],
            ymn=ylim[1],
            ymx=ylim[2])

#We can then summarize the data using the slices function. Here we “slice” the intermediate locations by day. 
#You can set breaks to day, month, etc.

S <- slices(type="intermediate",
            breaks="day",
            mcmc=fit,
            grid=r)

#Next if we are interested in certain days or a period of interest you may want to pull out locations that fall within a 
#certain time period. Here we extract the dates from the MCMC process.

DATES <- S$mcmc[[1]]$time[ which( S$mcmc[[1]]$rise==TRUE) ]

#We can create objects that correspond with certain dates of interest. These objects are the rows that correpsond 
#with the dates of interest.

# Find the rows that correspond with the dates of interest 

ReleaseDay <- 1
Aug01 <- which(strptime(DATES, format = "%Y-%m-%d", tz = "GMT")==as.POSIXct("2011-08-01",format="%Y-%m-%d",tz="GMT"))
Jan01 <- which(strptime(DATES, format = "%Y-%m-%d", tz = "GMT")==as.POSIXct("2012-01-01",format="%Y-%m-%d",tz="GMT"))
Mar01 <- which(strptime(DATES, format = "%Y-%m-%d", tz = "GMT")==as.POSIXct("2012-03-01",format="%Y-%m-%d",tz="GMT"))

#We can then use the sliceInterval function to extract and summarize the locations between the dates of interest. 
#For example, if we were interested in time of year when an individual is on the breeding grounds - we may specify the 
#dates from the ReleaseDay until the end of the breeding season - say Aug01. We may also be interested in the stationary 
#non-breeding period - say between Jan01 and Mar01. Note these dates were chosen haphazardly for this analysis. 
#When conducting your own analysis make sure to consider the ecology of your species.

# "Slice" the data and save all dates between Release date and August 1.
breed <-slice(S,k=c(ReleaseDay:Aug01)) # k = c(Start Date : End Date) using row numbers
nonbreed <- slice(S, k=c(Jan01:Mar01))

#We can plot the results on the same map using the following code

plot(wrld_simpl, 
     xlim = c(-96, -72),
     ylim = c(23,47),
     border = "gray",
     col = "gray88")

plot(breed,
     useRaster=TRUE,
     axes=FALSE, 
     add=TRUE,
     legend=FALSE,
     col=c("transparent",rev(bpy.colors(50))),
     cex.axis=0.7)

plot(nonbreed,
     useRaster=TRUE,
     axes=FALSE, 
     add=TRUE,
     legend=FALSE,
     col=c(rep("transparent",5),rev(bpy.colors(50))),
     cex.axis=0.7)
plot(wrld_simpl,border="gray",add=TRUE)
raster::scalebar(d = 500, xy = c(-96,25),divs=2,lonlat=T,label = c(0,250,500),below="km",type="bar",cex=0.5)
box()

#Getting the summary statistics is fairly easy. Here we collapse all chains using the locationSummary function.

zummary <- locationSummary(fit$z,
                           time=fit$model$time,
                           collapse=TRUE) 

head(zummary)

#Make a plot that shows the median Longitude and Latitude through time along with the 95% credible interval.

par(mfrow=c(2,1),mar=c(4,4,0,0))
plot(zummary$Time1,zummary$"Lon.50%",
     type="l",
     ylab = "Longitude",
     xlab = "",
     yaxt = "n")
axis(2, las = 2)
polygon(x=c(zummary$Time1,rev(zummary$Time1)),
        y=c(zummary$`Lon.2.5%`,rev(zummary$`Lon.97.5%`)),
        border="gray",
        col="gray")
lines(zummary$Time1,zummary$"Lon.50%")
plot(zummary$Time1,zummary$"Lat.50%",
     type="l",
     ylab = "Latitude",
     xlab = "",
     yaxt = "n",
     ylim = c(20,60))
axis(2, las = 2)
polygon(x=c(zummary$Time1,rev(zummary$Time1)),
        y=c(zummary$`Lat.2.5%`,rev(zummary$`Lat.97.5%`)),
        border="gray",
        col="gray")
lines(zummary$Time1,zummary$"Lat.50%")

#Migratory Routes
#Below are the approximate dates at each ‘stationary’ period determined via changes in the sunrise / sunset times. 
#The summary map shows approximate location. These are just to get a general idea of when and where they stopped.

library(GeoLight)

twl.geolight <- data.frame(Twilight=twilight(fit$model$time[-length(fit$model$time)],
                                             lon= zummary$`Lon.50%`,
                                             lat=zummary$`Lat.50%`,
                                             fit$model$rise[-length(fit$model$time)],
                                             zenith=Zenith,
                                             iters=3),
                           Rise=fit$model$rise[-length(fit$model$time)])

twl.GL <- data.frame(tFirst=twl.geolight[-nrow(twl.geolight),1],
                     tSecond=twl.geolight[-1,1],
                     type=ifelse(twl.geolight[,2],1,2)[-nrow(twl.geolight)])



stops <- changeLight(twl = twl.GL[complete.cases(twl.GL),],
                     quantile = 0.9,
                     days = 2,
                     plot = FALSE, 
                     summary = TRUE)

Sites <- mergeSites(tFirst = twl.GL$tFirst,
                    tSecond = twl.GL$tSecond, 
                    type = twl.GL$type,
                    site = stops$site,
                    degElevation = 90-Zenith,
                    distThreshold = 300)

siteMap(cbind(zummary$`Lon.50%`, zummary$`Lat.50%`),
        map.range="America",
        site=Sites$site,
        xlim = xlim,
        ylim = ylim,
        type='cross',hull=F,legend=FALSE)

#End SGAT

#.......................................................

#FlightR
#For use with .lux (Intigeo) files

library(FLightR)

#Read in data
d.lux<-readMTlux("data/Godwit.lux")

#Plot data

lightImage( tagdata = d.lux, # light data
            offset = 12,     # adjusts the y-axis to put night (dark shades) in the middle
            zlim = c(0, 64), # y axis
            dt = 120)        # miniumn dark period

#Define twilights

seed <- as.POSIXct("2014-01-01 04:00", origin  = "1970-01-01", tz = "GMT")

twl.lux  <- findTwilights(tagdata = d.lux, 
                          threshold = 1.5, 
                          include = seed,
                          dark.min = 0) # minimum dark period in minutes

#Adjust/delete outlier twilights
twl.lux <- twilightEdit(twilights = twl.lux, 
                        window = 4,           # two days before and two days after
                        outlier.mins = 45,    # difference in mins
                        stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                        plot = TRUE)

#Save twilights
write.csv(twl.lux,"Edited_Twilights.csv", quote = FALSE, row.names = FALSE)

#Remove the transitions that were deleted by subsetting them out of the twl.lux data set.

twl.lux <- subset(twl.lux, Deleted == FALSE)

#Set the Twilight value used by default to be the edited twilights from above.

twl.lux$Twilight<-twl.lux$Twilight0

#Important Do not adjust time of twilights like we did above for SGAT if using FLightR because it does it automatically.
#Now we need to covert the data .lux file into a format that FLightR can use. Here we use the BAStag2TAGS function in the to do that.

lux.tags<-BAStag2TAGS(raw = d.lux, 
                      twl = twl.lux, 
                      threshold = 1.5)

#Now we write the lux.tags object to a file that is read in by FLightR to begin the analysis.
#NOTE Be sure to change the file path name to match your working directory.

write.csv(lux.tags,"data/Gotwit.lux.csv",quote=FALSE,row.names=FALSE)

# opens and formats data straight from either the file we just created or TAGS formatted data. 
Light.Data<-get.tags.data("data/Gotwit.lux.csv") 

#Set capture location

# start location longitude and latitude
CapLocs=c(5.43, 52.93) 

#Set calibration period
#The dawn (red) and dusk transitions (black) should be very similar (overlap). 
#Add two vertical lines to determine if 1) the end of the first calibration period and 2) the start
#of the second calibration period at the end the data before capture.

plot.slopes.by.location(Proc.data = Light.Data, 
                        location = CapLocs)

#Use abline to visualize potential calibration periods

# end of first calibration period
abline(v=as.POSIXct("2013-08-23"), col = "blue",lwd=2) 

# start of the second calibration period
abline(v=as.POSIXct("2014-05-05"), col = "green", lwd=2) 

#Create a data.frame that includes data from each calibration period.
#The columns are:
#1. start of the calibration period,
#1. end of the calibration period,
#1. longitude of the calibration period,
#1. latitude of the calibration period

#This will create two lines of data
Calibration.periods<-data.frame(calibration.start = as.POSIXct(c(NA, "2014-05-05")),   
                                calibration.stop = as.POSIXct(c("2013-08-23", NA)),
                                lon = CapLocs[1], 
                                lat = CapLocs[2]) 

#You can also use two geographic coordinates if you have more than one calibration location. 
#This may occur if you deploy the tag at one location and happen to capture it at a different 
#location. You could use the following code to use two different locations. lon = c(5.43, 6.00), 
#lat = c(52.93,52.94)

#Here is what the calibration period data frame looks like. 
#This object is telling FLightR to start calibration at the start of recording and stop on Aug 23
#2013 and the tag was located at the location specified by the latitude, longitude. 
#The second row tells FLightR that the second calibration period should start on May 05 2014 
#and continue until the end of the data set.

#View results
Calibration.periods

##   calibration.start calibration.stop  lon   lat
## 1              <NA>       2013-08-23 5.43 52.93
## 2        2014-05-05             <NA> 5.43 52.93

#Create a calibration object that FLightR will use to determine the relationship between 
#recorded light-levels and the expected light-levels.

#create a calibration object 
Calibration<-make.calibration(Proc.data = Light.Data, 
                              Calibration.periods = Calibration.periods,
                              model.ageing = FALSE,
                              plot.each = FALSE, 
                              plot.final = FALSE)

#save it for later use
save(Calibration, file = "data/FLightR_calibration")

#loads object called Calibration
load("data/FLightR_calibration") 

#Establishing spatial grid
#The default resolution is 50 x 50km grid cells. The inputs or terms for left, right, bottom and top 
#define your bounding box. The argument distance.from.land.allowed.to.use should be a vector with length of two. 
#The first number is a negative distance allowed to use while over land (restricts birds to flying only over 
#coastlines and water) and second is distance from land allowed to use while over water (restricts birds to flying 
#only over coastlines and land). The distance.from.land.allowed.to.stay should also be a vector of length two. 
#The first number is negative distance where the bird is allowed to be stationary (restricts birds to landing only 
#on coastlines and land). The second value is distance from land allowed to fly over during twilight while over water 
#(restricts birds to landing only on coastlines and water). Use infinity c(-Inf,Inf) to not use any restrictions. 
#We won't restrict paths for this example here.

Grid<-make.grid(left = -14, 
                bottom = 30,
                right = 13, 
                top = 57,
                #Use infinity to withold restrictions on migration paths
                distance.from.land.allowed.to.use = c(-Inf, Inf),  
                distance.from.land.allowed.to.stay = c(-Inf, Inf))

#Create a proposal
#Create an array of settings and data that incorporates all the objects created during earlier steps.
#1. the light data with the detected twilight events (Proc.data)
#2. the spatial parameters (grid)
#3. geo corrdinates of the starting location (start)
#4. the calibration parameters (Calibration)

Sys.time()
a<-Sys.time()

all.in<-make.prerun.object(Proc.data = Light.Data, 
                           Grid = Grid, 
                           start = CapLocs, # c(Longitude, Latitude)
                           Calibration = Calibration)

Sys.time()-a

#Running the particle filter
#The following parameters can be preset:
#1. number of particles (le4 is recommended for test and le6 for the analysis)
#2. known.last - TRUE if you know the track ends where it begins (FALSE is default)
#3. check.outliers - TRUE for the "on the fly" discard of outliers (only recommended to make pretty maps)

nParticles = 1000    #just a quick trial

a <- Sys.time()     

Result <- run.particle.filter(all.out = all.in, 
                              threads = -1,
                              nParticles = nParticles, 
                              known.last = TRUE,
                              precision.sd = 25, 
                              check.outliers = FALSE)
Sys.time() - a

#Save results
save(Result, file="data/FLightR_results.RData")

#Results

names(Result)

names(Result$Spatial)

names(Result$Results)

#various tables
print(Result$Results$Quantiles[1:5,],digits=3)

str(Result$Results$Quantiles)

head(Result$Results$Movement.results)

#Plot a map

map.FLightR.ggmap(Result)

#Plot and save a map

map.FLightR.ggmap(Result, save.options = list(filename = "data/FLightR.map.pdf"))

#Plot lat and long throughout tracking period

plot.lon.lat(Result)





