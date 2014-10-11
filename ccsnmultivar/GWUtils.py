#!/usr/bin/python

import numpy as np
import scipy as sp
from sys import *

######################## Define some vars #####################################
N_fd = 6145 			# Num freq samples
dF = 1./3. 			# Frequency resolution of data
lowBin = int(10./dF) + 1 	# Bin number corresponding to 10Hz - low 
				# frequency cut-off
AWGS84 = 6.378137e6 		# Semimajor axis of earth [m]
BWGS84 = 6.356752314e6 		# Semiminor axis of earth [m] 

EpochJ2000_0_UTC = 630763213 	# UTC time at the J2000.0 epoch
LeapSeconds_2012_EpochJ2000 = 3 # Number of leap seconds since J2000.0 epoch -
        			# Will need to update after next leap second!
secPerDay = 86400 		# Num secs/day
c = 299792458.	 		# Speed of light [ms^-1]
###############################################################################


######################## Detector function defs ###############################
def InitDetectorNetwork(numDets, detNetwork, LIGO3FLAG=0):
    """ InitDetectorNetwork - function to initialise desired detector network
        specified as input.

	numDets - Number of IFOs in the detector network to be considered, can
        take any value in the range [1,5].
        detNetwork - 5 character string to signify which detector network should
        be considered.  1st, 2nd, 3rd, 4th and 5th characters represent H1, L1,
        V1, I1 and K1 detectors respectively.
	LIGO3FLAG - Set to 1 to use LIGO3 PSD instead of aLIGO PSD for H1, L1
        and I1 detectors.

	Returns Noise - (N_fd,numDets) array containing simulated Gaussian detector
			noise for each detector in the network.
		PSD - (N_fd,numDets) array containing the noise PSD for each
		      detector in the network.
		detNames - (numDets) list containing names of detectors in the
			   network.
		detPosVector - (3,numDets) array containing position vector of
			       detector vertex for each detector in the network
			       in geocentric coordinates.
		detRespTensor - (3,3,numDets) array containing detector response
				tensor for each detector in the network. 

        Sarah Gossan 2012. Last updated 02/18/14. """

    # Create arrays for output
    Noise = np.zeros((N_fd,numDets),complex)
    PSD = np.zeros((N_fd,numDets),float)
    detPosVector = np.zeros((3,numDets),float)
    detRespTensor = np.zeros((3,3,numDets),float)

    # Create detector names array
    detNames = []
    if detNetwork[0] == 'H':
	detNames.append('H1')
    if detNetwork[1] == 'L':
	detNames.append('L1')
    if detNetwork[2] == 'V':
	detNames.append('V1')
    if detNetwork[3] == 'I':
	detNames.append('I1')
    if detNetwork[4] == 'K':
	detNames.append('K1')

    # Get detector PSD, Noise, position vector, response tensor and time delay
    # from geocenter dependent on working domain
    for detIdx in range(0,numDets):
	# PSDs
	PSD[:,detIdx] = GetDetectorPSD(detNames[detIdx])
        Noise[:,detIdx] = GenerateGaussianNoise(PSD[:,detIdx])
       	# Position vector and response tensor
	[detRespTensor[:,:,detIdx],detPosVector[:,detIdx]] = \
	GetDetectorData(detNames[detIdx])

    return Noise,PSD,detNames,detPosVector,detRespTensor


def GetDetectorPSD(detectorName, LIGO3FLAG=0):
    """ GetDetectorPSD - function to return the PSD of the detector described by
	detectorName.

	detectorName - Name of required GW IFO.  Can be 'H1', 'L1', 'V1', 'I1',
                       or 'K1'.
	LIGO3FLAG - Set to 1 to use LIGO3 PSD instead of aLIGO PSD for H1, L1
        and I1 detectors.

	Returns PSD - the noise PSD for the detector described by detectorName.

        Sarah Gossan 2012. Last updated 02/18/14. """

    # H1, L1 or I1 - use aLIGO ZERODET_HIGHP configuration (fLow = 9Hz)
    # Noise curve is super complicated so just load in PSD file for now
    if detectorName == 'H1' or detectorName == 'L1' or detectorName == 'I1':
        if LIGO3FLAG:
            # Read in PSD
            PSD = np.loadtxt('LIGO3_PSD.txt')
        else:
            # Read in PSD
            PSD = np.loadtxt('ZERO_DET_high_P_PSD.txt')
	# Only want second column of file
        PSD = PSD[:,1]
    # V1 - use analytical expression for AdVirgo (fLow = 10Hz)
    elif detectorName == 'V1':
        # Use analytical expression from arXiv:1202.4031v2
	x = np.linspace(0,N_fd-1,num=N_fd)*dF/300.
        x[0] = x[1] # Not going to use f=0Hz component anyway, but this stops 
		    # the log fn complaining
	x = np.log(x)
	xSq = x*x
     	asd = 1.259e-24*(0.07*np.exp(-0.142 - 1.437*x + 0.407*xSq) + \
                         3.1*np.exp(-0.466 - 1.043*x - 0.548*xSq) + \
                         0.4*np.exp(-0.304 + 2.896*x - 0.293*xSq) + \
                         0.09*np.exp(1.466 + 3.722*x - 0.984*xSq))
        PSD = asd**2 
    # K1 - use analytical expression for KAGRA (fLow = 10Hz) 
    elif detectorName == 'K1':
	# Use analytical expression from arXiv:1202.4031v2
	x = np.linspace(0,N_fd-1,num=N_fd)*dF/300.
        x[0] = x[1] # Not going to use f=0Hz component anyway, but this stops 
                    # the log fn complaining
        x = np.log(x)
        xSq = x*x
	asd = 6.499e-25*(9.72e-9*np.exp(-1.43 - 9.88*x - 0.23*xSq) + \
		    	 1.17*np.exp(0.14 - 3.10*x - 0.26*xSq) + \
			 1.70*np.exp(0.14 + 1.09*x - 0.013*xSq) + \
			 1.25*np.exp(0.071 + 2.83*x - 4.91*xSq))
        PSD = asd**2
          
    return PSD


def GenerateGaussianNoise(PSD):
    """ GenerateGaussianNoise - function to generate Fourier domain Gaussian 
        detector noise colored by a detector PSD.
        
        PSD - Noise power spectral density with which to color the Gaussian
        noise.

        Returns Noise - complex Fourier domain Gaussian detector noise colored 
                        by a detector PSD.

        Sarah Gossan 2012. Last updated 02/18/14. """

    Noise = np.zeros((N_fd),complex)
    # Generate noise from PSD 
    Real = np.random.randn(N_fd)*np.sqrt(PSD/(4.*dF))
    Imag = np.random.randn(N_fd)*np.sqrt(PSD/(4.*dF))
    Noise = Real + 1j*Imag

    return Noise


def GetDetectorData(detectorName):
    """ GetDetectorData - function to return the locations and detector response
 	tensor of the detector described by detectorName.

	detectorName - Name of required GW IFO.  Can be 'H1', 'L1', 'V1', 'I1'
	or 'K1'.

	Returns detRespTensor - Detector response tensor of the IFO.
	    	detPosVector - Position vector of the IFO.

        Sarah Gossan 2012. Last updated 02/18/14. """

    # Allocate direction vectors and position of IFO given detectorName
    detPosVector = np.zeros(3)
    nx = np.zeros(3)
    ny = np.zeros(3)

    ################# aLIGO HANFORD #############
    if detectorName == 'H1':
	nx[0] = -0.22389266154
	nx[1] = 0.79983062746
	nx[2] = 0.55690487831

	ny[0] = -0.91397818574
	ny[1] = 0.02609403989
	ny[2] = -0.40492342125

	detPosVector[0] = -2.16141492636e6
	detPosVector[1] = -3.83469517889e6
	detPosVector[2] = 4.60035022664e6
    ################# aLIGO LIVINGSTON ##########
    elif detectorName == 'L1':
	nx[0] = -0.95457412153
	nx[1] = -0.14158077340
	nx[2] = -0.26218911324

	ny[0] = 0.29774156894
	ny[1] = -0.48791033647 
	ny[2] = -0.82054461286

	detPosVector[0] = -7.42760447238e4
	detPosVector[1] = -5.49628371971e6
	detPosVector[2] = 3.22425701744e6
    ################# AdvVIRGO ##################
    elif detectorName == 'V1':
	nx[0] = -0.70045821479
        nx[1] = 0.20848948619
        nx[2] = 0.68256166277

        ny[0] = -0.05379255368 
        ny[1] = -0.96908180549
        ny[2] = 0.24080451708

        detPosVector[0] = 4.54637409900e6
        detPosVector[1] = 8.42989697626e5
        detPosVector[2] = 4.37857696241e6
    ################# INDIGO ####################
    elif detectorName == 'I1':
	# Not cached LIGO detector - construct x, nx and ny from other
	# quantities
	vxLat = (14. + 14./60.)*(np.pi/180.) # Vertex latitude
	vxLon = (76. + 26./60.)*(np.pi/180.) # Vertex longitude
	vxElev = 0. # Vertex elevation
	xAlt = 0. # Altitude of x-arm
	yAlt = 0. # Altitude of y-arm
	xAz = np.pi/2. # Azimuth of x-arm
	yAz = 0. # Azimuth of y-arm

	cosLat = np.cos(vxLat)
        sinLat = np.sin(vxLat)
        cosLon = np.cos(vxLon)
        sinLon = np.sin(vxLon)

        ellDenom = np.sqrt(AWGS84*AWGS84*cosLat*cosLat + \
                           BWGS84*BWGS84*sinLat*sinLat)
        locRho = cosLat*(AWGS84*AWGS84/ellDenom + vxElev)
        
        # Position
        detPosVector[0] = locRho*cosLon
        detPosVector[1] = locRho*sinLon
        detPosVector[2] = sinLat*(BWGS84*BWGS84/ellDenom + vxElev)

        # x-Arm vector
        cosxAlt = np.cos(xAlt)
        sinxAlt = np.sin(xAlt)
        cosxAz = np.cos(xAz)
        sinxAz = np.sin(xAz)
        uxNorth = cosxAlt*cosxAz
        uxEast = cosxAlt*sinxAz
        uxRho = -sinLat*uxNorth + cosLat*sinxAlt

        nx[0] = cosLon*uxRho - sinLon*uxEast
        nx[1] = sinLon*uxRho + cosLon*uxEast
        nx[2] = cosLat*uxNorth + sinLat*sinxAlt

        # y-Arm vector
        cosyAlt = np.cos(yAlt)
        sinyAlt = np.sin(yAlt)
        cosyAz = np.cos(yAz)
        sinyAz = np.sin(yAz)
        uyNorth = cosyAlt*cosyAz
        uyEast = cosyAlt*sinyAz
        uyRho = -sinLat*uyNorth + cosLat*sinyAlt

        ny[0] = cosLon*uyRho - sinLon*uyEast
        ny[1] = sinLon*uyRho + cosLon*uyEast
        ny[2] = cosLat*uyNorth + sinLat*sinyAlt

    ################# KAGRA #####################
    elif detectorName == 'K1':
	# Not cached LIGO detector - construct x, nx and ny from other
        # quantities
        vxLat = (36.25)*(np.pi/180.) # Vertex latitude
        vxLon = (137.18)*(np.pi/180.) # Vertex longitude
        vxElev = 0. # Vertex elevation
        xAlt = 0. # Altitude of x-arm
        yAlt = 0. # Altitude of y-arm
        xMid = 1500. # Midpoint of x-arm
        yMid = 1500. # Midpoint of y-arm
        xAz = 19.*(np.pi/180.) + np.pi/2. # Azimuth of x-arm
        yAz = 19.*(np.pi/180.) # Azimuth of y-arm

        cosLat = np.cos(vxLat)
	sinLat = np.sin(vxLat)
	cosLon = np.cos(vxLon)
	sinLon = np.sin(vxLon)

	ellDenom = np.sqrt(AWGS84*AWGS84*cosLat*cosLat + \
		           BWGS84*BWGS84*sinLat*sinLat)
	locRho = cosLat*(AWGS84*AWGS84/ellDenom + vxElev)

	# Position
	detPosVector[0] = locRho*cosLon
	detPosVector[1] = locRho*sinLon
	detPosVector[2] = sinLat*(BWGS84*BWGS84/ellDenom + vxElev)

	# x-Arm vector
	cosxAlt = np.cos(xAlt)
	sinxAlt = np.sin(xAlt)
	cosxAz = np.cos(xAz)
	sinxAz = np.sin(xAz)       
	uxNorth = cosxAlt*cosxAz
	uxEast = cosxAlt*sinxAz
	uxRho = -sinLat*uxNorth + cosLat*sinxAlt

	nx[0] = cosLon*uxRho - sinLon*uxEast
	nx[1] = sinLon*uxRho + cosLon*uxEast
	nx[2] = cosLat*uxNorth + sinLat*sinxAlt

	# y-Arm vector
	cosyAlt = np.cos(yAlt)
        sinyAlt = np.sin(yAlt)
        cosyAz = np.cos(yAz)
        sinyAz = np.sin(yAz)
	uyNorth = cosyAlt*cosyAz
        uyEast = cosyAlt*sinyAz
        uyRho = -sinLat*uyNorth + cosLat*sinyAlt

        ny[0] = cosLon*uyRho - sinLon*uyEast
        ny[1] = sinLon*uyRho + cosLon*uyEast
        ny[2] = cosLat*uyNorth + sinLat*sinyAlt

    # Compute detRespTensor as 
    # detRespTensor = 0.5(nx \outerprod nx - ny \outerprod ny)
    detRespTensor = np.zeros((3,3),float)
    for i in range(0,3):
	for j in range(0,3):
	    detRespTensor[i,j] = (nx[i]*nx[j] - ny[i]*ny[j])*0.5

    return detRespTensor,detPosVector


def ComputeAntennaResponse(theta, phi, psi, detRespTensor):
    """ ComputeAntennaResponse - function to compute F+,x of a detector
	characterised by detRespTensor to a source located at (theta,phi,psi) 
	in geocentric coordinates.

	theta - Polar angle of source's sky position in geocentric
                coordinates.
	phi -  Azimuthal angle of source's sky position in geocentric
               coordinates.
 	psi - Polarisation angle of gravitational waves, measured
	      counter-clockwise about direction of propagation from line of 
	      constant theta pointing to decreasing phi to positive x-axis of 
	      source coordinates.
        detRespTensor - (3,3) array containing detector response tensor.

	Returns Fp - Antenna response of detectorName to the plus GW
		polarisation of a source at (RA,Dec,psi).
		Fc - Antenna response of detectorName to the cross GW
		polarisation of a source at (RA,Dec,psi).

        Sarah Gossan 2012. Last updated 02/18/14. """

    # Define axes of wave frame in terms of axes of Earth-centered axes and
    # source angles
    X = np.zeros(3)
    Y = np.zeros(3)
    X[0] = np.sin(phi)*np.cos(psi) - np.sin(psi)*np.cos(phi)*np.cos(theta)
    X[1] = -np.cos(phi)*np.cos(psi) - np.sin(psi)*np.sin(phi)*np.cos(theta)
    X[2] = np.sin(psi)*np.sin(theta)
    Y[0] = -np.sin(phi)*np.sin(psi) - np.cos(psi)*np.cos(phi)*np.cos(theta)
    Y[1] = np.cos(phi)*np.sin(psi) - np.cos(psi)*np.sin(phi)*np.cos(theta)
    Y[2] = np.sin(theta)*np.cos(psi)
    
    # Construct polarisation vectors e+,x as
    # e+ = (X \outerprod X - Y \outerprod Y)
    # ex = (X \outerprod Y + Y \outerprod X)
    ep = np.zeros((3,3),float)
    ec = np.zeros((3,3),float)
    
    for i in range(0,3):
        for j in range(0,3):
	    ep[i,j] = X[i]*X[j] - Y[i]*Y[j]
	    ec[i,j] = X[i]*Y[j] + Y[i]*X[j]

    # Construct Fp,Fc as 
    # Fp = \sum_ij(D_ij ep_ij) 
    # Fc = \sum_ij(D_ij ec_ij)
    Fp = 0.
    Fc = 0.
    for i in range(0,3):
	for j in range(0,3):
	    Fp += detRespTensor[i,j]*ep[i,j]
 	    Fc += detRespTensor[i,j]*ec[i,j]

    return Fp,Fc


def ComputeArrivalTimeDifference(RA, Dec, GPSTime, Det1Pos, Det2Pos=([0,0,0])):
    """ ComputeArrivalTimeDifference - function to compute the difference in
	arrival times in a GW signal from a source located at celestial coordinates
	(RA,Dec) that reaches the geocenter at GPSTime, between two detectors 
	located at Det1Pos and Det2Pos.

	RA - Right ascension of source [radians].
	Dec - Declination of source [radians].
        GPSTime - GPS time that the GW signal from the source reached the
		  geocenter.
	Det1Pos - (3,) array. Position vector of detector 1.

	################ Optional inputs #####################
	Det2Pos - (3,) array. Position vector of detector 2 (default =
	          ([0,0,0]), the position vector of the geocenter.  In this
		  case, ComputeArrivalTimeDifference returns the time delay of 
		  the source arrival time from the geocenter).

	Returns deltaT - the difference in arrival times between detectors
			 located at Det1Pos and Det2Pos.  This is positive when
			 the signal arrives at Det1Pos *after* Det2Pos, and negative 
			 when the signal arrives at Det1Pos before Det2Pos.

        Sarah Gossan 2012.  Adapted from TimeDelay.c, written by Jolien Creighton, 
	David Chin, Steven Fairhurst, Kipp Cannon, Alexander Dietz, Drew Keppel
	2007. Last updated 02/18/14. """

    # Compute GMST for GPSTime
    GMST = ComputeGMST(GPSTime)

    # Get Greenwich hour angle from RA and GMST
    GHA = GMST - RA

    # Compute unit vector pointing from geocenter to source
    eHatSrc = np.zeros(3)
    eHatSrc[0] = np.cos(Dec)*np.cos(GHA)
    eHatSrc[1] = -np.cos(Dec)*np.sin(GHA)
    eHatSrc[2] = np.sin(Dec)

    # Get position of detector 2 with respect to detector 1
    Det21Pos = Det2Pos - Det1Pos

    # Dot source vector into relative detector position and divide by speed of
    # light to get difference in arrival time
    deltaT = np.dot(eHatSrc,Det21Pos)/c

    return deltaT


def ComputeLightTravelTime(Det1Pos, Det2Pos):
    """ ComputeLightTravelTime - function to compute light travel time between
	two GW detectors at positions Det1Pos and Det2Pos.

	Det1Pos - (3,) array. Position vector of detector 1.
	Det2Pos - (3,) array. Position vector of detector 2.

	Returns travelTime - Light travel time between two GW detectors at
			     positions Det1Pos and Det2Pos [s].

        Sarah Gossan 2012. Adapted from TimeDelay.c, written by Jolien
	Creighton, David Chin, Steven Fairhurst, Kipp Cannon, Alexander Dietz, 
	Drew Keppel 2007."""

    # Get relative position vector
    Det21Pos = Det2Pos - Det1Pos
  
    # Dot difference vector into itself to get magnitude of detector separation
    dist = np.sqrt(np.dot(Det21Pos,Det21Pos))

    # Normalise with speed of light
    travelTime = dist/c

    return travelTime


def ConvertRADecToThetaPhi(RA, Dec, GPSTime): 
    """ ConvertRADecToThetaPhi - function to convert right ascension and 
        declination of an astronomical source to theta and phi in the geocentric 
        coordinate system. 
            
        RA - Right ascension of source. 
        Dec - Declination of source. 
        GPSTime - GPS time that the GW signal from the source reached the
		  geocenter. 
                
        Returns Theta - Polar angle of source's sky position in geocentric 
                        coordinates. 
                Phi - Azimuthal angle of source's sky position in geocentric 
                      coordinates. 
            
        Sarah Gossan 2012. """ 
                
    # Get Greenwich mean sidereal time corresponding to GPSTime 
    GMST = ComputeGMST(GPSTime) 
            
    # RA -> Phi
    Phi = RA - GMST
            
    # Dec -> Theta
    Theta = np.pi/2. - Dec
            
    return Theta,Phi


def ComputeGMST(GPSTime):
    """ ComputeGMST - function to compute the Greenwich mean sidereal time from
                      the GPS time.

        GPSTime - GPS time that the GW signal from the source reached the
                  geocenter.

        Returns GMST - the Greenwich mean sidereal time corresponding to
                       GPSTime.

        Sarah Gossan 2012. """

    # Difference in Julian Date between GPSTime and the J2000.0 epoch
    # Subtract half a day as Julian Days start at noon
    D = np.round((GPSTime - EpochJ2000_0_UTC)/secPerDay) - 0.5

    # Number of Julian centuries since J2000.0 epoch
    T = D/36525

    # 1st approximation to GMST without corrections for leap seconds
    GMST0 = 6.697374558 + 2400.051336*T + 0.000025862*T*T
    # Keep within range [0,24]hrs
    GMST0 = np.mod(GMST0,24)

    # Corrections for leap seconds
    UTCSec = GPSTime - EpochJ2000_0_UTC - secPerDay/2 - \
             LeapSeconds_2012_EpochJ2000
    UTCHr = np.mod(UTCSec/3600,24)

    # Add corrections and keep within [0,24]hr range
    GMST = GMST0 + UTCHr*1.002737909
    GMST = np.mod(GMST,24)

    # Convert from hours to degrees to radians
    GMST *= 15.*(np.pi/180.)

    return GMST
###############################################################################

############################# Stats function defs #############################
def LikelihoodFunction(Template, Data, PSD, detRespP, detGCDelay=0):
    """ LikelihoodFunction - function to calculate the likelihood of livePoint,
        given Data.

        Template - (N_fd) complex array containing Fourier domain trial signal.
        Data - (N_fd) complex array containing Fourier domain GW data.
        PSD - Noise power spectral density for a gravitational wave detector.
        detRespP - Antenna response to the plus GW polarisation for the
                   detector.
	detGCDelay - Time delay of detector from geocenter (default = 0, use
		     detGCDelay only if computing logL for more than one 
		     detector.

        Returns logL of Template.

        Sarah Gossan 2012. Last updated 02/18/14. """

    # Correct template for geocenter delay and antenna response function
    if detGCDelay:
        phaseGCDelay = -2.*np.pi*np.linspace(0,N_fd-1,num=N_fd)*dF*detGCDelay*1j 
	Template *= phaseGCDelay
    Template *= detRespP
    # Calculate logL - simple Gaussian
    logL = -2.*dF*np.sum(pow(abs(Data[lowBin:] - Template[lowBin:]),2.)/\
	   PSD[lowBin:])

    return logL
###############################################################################
