# Converts a light wavelength (in nm) into an RGB color name
# Author	 : RL <http://codingmess.blogspot.fr/2009/05/conversion-of-wavelength-in-nanometers.html>
# Adapted by : Sylvain Mareschal <maressyl@gmail.com>
wav2RGB <- function(wav) {
	# Initialize
	SSS <- R <- G <- B <- integer(length(wav))
	
	# Color
	zone <- wav >= 380 & wav < 440
	R[zone] <- - (wav[zone] - 440) / (440 - 350)
	B[zone] <- 1
	
	zone <- wav >= 440 & wav < 490
	G[zone] <- (wav[zone] - 440) / (490 - 440)
	B[zone] <- 1
	
	zone <- wav >= 490 & wav < 510
	G[zone] <- 1
	B[zone] <- - (wav[zone] - 510) / (510 - 490)
	
	zone <- wav >= 510 & wav < 580
	R[zone] <- (wav[zone] - 510) / (580 - 510)
	G[zone] <- 1
	
	zone <- wav >= 580 & wav < 645
	R[zone] <- 1
	G[zone] <- - (wav[zone] - 645) / (645 - 580)
	
	zone <- wav >= 645 & wav <= 780
	R[zone] <- 1
	
	# Intensity correction
	zone <- wav >= 380 & wav < 420
	SSS[zone] <- 0.3 + 0.7*(wav[zone] - 350) / (420 - 350)
	
	zone <- wav >= 420 & wav <= 700
	SSS[zone] <- 1
	
	zone <- wav > 700 & wav <= 780
	SSS[zone] <- 0.3 + 0.7*(780 - wav[zone]) / (780 - 700)
	
	# Color name
	out <- rgb(SSS*R, SSS*G, SSS*B)
	
	return(out)
}
