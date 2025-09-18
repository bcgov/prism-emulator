# Code for creating a custom legend that meaningfully scales differences/biases in log-transformed precipitation 
# also accomodates absolute differences in temperature
# the key line is line 22 (labels <- ...)
# Colin Mahony September 10th, 2024

elements <- c("pr", "tasmax", "tasmin")
element.names <- c("\nprecipitation (%)", "mean daily\nmax. temperature (K)", "mean daily\nmin. temperature (K)")

element="pr" # hard-coded for this example

lim.upper <- if(element=="pr") 1 else 7
lim.lower <- if(element=="pr") -1 else 0

inc=(lim.upper-lim.lower)/100
breaks=seq(lim.lower, lim.upper+inc, inc)
ColScheme=if(element=="pr") colorRampPalette(rev(hcl.colors(5,"Blue-Red 3")))(length(breaks)-1) else colorRampPalette(hcl.colors(5,"Blue-Red 3"))(length(breaks)-1)
pct <- if(element=="pr") 100 else 1

plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))  
xl <- 0.1; yb <- 0.4; xr <- 0.9; yt <- 0.6
rect(head(seq(xl,xr,(xr-xl)/length(ColScheme)),-1), yb,  tail(seq(xl,xr,(xr-xl)/length(ColScheme)),-1),  yt,  border=NA, col=ColScheme)
rect(xl,  yb,  xr,  yt)
labels <- if(element=="pr") paste(round(2^seq(lim.lower,lim.upper,(lim.upper-lim.lower)/4), 2)*pct-pct, "%", sep="") else round(seq(lim.lower,lim.upper,(lim.upper-lim.lower)/2), 2)*pct
text(seq(xl,xr,(xr-xl)/(length(labels)-1)),rep(yb,length(labels)),labels,pos=1,cex=1.5,font=1, offset=0.5)
text(mean(c(xl,xr)), yt+0.01, paste("Difference in", element.names[which(elements==element)]), pos=3, cex=1.5, font=2)

