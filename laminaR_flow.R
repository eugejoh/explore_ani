## Laminar Flow in Pipe
# Reynolds number (Re) < 2100 in horizontal pipe
# Hagen-Poiseuille's Law
# Qv = (pi*R^4*Pdelta) / (8*mu*l)

# flow rate is proportional to pressure drop, pipe radius (^4)
#     inversely proportional to viscosity, and pipe length

# Packages ####
library(pracma) #load practical math package
library(ggplot2) #data viz
library(scales) #plot scales
library(dplyr) #data wrangling
library(viridis) #data viz palette
library(animation) #animation
library(gganimate) #ggplot animation

## Initial Conditions ####
# Pressure difference (delta P)
Pi <- 10 #initial pressure [Pa]
Pf <- 50 #final pressure [Pa]

## Boundary Conditions #
# Pipe and fluid characteristics
mu <- 8.94e-4 #viscosity for water. units: [Pa*s] at temp = 25 C
R <- 0.1 #radius of pipe [m]
l <- 4 #length of pipe [m]

N <- 100 #radial resolution
r <- seq(-R,R,length.out=N) #diameter of pipe
S <- pi*(R^2) #area of pipe axial

# Pressure-time dependence
time <- seq(0,19) #length 20
Pit <- 0.2*time*Pi 
Pft <- 0.2*time*Pf
Pdiff <- abs(Pft-Pit) #pressure difference over time


# Create Data Viz ####
# plot shows V (average speed) and r (pipe radius)
# ((Pdiff)*(R^2)/(4*mu*l))*(1-(r^2)/R^2) #eqn for V
# out1 <- list()
# for (i in seq_along(Pdiff)) {
#   out1[[i]] <- ((Pdiff[i])*(R^2)/(4*mu*l))*(1-(r^2)/R^2)
# }

out1 <- lapply(seq_along(Pdiff), function(x) {
  ((Pdiff[x])*(R^2)/(4*mu*l))*(1-(r^2)/R^2)
})

out1_df <- data.frame(matrix(unlist(out1),ncol=length(time)),stringsAsFactors = FALSE)
names(out1_df) <- seq_along(names(out1_df))

df_t <- tidyr::gather(out1_df,key="time",value="V")

df_t2 <- cbind(df_t,r) %>% 
  mutate(time = as.integer(time)) %>% #change to integer type
  mutate(Qv = V/(1-(r^2)/R^2)) #add flowrate parameter

p1 <- ggplot(data=df_t2, aes(y=r,x=V,col=Qv,frame=time,cumulative=TRUE)) +
  geom_path() + geom_point(alpha=0.05) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_colour_viridis(option = "C", discrete = FALSE) +
  theme_dark() + theme(axis.title = element_text(size=14),
                       axis.text = element_text(size=12),
                       plot.title = element_text(size=16,face = "bold")) +
  labs(title = "Hagen-Poiseuille's Law: Velocity Profile t =", 
       y = "Pipe Radius (m)", x = "Velocity (m/s)")
p1

gganimate(p1,filename = "Vprofile1.gif", title_frame = TRUE,
          interval=0.2, ani.width = 1100, ani.height = 400)
