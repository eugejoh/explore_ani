# Explore animation package

library(pracma) #load practical math package
library(ggplot2)

library(animation)
library(gganimate)
ani.options(interval = 0.15)

last <- function(x) { return( x[length(x)] ) } #returns last element in vector

# Fick's First Law: Diffusion ####
# MATLAB code equivalent
xline <- seq(-20e-6,20e-6,length.out = 40) #x-axis vector
yline <- seq(-20e-6,20e-6,length.out = 40)
D <- 1e-14 #material diffusion coeffecient
t_end <- 200 #end of simulation time

# x1 <- xline[1]
# x_end <- last(xline)
# for (t in 1:t_end) {
#   c <- 50*(1-erf(xline/(2*sqrt(D*t)))) #use of erf() fnc, error of Phi function, variant of cumulative normal
#   plot(xline,c,type="b")
# }

time <- 0.1
diffusion_fnc <- function(D,t=time) { #use of erf() fnc, error of Phi function 
  return(
    (1-erf(xline/(2*sqrt(1e-4*D*t))))
  ) #variant of cumulative distribution function of normal distribution
}

df <- data.frame(x=xline,
                 cx=diffusion_fnc(D1,time),
                 cy=diffusion_fnc(D3,time))
# dataframe method in R
df <- data.frame(x=xline,
                 y=yline,
                 cx=50*(1-erf(xline/(2*sqrt(1e-13*time)))),
                 cy=50*(1-erf(yline/(2*sqrt(D*time))))
)

p1 <- ggplot(data=df,aes(x=x,y=cx,col=cx,frame=x,cumulative = TRUE)) +
  geom_point(aes(cumulative = TRUE)) +
  geom_line(aes(cumulative = TRUE)) +
  scale_x_continuous(breaks=pretty_breaks(n=6),labels=scales::comma)
p1
gganimate(p1, title_frame = FALSE, ani.height = 600, ani.width = 800)


## Laminar Flow in Pipe
# Reynolds number (Re) < 2100 in horizontal pipe
# Hagen-Poiseuille's Law
# Qv = (pi*R^4*Pdelta) / (8*mu*l)
# flow rate is proportional to pressure drop, pipe radius (^4)
#     inversely proportional to viscosity, and pipe length

## Initial Conditions ####
# Pressure difference (delta P)
Pi <- 10 #initial pressure [Pa]
Pf <- 50 #final pressure [Pa]

## Boundary Conditions ####
# Pipe and fluid characteristics
mu <- 8.94e-4 # visocity for water. units: [Pa*s] at temp = 25 C
R <- 0.1 #radius of pipe [m]
l <- 1 #length of pipe [m]

# initial conditions
Pi <- 10 # initial pressure
Pf <- 50 # final pressure
mu <- 8.94e-4 #mu dynamic visocity
R <- 0.10 #radius of pipe
l <- 1.00 #length of pipe

# time dependent
t <- seq(0,19)
Pit <- 0.2*t*Pi
Pft <- 0.2*t*Pf
Pdiff <- abs(Pft-Pit)

# boundary conditions
N <- 100 #radial resolution
r <- seq(-R,R,length.out=N) #diameter of pipe
S <- pi*(R^2) #area of pipe axial
Q <- Vmean*S #volumetric flow rate

# Create Data Viz
# plot shows V (average speed) and r (pipe radius)
((Pdiff)*(R^2)/(4*mu*l))*(1-(r^2)/R^2) #eqn for V

out1 <- list()
for (i in seq_along(Pdiff)) {
  out1[[i]] <- ((Pdiff[i])*(R^2)/(4*mu*l))*(1-(r^2)/R^2)
}

out1_df <- data.frame(matrix(unlist(out1),ncol=25),stringsAsFactors = FALSE)
names(out1_df) <- seq_along(names(out1_df))

df_t <- tidyr::gather(out1_df,key="t",value="V")

df_t2 <- cbind(df_t,r) %>% mutate(t = as.integer(t))

ggplot(data=df_t2,aes(x=r,y=V,col=t)) +
  geom_path() + geom_point(alpha=0.1) +
  scale_color_viridis(option = "C",discrete = FALSE)

gganimate(
  ggplot(data=df_t2, aes(y=r,x=V,col=t,frame=t,cumulative=TRUE)) +
    geom_path() + geom_point(alpha=0.1) +
    scale_colour_viridis(option = "C", discrete = FALSE),
  title_frame = FALSE, interval=0.2,
  ani.width = 600, ani.height = 250
)


out2 <- ggplot(test1_df,aes(y=r,x=V,frame=.frame)) +
  geom_path() + geom_point(alpha=0.1) +
  scale_colour_viridis(option = "C", discrete = FALSE)

gganimate(out2, title_frame = FALSE, interval=0.2)

x1 <- -10:10

df1 <- data.frame(x=x1,
                  y=x1^2,
                  t=1)

df2 <- data.frame(x=x1,
                  y=0.5*x1^2,
                  t=2)

df3 <- data.frame(x=x1,
                  y=0.25*x1^2,
                  t=3)

df4 <- data.frame(x=x1,
                  y=0.1275*x1^2,
                  t=4)

df5 <- data.frame(x=x1,
                  y=0.06375*x1^2,
                  t=5)

df_t <- rbind(df1,df2,df3,df4,df5)

out1 <- ggplot(data=df_t,aes(x=x,y=y,frame=t)) +
  geom_point(aes(col=t))

gganimate(out1,title_frame = FALSE, interval=0.5)


# DEPRECATED ####
# boundary conditions
N <- 100 #radial resolution
r <- seq(-R,R,length.out=N) #diameter of pipe
S <- pi*(R^2) #area of pipe axial
Vmax  <- abs(Pi-Pf)*(R^2)/(4*mu*l) # pressure diff * r^2 / 4*mu*L
Vmean <- Vmax/2
V <- Vmax*(1-(r^2)/R^2) #
Q <- Vmean*S #volumetric flow rate
plot(V,r,type="b")

df_poiseuille <- data.frame(Pit,Pft,Pdiff,r,
                            Vmax = Pdiff*(R^2)/(4*mu*l),
                            V = Vmax*(1-(r^2)/(R^2)))


library(pracma) #load practical math package
library(ggplot2) #data viz
library(scales)
library(dplyr) #data wrangling
library(purrr)

library(animation) #animation
library(gganimate) #ggplot animation
ani.options(interval = 0.15)

# nifty function for last element in vector
last <- function(x) { return( x[length(x)] ) }

# Fick's First Law: Diffusion ####
# MATLAB code equivalent
xline <- seq(-20e-6,20e-6,length.out = 40) #x-axis vector
#material diffusion coeffecients
D1 <- 1.02e-5 #benzene dissolved and water
D2 <- 1.25e-5 #chlorine dissolved and water
D3 <- 0.97e-5 #propane dissolved and water
D4 <- 1.64e-5 #ammonia dissolved and water
D_all <- c("benzene"=D1,"chlorine"=D2,"propane"=D3,"ammonia"=D4)
time <- 20 #time length

diffusion_fnc <- function(D,t=time) { #use of erf() fnc, error of Phi function 
  return(
    50*(1-erf(xline/(2*sqrt(D*t))))
    ) #variant of cumulative distribution function of normal distribution
}

50*(1-erf(xline/(2*sqrt(D*t))))

df <- map_dfc(D_all,~data.frame(diffusion_fnc(xline,.x,time)))
names(df) <- names(D_all)
df

# create a data frame method in R
df <- data.frame(x=xline, 
                 c1=diffusion_fnc(1e-14,t),
                 c2=diffusion_fnc(D2,t)
) 

p1 <- ggplot(data=df) +
  geom_rect(aes(ymax=Inf, ymin=-Inf,
                         xmax=0, xmin=-2e-05), col="grey45", fill = "#5A5A5A", alpha=0.01) +
  geom_rect(aes(ymax=Inf, ymin=-Inf,
                         xmax=2e-05,xmin=0), col="grey45", fill = "#2262A3", alpha=0.01) +
  geom_point(aes(x=x, y=c1, col=c1, frame=x, cumulative = TRUE), size = 2) +
  geom_line(aes(x=x, y=c1, col=c1, frame=x, cumulative = TRUE)) +
  geom_vline(xintercept = 0,col="grey15",linetype=2) +
  scale_x_continuous(breaks=pretty_breaks(n=6),labels=scales::comma) +
  scale_colour_distiller(type = "seq",palette = 1) +
  labs(title = "Fick's 1st Law of Diffusion",
       y = "Concentration (mol)", x = "Distance (mm)")
p1
gganimate(p1, title_frame = FALSE, ani.height = 300, ani.width = 500)


## Laminar Flow in Pipe
# Reynolds number (Re) < 2100 in horizontal pipe
# Hagen-Poiseuille's Law
# Qv = (pi*R^4*Pdelta) / (8*mu*l)
# flow rate is proportional to pressure drop, pipe radius (^4)
#     inversely proportional to viscosity, and pipe length

## Initial Conditions ####
# Pressure difference (delta P)
Pi <- 10 #initial pressure [Pa]
Pf <- 50 #final pressure [Pa]

## Boundary Conditions ####
# Pipe and fluid characteristics
mu <- 8.94e-4 # visocity for water. units: [Pa*s] at temp = 25 C
R <- 0.1 #radius of pipe [m]
l <- 1 #length of pipe [m]

# initial conditions
Pi <- 10 # initial pressure
Pf <- 50 # final pressure
mu <- 8.94e-4 #mu dynamic visocity
R <- 0.10 #radius of pipe
l <- 1.00 #length of pipe

# time dependent
t <- seq(0,19)
Pit <- 0.2*t*Pi
Pft <- 0.2*t*Pf
Pdiff <- abs(Pft-Pit)

# boundary conditions
N <- 100 #radial resolution
r <- seq(-R,R,length.out=N) #diameter of pipe
S <- pi*(R^2) #area of pipe axial
Q <- Vmean*S #volumetric flow rate

# Create Data Viz
# plot shows V (average speed) and r (pipe radius)
((Pdiff)*(R^2)/(4*mu*l))*(1-(r^2)/R^2) #eqn for V

out1 <- list()
for (i in seq_along(Pdiff)) {
  out1[[i]] <- ((Pdiff[i])*(R^2)/(4*mu*l))*(1-(r^2)/R^2)
}

out1_df <- data.frame(matrix(unlist(out1),ncol=25),stringsAsFactors = FALSE)
names(out1_df) <- seq_along(names(out1_df))

df_t <- tidyr::gather(out1_df,key="t",value="V")

df_t2 <- cbind(df_t,r) %>% mutate(t = as.integer(t))

ggplot(data=df_t2,aes(x=r,y=V,col=t)) +
  geom_path() + geom_point(alpha=0.1) +
  scale_color_viridis(option = "C",discrete = FALSE)

gganimate(
  ggplot(data=df_t2, aes(y=r,x=V,col=t,frame=t,cumulative=TRUE)) +
    geom_path() + geom_point(alpha=0.1) +
    scale_colour_viridis(option = "C", discrete = FALSE),
  title_frame = FALSE, interval=0.2,
  ani.width = 600, ani.height = 250
)