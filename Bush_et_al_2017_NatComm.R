install.packages('deSolve')
install.packages('ggplot2')
install.packages('reshape')
install.packages('gridExtra')
install.packages('scales')
library('deSolve')
library('ggplot2')
library('reshape')
library('gridExtra')
library('scales')

# plot results function
plot_results <- function(init_df = init_df){
  
  # Populaton Density Dynamics
  rec_growth <- cbind.data.frame(PB=init_df$N_PB, 
                                 SB=init_df$N_SB, 
                                 CB=init_df$N_CB)
  ss_growth <- melt(rec_growth)
  ss_growth$time <- times
  
  p1 <- ggplot(data=ss_growth, mapping=aes(x=time/24, y=value, col=variable)) + 
    geom_line(lwd=1.5)  + 
    scale_color_manual(values=c("#ae017e", "black", "#238b45")) + 
    ggtitle('Population Density') +
    theme_bw() + 
    ylab("") +
    xlab('Time (days)') +
    theme(legend.position = "bottom", legend.key.size = unit(0.8, "cm"), legend.direction = "horizontal")+ 
    ylab('Population Density') +
    theme(legend.title=element_blank()) +
    scale_y_log10(limits = c(1,NA),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    annotation_logticks(sides = 'lr')
  
  # Substrate dynamics
  rec_substrate <- cbind.data.frame(SO=init_df$SO, SR=init_df$SR, O=init_df$O, P= init_df$P)
  ss_substrate <- melt(rec_substrate)
  ss_substrate$time <- times
  
  p2 <- ggplot(data=ss_substrate, mapping=aes(x=time/24, y=value, col=variable)) + 
    geom_line(lwd=1.5) +
    theme_bw() + 
    theme(legend.position = "bottom", legend.key.size = unit(0.8, "cm"))  + 
    scale_color_manual(values=c("red", "#bf812d", "blue", "#9ecae1")) +
    ggtitle('Substrates') + 
    ylab("Substrate (µM)") +
    xlab('Time (days)') +
    theme(legend.title=element_blank())
  
  
  return(grid.arrange(p1, p2, nrow=2))
}


## multi-species dynamic growth model
community_model <- function(time, y, parms){
  with(as.list(c(y, parms)), {
    g_CB_P     <- gmax_CB * ( P / (K_CB_P + P) ) # (8.1)
    g_PB_P_SR  <- gmax_PB * ( P / (K_PB_P + P) ) *  ( SR / (K_PB_SR + SR) ) # (8.2)
    g_SB_P_SO  <- gmax_SB * ( P / (K_SB_P + P) ) *  ( SO / (K_SB_SO + SO) ) # (8.3)
    
    h_CB_SR <- 1 / ( 1 + (SR / KH_CB_SR) ) # (9.1)
    h_PB_O  <- 1 / ( 1 + (O  / KH_PB_O) )  # (9.2)
    h_SB_O  <- 1 / ( 1 + (O  / KH_SB_O) )  # (9.3)

    dN_CB  <- (g_CB_P    * h_CB_SR * N_CB) - (m_CB * N_CB)  # (1)
    dN_PB  <- (g_PB_P_SR * h_PB_O  * N_PB) - (m_PB * N_PB)  # (2)
    dN_SB  <- (g_SB_P_SO * h_SB_O  * N_SB) - (m_SB * N_SB)  # (3)
    
    dSO <-   (1/y_S_PB * g_PB_P_SR * h_PB_O  * N_PB) - (1/y_S_SB * g_SB_P_SO * h_SB_O * N_SB) + (c * O * SR) + ( a_s * (S_O_b - SO) )  # (4)
    dSR <- - (1/y_S_PB * g_PB_P_SR * h_PB_O  * N_PB) + (1/y_S_SB * g_SB_P_SO * h_SB_O * N_SB) - (c * O * SR) + ( a_s * (S_R_b - SR) )  # (5)
    dO  <-   (p_CB     * g_CB_P    * h_CB_SR * N_CB) - (c * O *SR) + ( a_O * (O_b - O) )  # (6)
    dP  <- - (1/y_P_CB * g_CB_P    * h_CB_SR * N_CB) - (1/y_P_PB * g_PB_P_SR * h_PB_O * N_PB) - (1/y_P_SB * g_SB_P_SO * h_SB_O * N_SB) + ( a_P * (P_b - P) ) # (7)
    
    
    return(list(c(dN_CB, 
                  dN_PB, 
                  dN_SB,
                  
                  dSO, 
                  dSR, 
                  dO, 
                  dP)))
  })}



# growth parameters 
gmax_CB <- 0.05 # Maximum specific growth rate of CB [hr-1]
gmax_PB <- 0.07 # Maximum specific growth rate of PB [hr-1]
gmax_SB <- 0.1 # Maximum specific growth rate of SB [hr-1]

# Half-saturation constant
K_PB_SR <- 10.0  # Half-saturation constant of PB on reduced sulfur [μM]
K_SB_SO <- 5.00  # Half-saturation constant of SB on oxidized sulfur [μM]
K_CB_P  <- 0.20  # Half-saturation constant of CB on phosphorus [μM]
K_PB_P  <- 0.50  # Half-saturation constant of PB on phosphorus [μM]
K_SB_P  <- 0.50  # Half-saturation constant of SB on phosphorus [μM]

# Half-inhibition constant
KH_CB_SR <- 300  # Half-inhibition constant of CB on reduced sulfur [μM]
KH_PB_O  <- 100  # Half-inhibition constant of PB on oxygen [μM]
KH_SB_O  <- 100  # Half-inhibition constant of SB on oxygen [μM]

# Yields
y_S_SB <- 3.33E7 # Yield of SB on oxidized sulfur [cells/µM]
y_S_PB <- 1.25E7 # Yield of PB on reduced sulfur [cells/µM]
y_P_CB <- 1.67E8 # Yield of CB on phosphorus [cells/µM]
y_P_PB <- 1.67E8 # Yield of PB on phosphorus [cells/µM]
y_P_SB <- 1.67E8 # Yield of SB on phosphorus [cells/µM]

# Production
p_CB <- 6E-9 # Production of oxygen per cyanobacterial cell [µM/cell]

# Mortality rate
m_CB <- 0.020 # Mortality rate of CB [1/h]
m_PB <- 0.028 # Mortality rate of PB [1/h]
m_SB <- 0.040 # Mortality rate of SB [1/h]

# Diffusivity 
a_s <- 0.001 # Diffusivity of sulfur [1/h]
a_O <- 8E-4 # Diffusivity of oxygen (10E-6 - 10E-2) [1/h]
a_P <- 0.010 # Diffusivity of phosphorus* [1/h]

# Background concentration
S_R_b <- 300 # Background concentration of reduced sulfur [µM]
S_O_b <- 300 # Background concentration of oxidized sulfur [µM]
O_b   <- 300 # Background concentration of oxygen [µM]
P_b   <- 9.5 # Background concentration of phosphorus (2 − 10) [µM]

c <- 4E-5 # Oxidation rate of reduced sulfur [1/µM * 1/h]


parms <- c(gmax_CB = gmax_CB,
           gmax_PB = gmax_PB,
           gmax_SB = gmax_SB,
           
           # Half-saturation constant
           K_PB_SR = K_PB_SR,
           K_SB_SO = K_SB_SO,
           K_CB_P = K_CB_P,
           K_PB_P = K_PB_P,
           K_SB_P = K_SB_P,
           
           # Half-inhibition constant
           KH_CB_SR = KH_CB_SR,
           KH_PB_O = KH_PB_O,
           KH_SB_O = KH_SB_O,
           
           # Yields
           y_S_SB = y_S_SB,
           y_S_PB = y_S_PB,
           y_P_CB = y_P_CB,
           y_P_PB = y_P_PB,
           y_P_SB = y_P_SB,
           
           # Production
           p_CB = p_CB,
           
           # Mortality rate
           m_CB = m_CB,
           m_PB = m_PB,
           m_SB = m_SB,
           
           # Diffusivity 
           a_s = a_s,
           a_O = a_O,
           a_P = a_P,
           
           # Background concentration
           S_R_b = S_R_b,
           S_O_b = S_O_b,
           O_b = O_b,
           P_b = P_b,
           
           c = c
)

# Time
times <- seq(from = 0, to = 4000, by = 1)

# SCENARIO I

# Environment
y <- c(N_CB = 5E1,
       N_PB = 1E7,
       N_SB = 1E7,
       
       SO = 300,
       SR = 300,
       O  = 10,
       P  = 10)

init_batch <- ode(y, times, community_model, parms)
init_df_batch <- data.frame(init_batch)
plot_results(init_df = init_df_batch)

# SCENARIO II
y <- c(N_CB = 1E8,
       N_PB = 1E2,
       N_SB = 1E2,
       
       SO = 500,
       SR = 50,
       O  = 300,
       P  = 4)

init_batch <- ode(y, times, community_model, parms)
init_df_batch <- data.frame(init_batch)
plot_results(init_df = init_df_batch)
