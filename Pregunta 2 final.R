#Creador: Daniela Robles
#Analisis de datos de tiempo real 

install.packages('pacman') #pacman sirve para llamar e instalar otros paquetes

library(pacman) #library es para ejecutar

p_load('readr', #readr para llama a las bases de datos
       'dplyr') #para facilitar el manejo de datos

#####################

#llamar base de datos

datos_pcr <- read_csv(file = "https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/examen2")
head(datos_pcr)
#################################

#Obtencion de los genes referencia y de interes

actina <- datos_pcr %>% 
  slice(1)

actina

genes_interes <- datos_pcr %>% 
  slice(-1)

genes_interes

################################
#hasta aqui me da bien los datos
#Para ddct

#extraccion de genes

controles <- datos %>% 
  filter(Condicion=="Control")

head(controles)

#Para sacar promedios

promedio_controles <- controles %>% 
  summarise(Mean_Cx1 = mean(Cx1),
            Mean_Cx2 = mean(Cx2),
            Mean_Cx3 = mean(Cx3),
            Mean_T1 = mean(T1),
            Mean_T2 = mean(T2),
            Mean_T3 = mean(T3)) %>% 
  mutate(Gen="promedio controles") %>% #generar columna de nombres
  select(7,1,2,3,4,5,6)

promedio_controles

#extraer genes de tabla datos

genes <- datos %>% 
  filter(Condicion=="Target") %>% 
  select(-2)
genes            

#sacar el 2^-DCT

DCT <- genes %>% 
  mutate(DCT_C1=2^-(Cx1-promedio_controles$Mean_C1),
         DCT_C2=2^-(Cx2-promedio_controles$Mean_C2),
         DCT_C3=2^-(Cx3-promedio_controles$Mean_C3),
         DCT_T1=2^-(T1-promedio_controles$Mean_T1),
         DCT_T2=2^-(T2-promedio_controles$Mean_T2),
         DCT_T3=2^-(T3-promedio_controles$Mean_T3)) %>% 
  select(-2,-3,-4,-5,-6,-7)
DCT

promedio_genes <- DCT %>% 
  mutate(Mean_DCT_Cx=(DCT_C1+DCT_C2+DCT_C3)/3,
         Mean_DCT_Tx=(DCT_C1+DCT_C2+DCT_C3)/3)
promedio_genes

#######
#Para sacar el 2^-DCT

DCT

promedio_genes <- DCT %>% 
  mutate(Mean_DCT_Cx=(DCT_C1+DCT_C2+DCT_C3)/3,
         Mean_DCT_Tx=(DCT_T1+DCT_T2+DCT_T3)/3)

promedio_genes

###########
DCT <- genes_interes %>% 
  mutate(DCT_Cx1=2^-(Cx1-promedio_controles$Mean_Cx1),
         DCT_Cx2=2^-(Cx2-promedio_controles$Mean_Cx2),
         DCT_Cx3=2^-(Cx3-promedio_controles$Mean_Cx3),
         DCT_T1=2^-(T1-promedio_controles$Mean_T1),
         DCT_T2=2^-(T2-promedio_controles$Mean_T2),
         DCT_T3=2^-(T3-promedio_controles$Mean_T3)) %>% 
  select(-2,-3,-4,-5,-6,-7)
DCT

promedio_genes <- DCT %>% 
  mutate(Mean_DCT_Cx=(DCT_C1+DCT_C2+DCT_C3)/3,
         Mean_DCT_Tx=(DCT_C1+DCT_C2+DCT_C3)/3)
promedio_genes

DDCT <- DCT %>% 
  mutate(DDCT=(Mean_DCT_Tx-Mean_DCT_Cx)) %>% 
  mutate("2^-DDCT"=(2^(-DDCT))) %>% 
  select(1, 11)

DDCT

