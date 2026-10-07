#PROGETTO LABORATORIO DI STATISTICA AMERIGO CIRIELLO
library(dplyr)
library(ggplot2)
library(imputeTS)
library(forecast)
library(psych)
library(dygraphs)
library(seasonal)
library(gridExtra)


# 1. Introduzione

### 1.1 Analisi strutturale

dati <- read.csv("MRTSSM4481USN.csv")
dati$observation_date <- as.Date(dati$observation_date)
dati$MRTSSM4481USN<- as.numeric(dati$MRTSSM4481USN)

glimpse(dati)

head(dati)


# 2. Analisi esplorativa

### 2.1 Analisi esplorativa generale

par(mfrow=c(1,3))
hist(dati$MRTSSM4481USN,col="green",main="Istogramma di vendite retail",freq=FALSE,
     xlab="Vendite In Milioni",ylab="Densità")
lines(density(dati$MRTSSM4481USN))
boxplot(dati$MRTSSM4481USN,main="Boxplot di vendite retail")

hist(dati$MRTSSM4481USN,main="Istogramma con frequenza",freq=TRUE,
     xlab="Vendite In Milioni",ylab="Frequenza")

psych::describe(dati$MRTSSM4481USN)
dev.off()
hist(dati$MRTSSM4481USN, main="Istogramma con Media e Mediana", freq=TRUE,
     xlab="Vendite In Milioni", ylab="Frequenza")
abline(v=mean(dati$MRTSSM4481USN), col="black", lwd=2)
abline(v=median(dati$MRTSSM4481USN), col="red", lwd=2, lty=2)
legend("topright", legend=c("Media", "Mediana"),
       col=c("black", "red"), lwd=2, lty=c(1, 2))

### 2.2 Dati mancanti

Vendite_Retail_USA <- ts(dati$MRTSSM4481USN, start=c(1992, 1), frequency=12)

statsNA(Vendite_Retail_USA)

### 2.3 Time plot

autoplot(Vendite_Retail_USA) +
  theme_linedraw() +
  labs(title="Vendite di abbigliamento al dettaglio negli USA",
       x="Anno", y="Milioni di dollari")+
  scale_x_continuous(breaks=seq(1992, 2026, by = 2))

dyRangeSelector(dygraph(Vendite_Retail_USA, main = "Vendite di abbigliamento al dettaglio negli USA (Grafico interattivo)"))

diff <- tapply(Vendite_Retail_USA, floor(time(Vendite_Retail_USA)), function(x) max(x) - min(x))
anni <- as.numeric(names(diff))
plot(anni, diff, type = "l",
     main = "Progressione differenze tra picchi per intero arco temporale",
     xlab = "Anno", ylab = "Escursione annua (max - min)")

### 2.4 Identificazione preliminare degli Outlier

tsoutliers(Vendite_Retail_USA)

OUTL <- tsoutliers(Vendite_Retail_USA)
df_out <- data.frame(
  tempo  = as.numeric(time(Vendite_Retail_USA))[OUTL$index],
  valore = as.numeric(Vendite_Retail_USA)[OUTL$index]
)

autoplot(Vendite_Retail_USA) +
  geom_point(data = df_out, aes(x = tempo, y = valore),
             color = "red") +
  labs(title = "Vendite abbigliamento USA: potenziali outlier su serie grezza",
       x = "Anno", y = "Milioni di dollari") +
  scale_x_continuous(breaks=seq(1992, 2026, by = 2))+
  theme_linedraw()

### 2.5 Analisi della stagionalità

ggseasonplot(Vendite_Retail_USA, year.labels=TRUE, year.labels.left=TRUE) +
ylab("$ milioni") +
ggtitle("Seasonal plot: Vendite retail abbigliamento USA")+theme_linedraw()

ggsubseriesplot(Vendite_Retail_USA) +
ylab("$ milioni") +
ggtitle("Seasonal subseries plot: Vendite retail abbigliamento USA")+theme_linedraw()

### 2.6 Autocorrelazione (ACF su dati grezzi)

ggAcf(Vendite_Retail_USA)+theme_linedraw()


# 3. Trasformazione, decomposizione e diferenziazione

### 3.1 Aggiustamento rispetto all' inflazione

CPI <- read.csv("CPIAUCNS.csv")
cpi_full <- CPI[CPI$observation_date >= "1992-01-01", ]
cpi_ts <- ts(cpi_full$CPIAUCNS, start = c(1992, 1), frequency = 12)
finestra <- ts.intersect(Vendite_Retail_USA, cpi_ts)
cpi_base <- tail(as.numeric(finestra[, "cpi_ts"]), 1)
Vendite_Reali_CPI <- (finestra[, "Vendite_Retail_USA"] / finestra[, "cpi_ts"]) * cpi_base

statsNA(Vendite_Reali_CPI)

Vendite_Reali_CPI2 <- na_interpolation(Vendite_Reali_CPI)#help na_interpolation
#oppure imputazione sulla serie del CPI
ggplot_na_imputations(Vendite_Reali_CPI, Vendite_Reali_CPI2, theme=ggplot2::theme_classic(), size_points=3)

par(mfrow=c(1,2))
plot(Vendite_Retail_USA,main="Serie grezza")
plot(Vendite_Reali_CPI2,main="Serie corretta CPI")
dev.off()

diff_per_anno <- tapply(Vendite_Reali_CPI2, floor(time(Vendite_Reali_CPI2)), function(x) max(x) - min(x))
anni <- as.numeric(names(diff_per_anno))
plot(anni, diff_per_anno, type = "l",
     main = "Progressione differenze tra picchi per intero arco temporale",
     xlab = "Anno", ylab = "Escursione annua (max - min)")

### 3.2 Trasformata Box-Cox

lambda<-BoxCox.lambda(Vendite_Reali_CPI2,method = "loglik")
lambda

par(mfrow = c(1, 2))
plot(Vendite_Reali_CPI2, main = "Serie reale (nessuna trasf.)", ylab = "Vendite reali")
plot(BoxCox(Vendite_Reali_CPI2, lambda),
     main = paste0("Box-Cox lambda = ", round(lambda, 2)), ylab = "trasformata")
dev.off()
### 3.3 Decomposizione

#### 3.3.1 Decomposizione additiva classica

additive<-decompose(Vendite_Reali_CPI2, type="additive")
autoplot(additive)+ggtitle("Decomposizione additiva per vendite retail")+theme_linedraw()

#### 3.3.2 Decomposizione X11

x11.fit <- seas(Vendite_Reali_CPI2, transform.function="none", x11="")
autoplot(x11.fit) + ggtitle("Decomposizione X11 per vendite retail")+theme_linedraw()

#### 3.3.3 Decomposizione X13 SEATS

x13.seats <- seas(Vendite_Reali_CPI2,transform.function = "none")
autoplot(x13.seats) + ggtitle("Decomposizione X13 SEATS per vendite retail")+theme_linedraw()

#### 3.3.4 Decomposizione STL robusta

stl.fit <- stl(Vendite_Reali_CPI2, s.window="periodic", robust=TRUE) #twindow 21 almeno
autoplot(stl.fit) +
  ggtitle("Decomposizione STL robusta per vendite retail")+theme_linedraw()

autoplot(Vendite_Reali_CPI2, series = "Serie reale", alpha = 0.4) +
  autolayer(trendcycle(additive),series = "Classica",size = 1) +
  autolayer(trendcycle(stl.fit),series = "STL robusta", size = 1) +
  autolayer(trendcycle(x11.fit),series = "X11",         size = 1) +
  autolayer(trendcycle(x13.seats),series = "X13SEATS",size = 1) +
  xlab("Anno") + ylab("Milioni $") +
  ggtitle("Confronto Decomposizioni") +
  guides(colour = guide_legend(title = "Tipologia")) +
  theme_linedraw()


### 3.4 Outlier riesame

R_t <- stl.fit$time.series[, "remainder"]
tsoutliers(R_t)

OUTL_R <- tsoutliers(R_t)
df_out_R <- data.frame(
  tempo  = as.numeric(time(R_t))[OUTL_R$index],
  valore = as.numeric(R_t)[OUTL_R$index]
)

autoplot(R_t) +geom_point(data = df_out_R, aes(x = tempo, y = valore),
             color = "red") +
  labs(title = "Riesame outlier sulla componente residua (STL)",
       x = "Anno", y = "Residuo") +
  scale_x_continuous(breaks = seq(1992, 2026, by = 2)) +
  theme_linedraw()

### 3.5 Stazionarietà e Differenziazione

nsdiffs(Vendite_Reali_CPI2)

Vendite_Reali_DIFFstag <- diff(Vendite_Reali_CPI2, lag = 12)
ndiffs(Vendite_Reali_DIFFstag)

Vendite_CPI2_STAZ <- Vendite_Reali_DIFFstag
autoplot(Vendite_CPI2_STAZ) +
  labs(title = "Serie reale dopo differenza stagionale (lag 12)",
       x = "Anno", y = "Variazioni") +
  theme_linedraw()

ggAcf(Vendite_CPI2_STAZ) +
  theme_linedraw() +
  ggtitle("ACF della serie reale differenziata")


# 4. Confronto con serie analoga

dati3 <- read.csv("MRTSSM4453USN.csv")
Liquor_nom <- ts(as.numeric(dati3$MRTSSM4453USN), start=c(1992,1), frequency=12)
fin_liquor <- ts.intersect(Liquor_nom, cpi_ts)
Liquor_Reali_pre <- (fin_liquor[,"Liquor_nom"] / fin_liquor[,"cpi_ts"]) * cpi_base

statsNA(Liquor_Reali_pre)

Liquor_Reali <- na_interpolation(Liquor_Reali_pre)
ggplot_na_imputations(Liquor_Reali_pre, Liquor_Reali, theme=ggplot2::theme_classic(), size_points=3)

### 4.1 Confronto grafico


g1 <- autoplot(Vendite_Reali_CPI2) +
  labs(title="Abbigliamento", x="Anno", y="Milioni $ reali") +
  theme_linedraw()
g2 <- autoplot(Liquor_Reali) +
  labs(title="Liquori", x="Anno", y="Milioni $ reali") +
  theme_linedraw()
grid.arrange(g1, g2, ncol=1)

### 4.2 Correlazione al netto di trend e stagionalità

rem_abbigliamento <- stl(Vendite_Reali_CPI2, s.window = "periodic")$time.series[, "remainder"]
rem_liquor        <- stl(Liquor_Reali,       s.window = "periodic")$time.series[, "remainder"]
df_rem <- data.frame(
  abbigliamento = as.numeric(rem_abbigliamento),
  liquori       = as.numeric(rem_liquor)
)
ggplot(df_rem, aes(abbigliamento, liquori)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", se=FALSE, colour="red") +
  labs(title="Residui STL: abbigliamento vs liquori",
       x="Residuo abbigliamento", y="Residuo liquori") +
  theme_linedraw()

cor(as.numeric(rem_abbigliamento), as.numeric(rem_liquor))

ccf(as.numeric(rem_abbigliamento), as.numeric(rem_liquor),
    lag.max = 12, main = "CCF: abbigliamento vs liquori (residui)")

#test di lljung box per residui

#confronto tra serie stazionarie senza residui della STL


# 5. Modellazione e previsione

### 5.1 Metodi Benchmark

#### 5.1.1 Average

meanf(Vendite_Reali_CPI2, h=12)

plot1 <- autoplot(Vendite_Reali_CPI2) +
  autolayer(meanf(Vendite_Reali_CPI2, h=12), series="Media", PI=TRUE) +
  theme_linedraw() +
  labs(title = "Previsione Average con Intervalli di Confidenza")
plot2 <- autoplot(Vendite_Reali_CPI2) +
  autolayer(meanf(Vendite_Reali_CPI2, h=12), series="Media", PI=FALSE) +
  theme_linedraw() +
  labs(title = "Previsione Average Puntuale")
grid.arrange(plot1, plot2)

#### 5.1.2 Naive

naive(Vendite_Reali_CPI2, h=12)

plot3<-autoplot(Vendite_Reali_CPI2) +
autolayer(naive(Vendite_Reali_CPI2, h=12), series="Naïve", PI=FALSE)+
  labs(title = "Previsione Naive puntuale")+
  theme_linedraw()

plot4<-autoplot(Vendite_Reali_CPI2) +
autolayer(naive(Vendite_Reali_CPI2, h=12), series="Naïve", PI=TRUE)+
  labs(title = "Previsione Naive con Intervalli di Confidenza")+
  theme_linedraw()

grid.arrange(plot3, plot4)

#### 5.1.3 Seasonal Naive

snaive(Vendite_Reali_CPI2, h=12)

plot5<-autoplot(Vendite_Reali_CPI2) +
autolayer(snaive(Vendite_Reali_CPI2, h=12), 
          series="Seasonal Naive", PI=TRUE)+
  labs(title = "Previsione Seasonal Naive con intervalli di confidenza")+
  theme_linedraw()

plot6<-autoplot(Vendite_Reali_CPI2) +
autolayer(snaive(Vendite_Reali_CPI2, h=12), 
          series="Seasonal Naive ", PI=FALSE)+labs(title = "Previsione Seasonal Naive puntuale")+
  theme_linedraw()

grid.arrange(plot5,plot6)

#### 5.1.4 Drift

rwf(Vendite_Reali_CPI2, h=12, drift=TRUE)

plot7<-autoplot(Vendite_Reali_CPI2) +
autolayer(rwf(Vendite_Reali_CPI2, h=12, drift=TRUE), 
          series="Drift", PI=TRUE)+
  labs(title = "Previsione Drift con intervalli di confidenza")+
  theme_linedraw()

plot8<-autoplot(Vendite_Reali_CPI2) +
autolayer(rwf(Vendite_Reali_CPI2, h=12, drift=TRUE), 
          series="Drift", PI=FALSE)+
  labs(title = "Previsione Drift puntuale")+
  theme_linedraw()

grid.arrange(plot7,plot8)

### 5.2 Previsione con decomposizione STL

stlf(Vendite_Reali_CPI2,method="rwdrift", h=12)
autoplot(stlf(Vendite_Reali_CPI2, h=12))

### 5.3 Regressione per serie storiche

tt <- time(Vendite_Reali_CPI2)
covid <- as.numeric(tt >= 2020 & tt < 2021)

fit.reg <- tslm(Vendite_Reali_CPI2 ~ trend + season + covid)
summary(fit.reg)

autoplot(Vendite_Reali_CPI2, series = "Osservato") +
  autolayer(fitted(fit.reg), series = "Stimato") +
  labs(title = "Regressione: valori osservati vs valori stimati",
       x = "Anno", y = "Milioni di dollari (reali)") +
  guides(colour = guide_legend(title = "Serie")) +
  theme_linedraw()

### 5.4 Diagnostica dei residui + 

checkresiduals(fit.reg)
qqnorm(residuals(fit.reg)); qqline(residuals(fit.reg), col = "red")

df.res <- data.frame(
  fitted = as.numeric(fitted(fit.reg)),
  resid  = as.numeric(residuals(fit.reg)),
  trend  = as.numeric(time(Vendite_Reali_CPI2))
)

p1 <- ggplot(df.res, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "red") +
  labs(title = "Residui vs valori adattati",
       x = "Valori adattati", y = "Residui") +
  theme_linedraw()

p2 <- ggplot(df.res, aes(x = trend, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "red") +
  labs(title = "Residui vs tempo",
       x = "Anno", y = "Residui") +
  theme_linedraw()

gridExtra::grid.arrange(p1, p2, nrow = 1)

### 5.5 Valutazione dell'accuratezza

h <- 12
n <- length(Vendite_Reali_CPI2)
train <- subset(Vendite_Reali_CPI2, end = n-h)
test  <- subset(Vendite_Reali_CPI2, start = n-h+1)

fc.snaive <- snaive(train, h = h)
fc.stl <- stlf(train, s.window = "periodic", robust = TRUE,
               method = "rwdrift", h = h)
tt.train    <- time(train)
covid.train <- as.numeric(tt.train >= 2020 & tt.train < 2021)
fit.reg.tr  <- tslm(train ~ trend + season + covid.train)

fc.reg <- forecast(fit.reg.tr, h = h,
                   newdata = data.frame(covid.train = rep(0, h)))

acc.snaive <- accuracy(fc.snaive, test)
acc.stl    <- accuracy(fc.stl,    test)
acc.reg    <- accuracy(fc.reg,    test)

tabella <- rbind(
  "Seasonal naïve" = acc.snaive["Test set", c("MAE","RMSE","MAPE","MASE")],
  "STL forecast"   = acc.stl["Test set",    c("MAE","RMSE","MAPE","MASE")],
  "Regressione"    = acc.reg["Test set",    c("MAE","RMSE","MAPE","MASE")]
)
tabella

autoplot(Vendite_Reali_CPI2, series = "Osservato") +
  autolayer(fc.snaive, series = "Seasonal naïve", PI = FALSE) +
  autolayer(fc.stl,    series = "STL forecast",   PI = FALSE) +
  autolayer(fc.reg,    series = "Regressione",    PI = FALSE) +
  coord_cartesian(xlim = c(2024, NA)) +
  labs(title = "Confronto delle previsioni sul test set (ultimi 12 mesi)",
       x = "Anno", y = "Milioni di dollari (reali)") +
  guides(colour = guide_legend(title = "Serie")) +
  theme_linedraw()

e.snaive <- tsCV(Vendite_Reali_CPI2, snaive, h = 1, initial = 205)

fstl <- function(y, h) stlf(y, s.window = "periodic", robust = TRUE,
                            method = "rwdrift", h = h)
e.stl <- tsCV(Vendite_Reali_CPI2, fstl, h = 1, initial = 205)
rmse.cv <- c(
  "Seasonal naïve" = sqrt(mean(e.snaive^2, na.rm = TRUE)),
  "STL forecast"   = sqrt(mean(e.stl^2,    na.rm = TRUE))
)
rmse.cv

### 5.6 Previsione finale + Scenario-Based

fc.finale <- forecast(stl.fit, method = "rwdrift", h = 12)
fc.finale

autoplot(fc.finale) +
  labs(title = "Previsione a 12 mesi - STL forecast (dollari reali)",
       x = "Anno", y = "Milioni di dollari (reali)") +
  theme_linedraw()

cpi_ultimo <- cpi_base

prev_reale <- as.numeric(fc.finale$mean)

scenari <- c(basso = 0.02, centrale = 0.03, alto = 0.05
             )
mesi <- 1:12

nominale <- sapply(scenari, function(tasso) {
  cpi_futuro <- cpi_ultimo * (1 + tasso)^(mesi/12)
  (prev_reale * cpi_futuro) / cpi_base
})

nominale <- data.frame(Mese = mesi, Reale = round(prev_reale, 0), round(nominale, 0))
nominale

col_sc <- c("blue","red","green")
matplot(mesi, nominale[, c("basso","centrale","alto")], type = "l", lty = 1,
        col = col_sc, lwd = 2,
        xlab = "Mese previsto", ylab = "Milioni di dollari (nominali)",
        main = "Previsione nominale sotto scenari di inflazione")
legend("topleft", legend = c("2% (basso)","3% (centrale)","5% (alto)"),
       col = col_sc, lwd = 2, bty = "n")

plot(CPI$CPIAUCNS,type='l')
summary(CPI$CPIAUCNS)


